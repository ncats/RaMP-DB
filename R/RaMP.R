#' @importFrom methods setClassUnion
#'
#' @importClassesFrom DBI DBIDriver
#'
#' @noRd
setClassUnion("DBIDriverOrNULL", c("DBIDriver", "NULL"))

#' Maybe have an additional slot of type `list` with additional information/
#' metadata retrieved from the database?
#'
#' @noRd
setClass(
    "RaMP",
    slots = c(
        driver = "DBIDriverOrNULL",
        dbname = "character",
        username = "character",
        conpass = "character",
        host = "character",
        port = "integer",
        dbSummaryObjCache = "list",
        versionSupport = "environment"
    ),
    prototype = prototype(
        driver = NULL,
        dbname = character(),
        username = character(),
        conpass = character(),
        host = character(),
        port = integer(),
        dbSummaryObjCache = list(),
        versionSupport = new.env()
    ))

#' Helper function to return the connection to the database, defined by the
#' internal settings of the RaMP object.
#'
#' @param x `RaMP` object.
#'
#' @return a DB connection object.
#'
#' @noRd
.dbcon <- function(x) {
  con <- dbConnect(drv = x@driver, dbname = .dbname(x), username = .username(x),
                            password = .conpass(x), host = .host(x), port = .port(x))
}

.dbname <- function(x) {
    if (length(x@dbname)) x@dbname
    else NULL
}

.username <- function(x) {
    if (length(x@username)) x@username
    else NULL
}

.conpass <- function(x) {
    if (length(x@conpass)) x@conpass
    else NULL
}

.host <- function(x) {
    if (length(x@host)) x@host
    else NULL
}

.port <- function(x) {
    if (length(x@port)) x@port
    else NULL
}

#' Helper function to check if the connection is/will be to a
#' SQLite database
#'
#' @noRd
.is_sqlite <- function(x) {
    inherits(x@driver, "SQLiteDriver")
}

#' @importMethodsFrom methods show
#'
#' @importFrom DBI dbDisconnect
#'
#' @exportMethod show
#' @param object a RaMP object
#' @rdname RaMP
setMethod("show", "RaMP", function(object) {
    if (is.null(object@driver))
        cat("Empty RaMP object")
    con <- .dbcon(x = object)
    on.exit(dbDisconnect(con))
    cat(class(object), "\n")
    ## Maybe get some additional information from the database with e.g.
    ## number of analytes or versions and list them.
})

setMethod(
  "initialize",
  "RaMP",
  function(.Object, ...) {
    .Object@versionSupport <- new.env()
    callNextMethod()
  }
)


#' @title Connection to a RaMP database
#'
#' @aliases show
#'
#' @description
#'
#' Connections to a *RaMP* database can be established and managed with the
#' `RaMP` function. The returned `RaMP` object provides the reference to the
#' database and it can be passed to the various functions to query that
#' specific database. RaMP databases are provided as self-contained SQLite
#' databases that are automatically downloaded and locally cached with the
#' `RaMP` function. The caching mechanism prevents repeated downloads of the
#' same database version.
#'
#' - `RaMP`: eventually download and connect to a RaMP database. Parameter
#'   `version` allows to specify the RaMP release version to which a connection
#'   should be established. If the specified version is not available locally,
#'   it will be downloaded and cached. Use `listRaMPVersions()` to list
#'   available local or remote databases. Alternatively, the connection to a
#'   RaMP database can be directly provided through parameter `dbcon`.
#'
#' - `listRaMPVersions`: list available local or remote RaMP database releases.
#'
#' @param version `character(1)` specifying the RaMP version to load. By
#'     default (`version = character()`), the most recent release will be
#'     used.
#'
#' @param local `logical(1)` for `listRaMPVersion`: whether remote
#'     (`local = FALSE`, default) or locally (`local = TRUE`) available RaMP
#'     versions should be listed.
#'
#' @name RaMP
#'
#' @importFrom methods new
#'
#' @importFrom RSQLite SQLite
#'
#' @importFrom DBI dbConnect
#'
#' @export
RaMP <- function(version = character(), branch = "main") {
    db_local <- listRaMPVersions(local = TRUE, branch = branch)
    if (!length(version)) {
        ## Get most recent remote version
        db_remote <- listRaMPVersions(local = FALSE, branch = branch)
        if (!length(db_remote))
            stop("Error getting available remote versions")

        # remote versions are returned in decreasing version order, take first.
        version <- db_remote[1]
    }
    if (!version %in% db_local) {
        ## Only check for remote versions if database not already cached
        db_remote <- listRaMPVersions(local = FALSE, branch = branch)
        if (!version %in% db_remote) {
            print(paste0("RaMP version '", version,"' not available. Use ",
                 "'listAvailableRaMPDbVersions()' to list available versions."))
            print(paste0("Checking for version '", version, "' on remote server."))
              avail <- .is_version_in_remote_lfs(version = version, branch = branch)
            if(!avail) {
              print("The spcified RaMP Database version is not available in local file cache OR in remote repository.")
              print("")
              print("Available versions [ running listAvailableRaMPDbVersions() ]:")
              listAvailableRaMPDbVersions(branch = branch)
              return(NULL)
            } else {
              print(paste0("Retrieving RaMP SQLite DB version '", version, "' from remote server."))
              .get_ramp_db(version = version, branch = branch)
            }
        }
    }
    db <- .RaMP(driver = SQLite(), dbname = .get_ramp_db(version = version, branch = branch))
    con <- .dbcon(x = db)

    # add the cache of summary data objects for enrichment
    db@dbSummaryObjCache <- setupRdataCache(db = db)

    on.exit(dbDisconnect(con))
    .valid_ramp_database(con = con, error = TRUE)
    db
}

#' Internal constructor - to also support MySQL/MariaDB connections.
#'
#' @importFrom RSQLite SQLite
#'
#' @noRd
.RaMP <- function(driver = SQLite(), dbname = character(),
                  username = character(), conpass = character(),
                  host = character(), port = integer()) {

    rampObj <- new("RaMP", driver = driver, dbname = dbname, username = username,
        conpass = conpass, host = host, port = port, dbSummaryObjCache = list())

    # creates the cache of R data objects
    rampObj@dbSummaryObjCache <- setupRdataCache(db = rampObj)

    setupVersionSupport(rampObj)

    return(rampObj)
}

#' simple validator function checking for validity of a RaMP database.
#'
#' @importFrom DBI dbListTables
#'
#' @noRd
.valid_ramp_database <- function(con, error = FALSE) {
    .required_tables <- c("db_version", "version_info", "analyte")
    msg <- character()
    if (!inherits(con, "DBIConnection"))
        msg <- "'con' is not a valid database connection."
    else {
        tbls <- dbListTables(con)
        if (!all(.required_tables %in% tbls))
            msg <- paste0("Database lacks required database tables. ",
                          "Is 'con' a connection to a RaMP database?")
        ## Possibly other validity tests
    }
    if (error && length(msg))
        stop(msg)
    else msg
}

#' @rdname RaMP
#'
listRaMPVersions <- function(local = FALSE, branch = "main") {
    if (local) {
        bfc <- BiocFileCache(cache = getBFCOption("CACHE"), ask = FALSE)
        ci <- bfcinfo(bfc)
        ramps <- ci$rname[grepl("RaMP", ci$rname)]
        sort(unname(vapply(ramps, .version_from_db_file, character(1))))
    } else {
       .get_remote_db_version_list(branch = branch)
    }
}

#' @importFrom BiocFileCache BiocFileCache getBFCOption bfcinfo bfcadd bfcremove bfcnew
#'
#' @description
#'
#' Check if a RaMP-DB for the specific version is available and download it
#' otherwise.
#'
#' @noRd
.get_ramp_db <- function(version, branch = "main") {
    bfc <- BiocFileCache(cache = getBFCOption("CACHE"), ask = FALSE)
    cacheInfo <- bfcinfo()
    cacheInfo <- cacheInfo[grepl(paste0("RaMP_SQLite_v", version, ".sqlite"),
                                 cacheInfo$rname), ]
    if (!nrow(cacheInfo)) {
        if (branch == "main") {
            message("Downloading RaMP-DB version ", version)
        } else {
            message("Downloading RaMP-DB version ", version, " from ", branch, " branch")
        }
        db_url <- paste0(
            "https://github.com/ncats/RaMP-DB/raw/", branch, "/db/RaMP_SQLite_v",
            version, ".sqlite.gz")
        path <- bfcadd(bfc, db_url, fname = "exact", archiveMethod='unzip')
        dbf <- sub(".gz", "", path, fixed = TRUE)
        if (file.exists(dbf))
            file.remove(dbf)
        R.utils::gunzip(path, remove = TRUE)
        bfcremove(bfc, names(path))
        db_file <- BiocFileCache::bfcnew(bfc, dbf)
    } else {
        message("Loading RaMP-DB version ", version, " from cache.")
        dbf <- cacheInfo$rname[1L]
    }
    dbf
}

#' Extract the version from a RaMP SQLite file name. Expected format:
#' RaMP_SQLite_v2.3.0.sqlite.gz.
#'
#' @noRd
.version_from_db_file <- function(x) {
    v <- sub(".*SQLite_v(.*?)\\.sqlite.*", "\\1", x)
    if (v == x)
        v <- ""
    v
}

#' @importFrom httr HEAD
#'
#' @description
#'
#' Tests if a version exists in remote github LFS
#'
#' @noRd
.is_version_in_remote_lfs <- function(version, branch="main") {
  db_url <- paste0(
    "https://github.com/ncats/RaMP-DB/raw/", branch, "/db/RaMP_SQLite_v",
    version, ".sqlite.gz")
  resp <- httr::HEAD(db_url)
  return(resp$status_code == 200)
}

#'
#' @description returns the lists of RaMP db versions available
#'
#' @noRd
.get_local_db_version_list <- function() {
  localVersions = c()
  bfc <- BiocFileCache(cache = getBFCOption("CACHE"), ask = FALSE)
  cacheInfo <- bfcinfo()
  suppressWarnings( {
    cacheInfo <- cacheInfo[grepl("RaMP_SQLite_v", cacheInfo$rname), ]$rname
    if(length(cacheInfo) > 0) {
      localVersions <- unlist(lapply(cacheInfo, FUN=.version_from_db_file))
    }
  })

  localVersions <- unique(localVersions)
  localVersions <- sort(localVersions, decreasing = T)

  return(localVersions)
}

#' @importFrom httr HEAD GET
#'
#' @description returns the list of RaMP db versions available
#'
#' @noRd
.get_remote_db_version_list <- function(branch = "main") {

  # a bit of a hack to parse html... and a regexpr could be cleaner...
  remoteURL = paste0("https://github.com/ncats/RaMP-DB/raw/", branch, "/db/")
  filenames = httr::GET(remoteURL, ftp.use.epsv = FALSE, dirlistonly = TRUE)
  filenames <- rawToChar(filenames$content)
  filelocs <- unlist(gregexpr(pattern = 'SQLite_v', filenames))
  filelocs <- filelocs + 8
  fileEnds <- filelocs + 18
  versions <- c()
  for(i in 1:length(filelocs)) {
    versions <- c(versions, substr(filenames, filelocs[i], fileEnds[i]))
  }
  versions <- substr(versions, 1, unlist(gregexpr(".sqlite", versions))-1)

  remoteVersions <- sort(versions, decreasing=T)
  remoteVersions <- unique(remoteVersions)

  return(remoteVersions)
}

#'
#' Lists local and remotely available RaMP SQLite DB versions and prompts with
#' message to download a new version if one exists.
#'
#' @export
listAvailableRaMPDbVersions <- function(branch = "main") {

  localVersions <- .get_local_db_version_list()
  remoteVersions <- .get_remote_db_version_list(branch = branch)

  newVersions <- setdiff(remoteVersions, localVersions)
  haveLocalVersions <- (length(localVersions)>0)

  print("Locally available versions of RaMP SQLite DB, currently on your computer:")
  if(haveLocalVersions) {
    print(localVersions)
  } else {
    print("No local versions of the RaMP Database were found.")
    print("Please use the command 'db <- RaMP()' to download the latest version into local file cache.")
    print("Alternatively you can use the command db <- RaMP(version = <remote_version_number>) using one of the versions listed below.")
  }

  if (branch == "main") {
    print("Available remote RaMP SQLite DB versions for download:")
  } else {
    print(paste0("Available remote RaMP SQLite DB versions (", branch, " branch) for download:"))
  }
  print(remoteVersions)

  if(length(newVersions) > 0 && haveLocalVersions) {
    print("The following RaMP Database versions are available for download:")
    print(newVersions)
    print("Use the command db <- RaMP(<new_version_number>) to download the specified version.")
  }

}


#'
#' Remove a RaMP database version from local file cache
#'
#' @param version a ramp version as a string argument, e.g. '2.4.0'. The version parameter can be set to 'all'.
#'
#' @export
#'
removeLocalRampDB <- function(version = 'none') {

  if(version == 'none') {
    message("Please specify a RaMP database version to remove from local cache, or specify version as 'all' to clear file cache.")
  } else {

    localVersions <- .get_local_db_version_list()
    cacheInfo <- BiocFileCache::bfcinfo()

    if(version == 'all') {

      hits <- grepl("RaMP", cacheInfo$rname)

      # if all entries in the cache contain 'RaMP' just remove the cache
      if(sum(hits) == nrow(cacheInfo)) {
        BiocFileCache::removebfc(ask=FALSE)
      } else {
        rids <- cacheInfo$rid[hits]
        BiocFileCache::bfcremove(rids=rids)
      }

      message("The RaMP database file cache has been cleared.")

    } else {

      if(!(version %in% localVersions)) {
        message("The specified version ramp DB version (", version, ") is not in the local file cache.")
        if(length(localVersions) == 0) {
          message("The RaMP database file cache is empty.")
        } else {
          message("Local RaMP database versions: ", paste(localVersions, collapse=", "))
        }
      } else {

        hits <- grepl(version, cacheInfo$rname)
        ids <- unlist(cacheInfo$rid[hits])
        if(length(ids) == 1) {
          BiocFileCache::bfcremove(rids = ids[1])
          message("Local RaMP database version removed: ", version)
        }
      }
    }
  }
}
