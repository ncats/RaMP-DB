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
        port = "integer"
    ),
    prototype = prototype(
        driver = NULL,
        dbname = character(),
        username = character(),
        conpass = character(),
        host = character(),
        port = integer()
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
    con <- dbConnect(x@driver, dbname = .dbname(x), user = .username(x),
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
#'
#' @rdname RaMP
setMethod("show", "RaMP", function(object) {
    if (is.null(object@driver))
        cat("Empty RaMP object")
    con <- .dbcon(object)
    on.exit(dbDisconnect(con))
    cat(class(object), "\n")
    ## Maybe get some additional information from the database with e.g.
    ## number of analytes or versions and list them.
})

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
RaMP <- function(version = character()) {
    db_local <- listRaMPVersions(local = TRUE)
    if (!length(version)) {
        ## Get most recent remote version
        db_remote <- listRaMPVersions(local = FALSE)
        if (!length(db_remote))
            stop("Error getting available remote versions")
        version <- db_remote[length(db_remote)]
    }
    if (!version %in% db_local) {
        ## Only check for remote versions if database not already cached
        db_remote <- listRaMPVersions(local = FALSE)
        if (!version %in% db_remote)
            stop("RaMP version '", version,"' not available. Use ",
                 "'listRaMPVersions()' to list available versions.")
    }
    db <- .RaMP(SQLite(), dbname = .get_ramp_db(version))
    con <- .dbcon(db)
    ## Maybe retrieve additional tables or information from the database
    ## and cache/store that in a slot within the RaMP object?
    on.exit(dbDisconnect(con))
    .valid_ramp_database(con, error = TRUE)
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
    new("RaMP", driver = driver, dbname = dbname, username = username,
        conpass = conpass, host = host, port = port)
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
#' @export
listRaMPVersions <- function(local = FALSE) {
    if (local) {
        bfc <- BiocFileCache(cache = getBFCOption("CACHE"), ask = FALSE)
        ci <- bfcinfo(bfc)
        ramps <- ci$rname[grepl("RaMP", ci$rname)]
        sort(unname(vapply(ramps, .version_from_db_file, character(1))))
    } else {
        ## Get all available (and supported versions) from e.g. figshare,
        ## zenodo, github?
        ## Placeholder until we figure out how to get available online
        ## version
        "2.3.0"
    }
}

#' @importFrom BiocFileCache BiocFileCache getBFCOption bfcinfo bfcadd bfcremove
#'
#' @description
#' 
#' Check if a RaMP-DB for the specific version is available and download it
#' otherwise.
#' 
#' @noRd
.get_ramp_db <- function(version) {
    bfc <- BiocFileCache(cache = getBFCOption("CACHE"), ask = FALSE)
    cacheInfo <- bfcinfo()
    cacheInfo <- cacheInfo[grepl(paste0("RaMP_SQLite_v", version, ".sqlite"),
                                 cacheInfo$rname), ]
    if (!nrow(cacheInfo)) {
        message("Downloading RaMP-DB version ", version)
        db_url <- paste0(
            "https://github.com/ncats/RaMP-DB/raw/sqlite/db/RaMP_SQLite_v",
            version, ".sqlite.gz")
        path <- bfcadd(bfc, db_url, fname = "exact")
        dbf <- sub(".gz", "", path, fixed = TRUE)
        if (file.exists(dbf))
            file.remove(dbf)
        R.utils::gunzip(path, remove = FALSE)
        bfcremove(bfc, names(path))
        db_file <- bfcadd(bfc, dbf, fname = "exact")
    } else {
        message("Loading RaMP-DB version ", version, " from cache.")
        db_file <- cacheInfo$rname[1L]
    }
    db_file
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
