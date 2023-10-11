

#' @title Get results from an SQL query to a RaMP database
#'
#' @description
#'
#' Utility function to execute the provided SQL query and get its results from
#' a RaMP-DB database.
#'
#' @param sql `character(1)` with the SQL query to run.
#'
#' @param db [RaMP()] object representing a RaMP database. By default
#'     (`db = RaMP()`) a connection to the most recent version is established,
#'     which will be downloaded first if it does not already exist in the local
#'     cache.
#'
#' @return The result from the query.
#'
#' @importMethodsFrom DBI dbGetQuery
#'
#' @export
runQuery <- function(sql, db = RaMP()) {
    con <- .dbcon(db)
    on.exit(dbDisconnect(con))
    dbGetQuery(con, sql)
}


#' A user utility method that returns the current DB Type
#' @returns MySQL or SQLite
#' @export
getRampDbType <- function() {

  dbType = 'MySQL'

  if(get("is_sqlite", pkg.globals)) {
    dbType = 'SQLite'
  }

  return(dbType)
}

verifySQLite <- function() {

  message("Checking for existing BiocFileCache entry for the RaMP SQLite Database.")
  bfc <- BiocFileCache::BiocFileCache(cache = BiocFileCache::getBFCOption("CACHE"), ask=F)
  cacheInfo <- BiocFileCache::bfcinfo()
  cacheInfo <- cacheInfo[grepl("RaMP", cacheInfo$rname),]

  if(nrow(cacheInfo) < 1) {
    message("")
    message("RaMP Database is not in file cache. Performing a one-time SQLite file download.")
    url = packageDescription("RaMP")$Config_ramp_db_url
    message("One time retrieval of RaMP Database Cache. This will take about 1 minute to download and unzip.")
    path <- BiocFileCache::bfcadd(bfc, url, fname='exact')
    cid <- names(path)
    R.utils::gunzip(path, remove=F)
    newpath <- gsub(".gz", "", path)
    BiocFileCache::bfcremove(bfc, cid)
    bfcEntry = BiocFileCache::bfcadd(bfc, newpath, fname='exact')
    pkg.globals$sqlite_file_path = bfcEntry
    message("SQLite has been initialized. Using file cache entry:")
    message(bfcEntry)
  } else {
    message("RaMP DB found in BiocFileCache, SQLite File:")
    message(cacheInfo$rpath[1])
    pkg.globals$sqlite_file_path = cacheInfo$rpath[1]
  }
}


setupRdata <- function(db = RaMP()) {

  sql = "select data_key, data_blob from ramp_data_object"

  objs <- RaMP:::runQuery(sql, db)

  dbSummaryData = list()

  for(i in 1:nrow(objs)) {
    varName = objs[i,1]
    blob = objs[i,2]
    blob = blob[[1]]
    obj = memDecompress(from=blob, type = 'gzip', asChar = T)
    data = data.frame(data.table::fread(obj, sep="\t"), row.names = 1)
    dbSummaryData[[varName]] <- data
    # assign(varName, data, envir = .GlobalEnv)
  }

  return(dbSummaryData)
}

