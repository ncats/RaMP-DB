

#' Utility wrapper method for queries handling mysql and sqlite
#' @param sql sql query to run
#' @returns returns the query result, Note: The returned data type is dependent on query
#' @export
runQuery <- function(sql) {
  conn = RaMP::connectToRaMP()
  if(get("is_sqlite",pkg.globals)) {
    rs <- RSQLite::dbGetQuery(conn, sql)
    RSQLite::dbDisconnect(conn)
  } else {
    rs <- RMariaDB::dbGetQuery(conn, sql)
    RMariaDB::dbDisconnect(conn)
  }
  return(rs)
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


setupRdata <- function() {
  if(!exists("kegg_gene")) {

    sql = "select data_key, data_blob from ramp_data_object"

    objs <- RaMP:::runQuery(sql)

    for(i in 1:nrow(objs)) {
      varName = objs[i,1]
      blob = objs[i,2]
      blob = blob[[1]]
      obj = memDecompress(from=blob, type = 'gzip', asChar = T)
      data = data.frame(data.table::fread(obj, sep="\t"), row.names=1)
      assign(varName, data, envir = .GlobalEnv)
    }
  }
}

