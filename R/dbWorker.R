

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


fetchSQLiteDB <- function() {
  getRaMP

  url = 'https://github.com/ncats/RaMP-DB/commit/054daa67a1a108436d2a22354cc3a7499913c516'

}


getRaMPLibDir <- function() {
  paths <- .libPaths()
  rampPath = NULL
  for(path in paths) {
    rampPath = paste0(path, "/RaMP")
    if(file.exists(rampPath)) {
      print("Hey got the ramp path")
      break
    }
  }
  print(rampPath)
  return(rampPath)
}
