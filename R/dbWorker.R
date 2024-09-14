

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
runQuery <- function(
    sql, db = RaMP()) {
    con <- .dbcon(x = db)
    on.exit(dbDisconnect(conn = con))
    dbGetQuery(conn = con, statement = sql)
}

setupRdataCache <- function(db = RaMP()) {

  sql = "select data_key, data_blob from ramp_data_object"

  objs <- runQuery(sql = sql, db = db)

  dbSummaryData = list()

  for(i in 1:nrow(objs)) {
    varName = objs[i,1]
    blob = objs[i,2]
    blob = blob[[1]]
    obj = memDecompress(from=blob, type = 'gzip', asChar = T)
    data = data.frame(data.table::fread(input = obj, sep="\t"), row.names = 1)
    dbSummaryData[[varName]] <- data
  }

  return(dbSummaryData)
}
