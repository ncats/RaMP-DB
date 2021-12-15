#' Retrieve RaMP version
#' @return current ramp databse version
#' @examples
#' \dontrun{
#' pkg.globals <- setConnectionToRaMP(dbname="ramp2",username="root",conpass="",host = "localhost")
#' getCurrentRaMPVersion()
#' }
#' @export
getCurrentRaMPVersion<-function(){
  con<-connectToRaMP()
  query1<-"select ramp_version from db_version where load_timestamp order by load_timestamp desc limit 1"
  results<-RMariaDB::dbGetQuery(con,query1)
  RMariaDB::dbDisconnect(con)  
  return(results)
}

#' Retrieve versions of each database within the current version of RaMP
#' @return database source version info
#' @examples
#' \dontrun{
#' pkg.globals <- setConnectionToRaMP(dbname="ramp2",username="root",conpass="",host = "localhost")
#' getCurrentRaMPDBVersions()
#' }
#' @export
getCurrentRaMPDBVersions<-function(){
  con <- connectToRaMP()
  query1<- "select * from version_info where status = 'current'"
  results<- RMariaDB::dbGetQuery(con,query1)
  RMariaDB::dbDisconnect(con)
  return(results)
}


#' Retrieve counts of entitites (e.g. Metabolites, Pathways, Metabolite-Pathway associations, etc.) for RaMP source databases
#' @return database sources and entity counts associated with each data source
#' @examples
#' \dontrun{
#' pkg.globals <- setConnectionToRaMP(dbname="ramp2",username="root",conpass="",host = "localhost")
#' getEntityCountsFromSourceDBs()
#' }
#' @export
getEntityCountsFromSourceDBs<-function(){
  entity_source_name <- entity_count <- c()
  con<-connectToRaMP()
  query1<-"select * from entity_status_info"
  results<-RMariaDB::dbGetQuery(con,query1)
  RMariaDB::dbDisconnect(con)
  results<-results[,-2]
  results<-results %>% tidyr::spread(unique(entity_source_name),entity_count,fill = "-")
  return(results)
}


