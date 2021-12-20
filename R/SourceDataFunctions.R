#' getCurrentDbVersion
#' @return current ramp database version
getCurrentDbVersion<-function(){
  con<-connectToRaMP()
  query1<-"select ramp_version from db_version where load_timestamp order by load_timestamp desc limit 1"
  results<-RMariaDB::dbGetQuery(con,query1)
  return(results)
}

#' getCurrentVersionInfo
#' @return database source version info
getCurrentVersionInfo<-function(){
  con <- connectToRaMP()
  query1<- "select * from version_info where status = 'current'"
  results<- RMariaDB::dbGetQuery(con,query1)
  return(results)
}


#' getEntityCountsFromSources
#' @return database sources and entity counts associated with each data source
getEntityCountsFromSources<-function(){
  con<-connectToRaMP()
  query1<-"select * from entity_status_info"
  results<-RMariaDB::dbGetQuery(con,query1)
  results<-results[,-2]
  results<-results %>% spread(unique(entity_source_name),entity_count)
  results[is.na(results)]=0
  results<- transform(results, HMDB = as.numeric(HMDB))
  return(results)
}



