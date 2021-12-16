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
getCurrentRaMPSourceDBVersions<-function(){
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

#' Retrieve RaMP Analyte Source Intersections, these indicate the level of analyte overlaps between our sources
#' @param analyteType returns analyte overlaps for 'metabolites' or 'genes'
#' @param format can be one of either 'json', 'upsetR_expression'
#' @return current analyte overlaps counts between current data sources, for the specified analyteType
#' @examples
#' \dontrun{
#' pkg.globals <- setConnectionToRaMP(dbname="ramp2",username="root",conpass="",host = "localhost")
#' jsonResult <- getRaMPAnalyteIntersections(analyteType='genes', format='json')
#' }
#' @export
getRaMPAnalyteIntersections<-function(analyteType='metabolites', format='json'){
  if(analyteType == 'metabolites') {
    query<-"select met_intersects_json from db_version where load_timestamp order by load_timestamp desc limit 1"
  } else if (analyteType == 'genes') {
    query<-"select gene_intersects_json from db_version where load_timestamp order by load_timestamp desc limit 1"
  } else {
    warning("The analyteType must be one of c('metabolites','genes')")
  }
  con<-connectToRaMP()
  results<-RMariaDB::dbGetQuery(con,query)
  RMariaDB::dbDisconnect(con)

  if(format == 'json') {
    if(nrow(results)>0) {
      results <- results[1,1]
    } else {
      results = ""
    }
  } else if(format == 'upsetR_expression') {
    #convert json to list
    if(nrow(results)>0) {
      resultList <- jsonlite::fromJSON(results[1,1])
      setList <- c()
      sizeList <- c()
      for(set in resultList$sets) {
        setList = c(setList, paste(set,collapse="&"))
      }
      for(size in resultList$size) {
        sizeList = c(sizeList, size)
      }

      print(length(sizeList))
      print(length(setList))

      df <- data.frame(matrix(sizeList,nrow=1))
      colnames(df) <- setList
      results = df
    }
  } else {
    warning("The format must be one of c('json','upsetR_expression')")
  }
  return(results)
}



