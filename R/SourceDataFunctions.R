#' Retrieve RaMP version
#' @param justVersion boolean value indicating if the method should just return the version id (default, justVersion = T),
#' or a table that includes db_version_id, load_timestamp (update time/date), version_notes, and the db_sql_url (a url for mysql schema download)
#' @param db a RaMP database object
#' @return current ramp databse version
#' @examples
#' \dontrun{
#' rampDB <- RaMP()
#' getCurrentRaMPVersion(db = rampDB)
#' }
#' @export
getCurrentRaMPVersion<-function(justVersion=T, db = RaMP()){
  results <- db@api$getRaMPVersion(justVersion = justVersion)
  return(results)
}

#' Retrieve versions of each database within the current version of RaMP
#' @param db a RaMP database object
#' @return database source version info
#' @examples
#' \dontrun{
#' getCurrentRaMPSourceDBVersions()
#' }
#' @export
getCurrentRaMPSourceDBVersions<-function(db = RaMP()){
  results<- db@api$getCurrentSourceVersion()
  return(results)
}

#' Retrieve counts of entitites (e.g. Metabolites, Pathways, Metabolite-Pathway associations, etc.) for RaMP source databases
#' @param db a RaMP database object
#' @return database sources and entity counts associated with each data source
#' @examples
#' \dontrun{
#' getEntityCountsFromSourceDBs()
#' }
#' @export
getEntityCountsFromSourceDBs<-function(db = RaMP()){
  entity_source_name <- entity_count <- c()
  results<- db@api$getEntityCountsFromSources()
  results<-results[,-2]
  results<-results %>% tidyr::spread(unique(entity_source_name),entity_count)
  results[is.na(results)]=0
  results<- with(results,{transform(results, HMDB = as.numeric(HMDB))})
  return(results)
}

#' Retrieve RaMP Analyte Source Intersections, these indicate the level of analyte overlaps between our sources
#' @param analyteType returns analyte overlaps for 'metabolites' or 'genes'
#' @param format can be one of either 'json', 'upsetR_expression'
#' @param scope value in c('global', 'mapped-to-pathway'), indicates all metabolite stats should be returned, or just those associated with pathways.
#' @param db a RaMP database object
#' @return current analyte overlaps counts between current data sources, for the specified analyteType
#' @examples
#' \dontrun{
#' jsonResult <- getRaMPAnalyteIntersections(analyteType='genes', format='json')
#' }
#' @export
getRaMPAnalyteIntersections<-function( analyteType='metabolites', format='json', scope='mapped-to-pathway', db = RaMP()){
  if (!(analyteType %in% c('metabolites', 'genes'))) {
    warning("The analyteType must be one of c('metabolites','genes')")
    #return an empty dataframe
    return(data.frame())
  }

  results <- db@api$getAnalyteIntersects(analyteType = analyteType, scope = scope)

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
    #return an empty dataframe
    return(data.frame())
  }
  return(results)
}



#' Retrieve list of pathway names
#' @return vector of unique pathway names (alphabetically ordered)
#' @param db a RaMP database object
#' @examples
#' \dontrun{
#' getPathwayNameList(db = rampDB)
#' }
#' @export
getPathwayNameList <- function(db = RaMP()){
  pathwayNames <- db@api$getPathwayNames()
  return(pathwayNames)
}
