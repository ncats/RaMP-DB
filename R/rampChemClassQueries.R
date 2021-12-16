# queries to retrieve and analyze chemical classes

#' Returns chemical class information comparing a metabolite subset to a metabolite population
#'
#' @param mets a list object of source prepended metaboite ids, representing a metabolite set of interest
#' @param pop an optional list object of source prepended metaboite ids, represenbting a larger list of metabolites from which the mets were selected this list serves as
#' the backround reference population of metabolites for comparision and enrichment. If NULL, the background population is taken as all RaMP DB metabolites.
#' @param includeRaMPids include internal RaMP identifiers (default is "FALSE")
#' @return Returns chemcial class information data including class count tallies and comparisons between metabolites of interest and the metabolite population,
#' metabolite mappings to classes, and query summary report indicating the number of input metabolites that were resolve and listing those metabolite ids
#' that are not found in the database.
#'
#' The returned object (return_obj below) contains three or four main result areas. Use str(return_ob) to see the structure described here.
#'
#' \strong{return_obj$count_summary} contains a dataframe for each category of class annotations (e.g. class, sub_class)
#' This count summary contains:
#' \describe{
#'   \item{category or level name}{The chemical class category, e.g. class, sub_class, super_class}
#'   \item{class name}{The chemical class name, e.g. Organooxygen compounds}
#'   \item{pop_counts and met_counts}{Population and metabolite list counts for each class}
#'   \item{fract_of_pop}{The fraction of the met list's metabolites compared to the population class hits}
#'   \item{fract_within_pop}{The fraction of population metabolites within the class describing how common the class is within the population.}
#'   \item{fract_within_mets}{The fraction of met list metabolites witin the class describing how common the class is within the metabolite list.}
#'}
#'
#' \strong{return_obj$met_classes} metabolite classes for each input metabolite list id
#'
#' \strong{return_obj$pop_classes} metabolite classes for each input population metabolite if a population is provided.
#'
#' \strong{return_obj$query_report} this reports on the query list size, the number of input ids
#' that were found in the database and a list of metabolite ids that were not found in the database.
#' There are two sections, one for the metabolite list and a second when an optional population list is provided.
#'@examples
#'\dontrun{
#' # metabolite list of interest
#' metList = c('hmdb:HMDB0000056',
#'             'hmdb:HMDB0000439',
#'             'hmdb:HMDB0000479',
#'             'hmdb:HMDB0000532',
#'             'hmdb:HMDB0001015',
#'             'hmdb:HMDB0001138',
#'             'hmdb:HMDB0029159',
#'             'hmdb:HMDB0029412',
#'             'hmdb:HMDB0034365',
#'             'hmdb:HMDB0035227',
#'             'hmdb:HMDB0007973',
#'             'hmdb:HMDB0008057',
#'             'hmdb:HMDB0011211')
#'
#' # the background population can be a separate ID list (preferred) or all database entries 
#' # (skip pop parameter).
#' pkg.globals <- setConnectionToRaMP(dbname="ramp2",username="root",conpass="",host = "localhost")
#' metClassResult <- chemicalClassSurvey(mets = mets)
#'
#' # show structure
#' utils::str(metClassResult)
#'
#' # show a count summary, metabolite class mappings and query report
#' metClassResult$count_summary$class
#' metClassResult$met_classes
#' metClassResult$query_report
#'}
#' @export
chemicalClassSurvey <- function(mets, pop = NULL, includeRaMPids = FALSE){
  conn <- connectToRaMP()
  print("Starting Chemical Class Survey")

  if(is.null(pop)) {
      res <- chemicalClassSurveyRampIdsFullPopConn(mets, conn)
  } else {
    res <- chemicalClassSurveyRampIdsConn(mets, pop, conn)
  }
  RMariaDB::dbDisconnect(conn)

  print("Finished Chemical Class Survey")
  if(includeRaMPids){
      return(res)
  }else{
      res$met_classes<-res$met_classes %>% cleanup
      return(res)
  }
}


#' returns chemical class information comparing a metabolite subset to a metabolite population
#'
#' @param mets a list object of source prepended metaboite ids, representing a metabolite set of interest
#' @param pop an optional list object of source prepended metaboite ids, represenbting a larger list of metabolites from which the mets were selected this list serves as
#' the backround reference population of metabolites for comparision and enrichment. If NULL, the background population is taken as all RaMP DB metabolites.
#' @return a data frame containing chemical class enrichment statistics
#'
#'@examples
#'\dontrun{
#' # metabolite list of interest
#' metList = c('hmdb:HMDB0000056',
#'             'hmdb:HMDB0000439',
#'             'hmdb:HMDB0000479',
#'             'hmdb:HMDB0000532',
#'             'hmdb:HMDB0001015',
#'             'hmdb:HMDB0001138',
#'             'hmdb:HMDB0029159',
#'             'hmdb:HMDB0029412',
#'             'hmdb:HMDB0034365',
#'             'hmdb:HMDB0035227',
#'             'hmdb:HMDB0007973',
#'             'hmdb:HMDB0008057',
#'             'hmdb:HMDB0011211')
#'
#' # the background population can be a separate ID list (preferred) or all database entries 
#' # (skip pop parameter).
#' pkg.globals <- setConnectionToRaMP(dbname="ramp2",username="root",conpass="",host = "localhost")
#' metClassResult <- chemicalClassSurvey(mets = mets)
#'
#' enrichedClassStats <- chemicalClassEnrichment(metClassResult)
#'}
#' @export
chemicalClassEnrichment <- function(mets, pop = NULL) {
    print("Starting Chemical Class Enrichment")
    classData <- chemicalClassSurvey(mets = mets, pop = pop,
                                     includeRaMPids = TRUE)
  enrichmentStat <- list()

  totalCountInfo <- getTotalFoundInCategories(classData)

  for (categoryData in classData$count_summary) {
    if(nrow(categoryData) > 0) {
      resultRow <- 1
      classCat <- categoryData[1,'class_level']
      totMetCnt <- totalCountInfo$mets[[classCat]]
      totPopCnt <- as.integer(totalCountInfo$pop[[classCat]])
      contingencyMat <- matrix(nrow=2, ncol=2)
      resultMat <- data.frame(matrix(ncol=7))
      colnames(resultMat) <- c("category", "class_name", "met_hits", "pop_hits",
                               "met_size", "pop_size", "p-value")
      for (i in 1:nrow(categoryData)) {
        if(categoryData[i,'mets_count'] >= 1) {
          contingencyMat[1,1] <- categoryData[i,'mets_count']
          contingencyMat[1,2] <- totMetCnt - contingencyMat[1,1]
          contingencyMat[2,1] <- as.integer(categoryData[i,'pop_count']) - contingencyMat[1,1]
          contingencyMat[2,2] <- totPopCnt - contingencyMat[2,1] - contingencyMat[1,2]
          className <- categoryData[i,'class_name']

          p <- stats::fisher.test(contingencyMat, alternative = "greater")
          p <- p$p.value

          row = list(as.character(classCat), as.character(className), contingencyMat[1,1], contingencyMat[1,1] + contingencyMat[2,1],
                     totMetCnt, totPopCnt, p)

          resultMat[resultRow, ] <- row

          resultRow <- resultRow + 1
        }
      }
      resultMat
      resultMat <- bhCorrect(resultMat)
      enrichmentStat[[as.character(classCat)]] <- resultMat
    }
  }
  print("Finished Chemical Class Enrichment")
  return(enrichmentStat)
}

