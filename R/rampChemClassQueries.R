# queries to retrieve and analyze chemical classes

#' Returns chemical class information comparing a metabolite subset to a metabolite population
#'
#' @param mets a list object of source prepended metaboite ids, representing a metabolite set of interest
#' @param pop an optional list object of source prepended metaboite ids, represenbting a lareger list of metabolites from which the mets were selected this list serves as
#' the backround reference population of metabolites for comparision and enrichment. If NULL, the background population is taken as all RaMP DB metabolites.
#' @param conpass the ramp database password
#' @param dbname the ramp database name
#' @param host the ramp database host name
#' @param username the ramp database user name
#' @return Returns chemcial class information data including class count tallies and comparisons between metabolites of interest and the metabolite population,
#' metabolite mappings to classes, and query summary report indicating the number of input metabolites that were resolve and listing those metabolite ids
#' that are not found in the database.
#'
#' The returned object (return_obj below) contains three or four main result areas. Use str(return_ob) to see the structure described here.
#'
#' \strong{return_obj$count_summary} contains a dataframe for each category of class annotations (e.g. class, sub_class)
#' This count summary contains:
#' \describe
#'   \item{category or level name}{The chemical class category, e.g. class, sub_class, super_class}
#'   \item{class name}{The chemical class name, e.g. Organooxygen compounds}
#'   \item{pop_counts and met_counts}{Population and metabolite list counts for each class}
#'   \item{fract_of_pop}{The fraction of the met list's metabolites compared to the population class hits}
#'   \item{fract_within_pop}{The fraction of population metabolites within the class describing how common the class is within the population.}
#'   \item{fract_within_mets}{The fraction of met list metabolites witin the class describing how common the class is within the metabolite list.}
#'
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
#' # the background population can be a separate ID list (preferred) or all database entries (skip pop parameter).
#' metClassResult <- chemicalClassSurvey(mets = mets, conpass, dbname, host, username)
#'
#' # show structure
#' str(metClassResult)
#'
#' # show a count summary, metabolite class mappings and query report
#' metClassResult$count_summary$class
#' metClassResult$met_classes
#' metClassResult$query_report
#'}
#' @export
chemicalClassSurvey <- function(mets, pop = NULL,
                                conpass,
                                dbname,
                                host,
                                username) {

  conn <- connectToRaMP(conpass=conpass, dbname = dbname, host=host, username=username)

  print("Starting Chemical Class Survey")

  if(is.null(pop)) {
    res <- chemicalClassSurveyRampIdsFullPopConn(mets, conn)
  } else {
    res <- chemicalClassSurveyRampIdsConn(mets, pop, conn)
  }
  RMariaDB::dbDisconnect(conn)

  print("Finished Chemical Class Survey")
  return(res)
}


#' returns chemical class information comparing a metabolite subset to a metabolite population
#'
#' @param classData a chemical class result object from chemicalClassSurvey
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
#' # the background population can be a separate ID list (preferred) or all database entries (skip pop parameter).
#' metClassResult <- chemicalClassSurvey(mets = mets, conpass, dbname, host, username)
#'
#' enrichedClassStats <- chemicalClassEnrichment(metClassResult)
#'}
#' @export
chemicalClassEnrichment <- function(classData) {
  print("Starting Chemical Class Enrichement")
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

          p <- fisher.test(contingencyMat, alternative = "greater")
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


###########
#
# Supporting functions
#
###########


# check for id prefixes
checkIdPrefixes <- function(idList) {
  idCount <- length(idList)
  prefixCount <- 0
  for(id in idList) {
    if(grepl(":",id, fixed = TRUE)) {
      prefixCount <- prefixCount + 1
    }
  }
  if(prefixCount/idCount < 0.9) {
    warn <- paste("RaMP expects ids to be prefixed with the source database." + (idCount-prefixCount) + " of " + idCount + " ids lack prefixes.\n", sep="")
    warnObj <- Warning(warn, call=TRUE, immediate=TRUE)
    print("Common metabolite prefixes: CAS:, chebi:, chemspider:, hmdb:, kegg:, LIPIDMAPS:, pubchem:")
    print("Examples: kegg:C02712, hmdb:HMDB04824, CAS:2566-39-4. The input list may contain a variety of id types.")
  }
}


# runs and tallies chemical class information *when there is a user defined population
chemicalClassSurveyRampIdsConn <- function(mets, pop, conn) {

  mets <- unique(mets)

  checkIdPrefixes(mets)

  pop <- unique(pop)

  checkIdPrefixes(pop)

  result <- list()

  # first handle metabolites of interest
  metStr <- paste(mets, collapse = "','")
  metStr <- paste("'" ,metStr, "'", sep = "")

  sql <- paste("select distinct a.ramp_id, b.sourceId, a.class_level_name, a.class_name, a.source from metabolite_class a, source b
          where b.rampId = a.ramp_id and b.sourceId in (",metStr,")")

  metsData <- RMariaDB::dbGetQuery(conn, sql)

  # need to filter for our specific source ids
  metsData <- subset(metsData, sourceId %in% mets)

  # get query summary
  metQueryReport <- queryReport(mets, metsData$sourceId)

  metsCountData <- data.frame(table(metsData$class_level_name,metsData$class_name))
  colnames(metsCountData) <- c("class_level", "class_name", "freq")
  metsCountData <- metsCountData[metsCountData$freq != 0,]
  metsCountData <- metsCountData[order(-metsCountData$freq),]

  print("...finished metabolite list query...")

  # Population info
  popStr <- paste(pop, collapse = "','")
  popStr <- paste("'" ,popStr, "'", sep = "")

  sql <- paste("select distinct a.ramp_id, b.sourceId, a.class_level_name, a.class_name, a.source from metabolite_class a, source b
          where b.rampId = a.ramp_id and b.sourceId in (",popStr,")")

  popData <- RMariaDB::dbGetQuery(conn, sql)

  #need to filter for our source ids
  popData <- subset(popData, sourceId %in% pop)

  # get query summary
  popQueryReport <- queryReport(pop, popData$sourceId)

  popCountData <- data.frame(table(popData$class_level_name, popData$class_name))
  colnames(popCountData) <- c("class_level", "class_name", "freq")
  popCountData <- popCountData[popCountData$freq != 0,]

  print("...finished population list query...")
  print("...colating data...")

  # merge count data
  mergeCountData <- merge(popCountData, metsCountData, by=c("class_level", "class_name"), all=TRUE)
  mergeCountData[is.na(mergeCountData)] <- 0
  colnames(mergeCountData)[3:4] = c("pop_count", "mets_count")
  mergeCountData$fract_of_pop <- mergeCountData$mets_count/mergeCountData$pop_count
  mergeCountData <- mergeCountData[order(-mergeCountData$pop_count),]

  classes <- sort(unique(mergeCountData$class_level))
  resultSummary <- list()

  for (className in classes) {
    subTable <- mergeCountData[mergeCountData$class_level == className,]
    subTable$fract_within_pop <- subTable$pop_count / sum(subTable$pop_count)
    subTable$fract_within_mets <- subTable$mets_count / sum(subTable$mets_count)
    resultSummary[[className]] <- subTable
  }

  print("...creating query efficiency summary...")
  result <- list()
  result[["count_summary"]] <- resultSummary
  result[["met_classes"]] <- metsData
  result[["pop_classes"]] <- popData

  result[["query_report"]] = list()
  result[["query_report"]][["met_query_report"]] = metQueryReport
  result[["query_report"]][["pop_query_report"]] = popQueryReport

  return(result)
}


# runs and tallies chemical class information when there is NOT a user defined population, uses the DB population
chemicalClassSurveyRampIdsFullPopConn <- function(mets, conn) {

  mets <- unique(mets)

  checkIdPrefixes(mets)

  result <- list()

  # first handle metabolites of interest
  metStr <- paste(mets, collapse = "','")
  metStr <- paste("'" ,metStr, "'", sep = "")

  sql <- paste("select distinct a.ramp_id, b.sourceId, a.class_level_name, a.class_name, a.source from metabolite_class a, source b
          where b.rampId = a.ramp_id and b.sourceId in (",metStr,")")

  metsData <- RMariaDB::dbGetQuery(conn, sql)

  # need to filter for our specific source ids
  metsData <- subset(metsData, sourceId %in% mets)

  # get query summary
  metQueryReport <- queryReport(mets, metsData$sourceId)

  metsCountData <- data.frame(table(metsData$class_level_name,metsData$class_name))
  colnames(metsCountData) <- c("class_level", "class_name", "freq")
  metsCountData <- metsCountData[metsCountData$freq != 0,]
  metsCountData <- metsCountData[order(-metsCountData$freq),]

  print("...finished metabolite list query...")

  sql <- paste("select class_level_name, class_name, count(1) as pop_hits from metabolite_class
                 group by class_level_name, class_name")

  popCountData <- RMariaDB::dbGetQuery(conn, sql)

  colnames(popCountData) <- c("class_level", "class_name", "freq")
  popCountData <- popCountData[popCountData$freq != 0,]

  print("...finished DB population query...")
  print("...colating data...")

  # merge count data
  mergeCountData <- merge(popCountData, metsCountData, by=c("class_level", "class_name"), all=TRUE)
  mergeCountData[is.na(mergeCountData)] <- 0
  colnames(mergeCountData)[3:4] = c("pop_count", "mets_count")
  mergeCountData$fract_of_pop <- mergeCountData$mets_count/mergeCountData$pop_count
  mergeCountData <- mergeCountData[order(-mergeCountData$pop_count),]

  classes <- sort(unique(mergeCountData$class_level))
  resultSummary <- list()

  for (className in classes) {
    subTable <- mergeCountData[mergeCountData$class_level == className,]
    subTable$fract_within_pop <- subTable$pop_count / sum(subTable$pop_count)
    subTable$fract_within_mets <- subTable$mets_count / sum(subTable$mets_count)
    resultSummary[[className]] <- subTable
  }

  print("...creating query efficiency summary...")
  result <- list()
  result[["count_summary"]] <- resultSummary
  result[["met_classes"]] <- metsData

  result[["query_report"]] = list()
  result[["query_report"]][["met_query_report"]] = metQueryReport

  return(result)
}


# reports on how well the query performed
queryReport <- function(queryList, foundList) {
  querySummary = list()
  querySummary[["query_list_size"]] <- length(unique(queryList))
  querySummary[["found_list_size"]] <- length(unique(foundList))
  querySummary[["missed_query_elements"]] <- setdiff(queryList, foundList)
  querySummary
}


# reports on the total number of metabolites found within a collection of classes
# the tallies are for the population list and the met list
getTotalFoundInCategories <- function(classData) {
  counts <- list()

  # met list values for each class category
  counts[["mets"]] <- table(classData$met_classes$class_level_name)

  # pop list values for each class category
  # two routes depending on wether we have a user provide population or the all-DB population
  if(!is.null(classData$pop_classes)) {
    # if we have a population classes object use 'table' to grab the tally
    counts[["pop"]] <- table(classData$pop_classes$class_level_name)
  } else {

    # merge the results for different class categories
    cntSum <- classData$count_summary[[1]]
    if(length(classData$count_summary) >1) {
      for(m in 2:length(classData$count_summary)) {
        cntSum <- rbind(cntSum,classData$count_summary[[m]])
      }
    }

    # use aggregate to get the count tally for each class category (instead of 'table')
    cSum <- aggregate(cntSum$pop_count, by=list(cntSum$class_level), FUN=sum)
    cNames <- cSum$Group.1
    cSum <- data.frame(t(cSum$x))
    colnames(cSum) <- cNames
    counts[["pop"]] <- cSum
  }
  return(counts)
}


# returns a p-value ordered matrix with adjP_BH column added
bhCorrect <- function(resultMat) {
  resultMat <- resultMat[order(resultMat$`p-value`),]
  bhPvals <- p.adjust(resultMat$`p-value`, method = "BH")
  resultMat$adjP_BH <- bhPvals
  return(resultMat)
}



