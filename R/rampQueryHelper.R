#' Find all synonym from a given metabolite's name
#' This function is used to filter out some super common synonyms like glyceride
#' Now, this function only format the user input, so the user vector, dataframe,
#' and entire string separated by comma are working.
#' @param synonym name to search for
#' @param full bool if return whole data.frame
#' @param find_synonym bool if find all synonyms or just return same synonym
#' as input (there are some common synonyms that will mess up whole searching)
#' @return a data frame that contains synonym in the first column rampId in the second column
rampFindSynonymFromSynonym <- function(synonym,full = FALSE,
	find_synonym = FALSE){

  if(is.character(synonym)){
    if(grepl("\n",synonym)[1]){
      list_metabolite <- strsplit(synonym,"\n")
      list_metabolite <- unlist(list_metabolite)
    } else if(grepl(",",synonym)[1]){
      list_metabolite <- strsplit(synonym,",")
      list_metabolite <- unlist(list_metabolite)
    } else {
      list_metabolite <- synonym
    }
  } else if(is.data.frame(synonym)){
    list_metabolite <- unlist(synonym)
  } else{
    message("Wrong Format of argument")
    return(NULL)
  }
  if(!find_synonym){
    #message("Dont Find synonym due to common synonyms (Triglyceride?).")
    list_metabolite <- unique(list_metabolite)
    list_metabolite <- sapply(list_metabolite,shQuote)
    list_metabolite <- paste(list_metabolite,collapse = ",")
    query <- paste0("select Synonym as origins,rampId from analytesynonym where Synonym in(",
                    list_metabolite,
                    ");")
  con <- connectToRaMP()
    df1 <- DBI::dbGetQuery(con,query)
    DBI::dbDisconnect(con)
    return(df1)
  }
  list_metabolite <- unique(list_metabolite)
  list_metabolite <- sapply(list_metabolite,shQuote)
  list_metabolite <- paste(list_metabolite,collapse = ",")
  query <- paste0("select Synonym as origins,rampId from analytesynonym where Synonym in(",
                  list_metabolite,
                  ");")
  con <- connectToRaMP()

  df1 <- DBI::dbGetQuery(con,query)
  DBI::dbDisconnect(con)
  rampid <- df1$rampId
  rampid <- sapply(rampid,shQuote)
  rampid <- paste(rampid,collapse = ",")
  query <- paste0("select * from analytesynonym where rampId in(",rampid,");")
  con <- connectToRaMP()
  df2 <- DBI::dbGetQuery(con,query)
  DBI::dbDisconnect(con)
  df2 <- merge(df1,df2)
  if(full){
    return(df2)
  }
  synonym <- df2$Synonym
  synonym
}
#' Find all source from given list of RaMP Ids
#' @param rampId could be a data frame return by rampFindSynonymFromSynonym
#' containing all information related to synonym. Or can be a list of
#' rampId
#' @param full return whole searching result or not (TRUE/FALSE)
#' @return a data frame that has all source Id in the column or the source table that has metaoblites entry
rampFindSourceFromId <- function(rampId=NULL,full = TRUE){

  if(is.data.frame(rampId)){
    list_id <- rampId$rampId
  } else if(is.character(rampId)){
    if(grepl("\n",rampId)[1]){
      list_id <- strsplit(rampId,"\n")
      list_id <- unlist(list_id)
    } else if(grepl(",",rampId)[1]){
      list_id <- strsplit(rampId,",")
      list_id <- unlist(list_id)
    } else {
      list_id <- rampId
    }
  } else{
    message("Wrong format of input")
    return(NULL)
  }
  list_id <- unique(list_id)
  list_id <- sapply(list_id,shQuote)
  list_id <- paste(list_id,collapse = ",")
  query <- paste0("select * from source where rampId in (",list_id,");")
  con <- connectToRaMP()
  df <- DBI::dbGetQuery(con,query)
  DBI::dbDisconnect(con)
  if(full){
    return(df)
  } else{
    return(df[,1])
  }
}


#' Fast search given a list of metabolites source Id
#' @param sourceid a vector of synonym that need to be searched
#' @param find_synonym bool if find all synonyms or just return same synonym
#' @return a list contains all metabolits as name and pathway inside.
rampFastPathFromSource<- function(sourceid,find_synonym = FALSE){
  # progress<- shiny::Progress$new()
  # progress$set(message = "Querying databases ...",value = 0)
  now <- proc.time()
  # on.exit(dbDisconnect(con))
  # find synonym

  #synonym <- rampFindSynonymFromSynonym(synonym,find_synonym=find_synonym)

  list_metabolite <- unique(sourceid)
  list_metabolite <- sapply(list_metabolite,shQuote)
  list_metabolite <- paste(list_metabolite,collapse = ",")
  query1 <- paste0("select * from source where sourceid in (",
                   list_metabolite,");")
  con <- connectToRaMP()
  df1<- DBI::dbGetQuery(con,query1)
  DBI::dbDisconnect(con)
  colnames(df1)[1] <-"sourceId2"
  #return(df1)
  rampid <- df1$rampId
  rampid <- sapply(rampid,shQuote)
  rampid <- paste(rampid,collapse = ",")
  query2 <- paste0("select * from analytehaspathway where
                   rampId in (",rampid,");")
  con <- connectToRaMP()
  df2 <- DBI::dbGetQuery(con,query2)
  DBI::dbDisconnect(con)
  #return(df2)
  id_list <- unique(df2$pathwayRampId)
  id_list <- sapply(id_list,shQuote)
  id_list <- paste(id_list,collapse = ",")
  print(id_list)
  query3 <- paste0("select * from pathway where pathwayRampId in (",
                   id_list,");")
  con <- connectToRaMP()
  df3 <- DBI::dbGetQuery(con,query3)
  DBI::dbDisconnect(con)
  #return(df3)
  mdf <- merge(df3,df2,all.x=T)
  mdf <- merge(mdf,df1,all.x = T)
  mdf <- unique(mdf)
  print("timing ...")
  print(proc.time()- now)
  return(unique(mdf[,c(6,4,3,7)]))
}
#' Find rampId from given source ID
#' The rampId can be plugged in other functions to continue query
#' @param sourceId a data frame or string separated by comma or string
#' separated by new line
#' @return data.frame that has sourceId and rampId and source as columns
rampFindSourceRampId <- function(sourceId){

  if(is.character(sourceId)){
    if(grepl("\n",sourceId)[1]){
      list_metabolite <- strsplit(sourceId,"\n")
      list_metabolite <- unlist(list_metabolite)
    } else if(grepl(",",sourceId)[1]){
      list_metabolite <- strsplit(sourceId,",")
      list_metabolite <- unlist(list_metabolite)
    } else {
      list_metabolite <- sourceId
    }
  } else if(is.data.frame(sourceId)){
    list_metabolite <-as.character(unlist(sourceId$sourceId))
  } else{
    message("Wrong Format of argument")
    return(NULL)
  }
  list_metabolite <- sapply(list_metabolite,shQuote)
  list_metabolite <- paste(list_metabolite,collapse = ",")
  con <- connectToRaMP()
  query <- paste0("select sourceId,IDtype as analytesource, rampId from source where sourceId in (",list_metabolite,");")
  df <- DBI::dbGetQuery(con,query)
  DBI::dbDisconnect(con)
  return(df)
}

#' Check for id prefixes
#' By default, the function checks on whether there are at least 90% of input user analytes with
#' appropriate RaMP-supported prefixes.
#' @param idList list of ids (typically input by user)
#' @param perc_cutoff percent cutof used to throw a warning that many ids are missing
checkIdPrefixes <- function(idList, perc_cutoff=0.9) {
  idCount <- length(idList)
  prefixCount <- 0
  for(id in idList) {
    if(grepl(":",id, fixed = TRUE)) {
      prefixCount <- prefixCount + 1
    }
  }
  if(prefixCount/idCount < perc_cutoff) {
    print(paste("RaMP expects ids to be prefixed with the source database." + (idCount-prefixCount) + " of " + idCount + " ids lack prefixes.\n", sep=""))
    #warnObj <- Warning(warn, call=TRUE, immediate=TRUE)
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
  sourceId <- NULL
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
  sourceId <- NULL
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
    cSum <- stats::aggregate(cntSum$pop_count, by=list(cntSum$class_level), FUN=sum)
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
  bhPvals <- stats::p.adjust(resultMat$`p-value`, method = "BH")
  resultMat$adjP_BH <- bhPvals
  return(resultMat)
}


