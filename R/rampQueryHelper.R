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
#' @return a data frame that has all source Id in the column or the source table that has metabolites entry
rampFindSourceFromId <- function(rampId = "",full = TRUE){
    if(rampId == ""){
        stop("Data must be a list or dataframe")
    }
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


#'chemicalClassSurveyRampIdsConn is a helper function that takes a list of metabolite ids, a list of 'population' metabolite ids
#' and a MariaDB Connection object. The method returns metabolite class information for the metabolite list and a population of all ramp metabolites.
#' @param mets a list object of prefixed metabolite ids of interest
#' @param pop a list object of prefixed metabolite ids, representing a larger population of metabolites from which the mets were selected.
#' @param conn a MariaDB Connection object to support queries
#' @returns a list object containing three objects 'count_summary', 'met_classes' and 'met_query_report'.
#' The count_summary is a dataframe containing metabolite classes and number of metabolites in each class.
#' The met_classes is a detailed listing of compound classes associated with each input metabolite
#' The met_query_report indicates the number of input metabolites, how many were found in the DB and the list of metabolites not found in RaMP DB.
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


#'chemicalClassSurveyRampIdsFullPopConn is a helper function that takes a list of metabolite ids and a MariaDB Connection object
#'and returns metabolite class information for the metabolite list and a population of all ramp metabolites.
#' @param mets a list object of prefixed metabolite ids of interest
#' @param conn a MariaDB Connection object to support queries
#' @returns a list object containing three objects 'count_summary', 'met_classes' and 'met_query_report'.
#' The count_summary is a dataframe containing metabolite classes and number of metabolites in each class.
#' The met_classes is a detailed listing of compound classes associated with each input metabolite
#' The met_query_report indicates the number of input metabolites, how many were found in the DB and the list of metabolites not found in RaMP DB.
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


#' queryReport is a helper function to report on the number of query items that were found and missed, and the list of missed query values.
#' @param queryList is a list object that contains all user input query values
#' @param foundList is a list object of all user input query values that were retrieved during the query
#' @returns returns a list object with three return values, 'query_list_size', 'found_list_size', 'missed_query_elements'
#' The 'size' values are integers for the size of the input query and the number of input query values found.
#' The missed_query_elements is a list containing the subset of query values that are not found during the query.
queryReport <- function(queryList, foundList) {
  querySummary = list()
  querySummary[["query_list_size"]] <- length(unique(queryList))
  querySummary[["found_list_size"]] <- length(unique(foundList))
  querySummary[["missed_query_elements"]] <- setdiff(queryList, foundList)
  querySummary
}

#' Utility method to return metabolite counts found in compound class categories
#' based on an input data compound class data object from the chemicalClassSurvey function
#' The returned counts for each class category are for both the metabolite id query list
#' and for the larger full or user-defined population of metabolite ids.
#' This method is used in the exported chemicalClassEnrichment function
#' @param classData Data object returned from a call to chemicalClassSurvey
#' This input contains lists of chemical classes that pertain to a query list of metabolites and pertaining to
#' metabolites in a larger metabolite population.
#' @returns a list object with two keys, 'mets' and 'pop' that each has a table of metabolite or population
#' chemical classes and metabolite counts per class. This supports the chemicalClassEnrichment function.
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


#' Helper function to return a Benjamini-Hochberg (BH) corrected p-values. The method supports
#' chemical class enrichment p-value corrections.
#' @param resultMat is a dataframe containing p-values in a column named 'p-value'
#' @returns the input dataFrame with BH corrected p-values with column name 'adjP_BH'
bhCorrect <- function(resultMat) {
  resultMat <- resultMat[order(resultMat$`p-value`),]
  bhPvals <- stats::p.adjust(resultMat$`p-value`, method = "BH")
  resultMat$adjP_BH <- bhPvals
  return(resultMat)
}

#' Get class info for an input of metabolite source Ids
#' @param sourceIds a vector of analytes (genes or metabolites) that need to be searched
#' @return a dataframe of chemClass info
rampFindClassInfoFromSourceId<-function(sourceIds){
    sourceIds <- unique(sourceIds)
    checkIdPrefixes(sourceIds)
    idsToCheck <- sapply(sourceIds,function(x){
        if(!grepl("hmdb|chebi|LIPIDMAPS",x)){
            return(x)
        }
    })
    idsToCheck <- paste(idsToCheck, collapse = "','")
    idsToCheck <- paste("'" ,idsToCheck, "'", sep = "")
    conn <- connectToRaMP()
    sql <- paste("select * from source where sourceId in (",idsToCheck,")")

    potentialMultiMappings <- RMariaDB::dbGetQuery(conn, sql)
    potentialMultiMappings <- potentialMultiMappings %>%
        dplyr::select("sourceId","rampId") %>%
        dplyr::distinct()

    multimapped<-duplicated(potentialMultiMappings$sourceId)
    sourceIds<-sapply(sourceIds,function(x){
        ifelse(x %in% multimapped, return("Ambiguous"),return(x))
    })
    if("Ambiguous" %in% sourceIds){
        noAmbiguous = length(which(sourceIds=="Ambiguous"))
        print(paste0(noAmbiguous,
                     " metabolite(s) could not be unambiguously mapped to a chemical structure and have been discarded"))
    }

                                        # first handle metabolites of interest
    metStr <- paste(sourceIds, collapse = "','")
    metStr <- paste("'" ,metStr, "'", sep = "")

    conn <- connectToRaMP()

    sql <- paste("select distinct a.ramp_id, b.sourceId, a.class_level_name, a.class_name, a.source from metabolite_class a, source b
          where b.rampId = a.ramp_id and b.sourceId in (",metStr,")")

    metsData <- RMariaDB::dbGetQuery(conn, sql)

                                        # need to filter for our specific source ids
    metsData <- subset(metsData, "sourceId" %in% sourceIds)

    DBI::dbDisconnect(conn)
    return(metsData)
}

#' Internal function for extracting annotations, used by pathway and chemical enrichment test functions
#' @param analytes a vector of analytes (genes or metabolites) that need to be searched
#' @param PathOrChem return "path" information for pathways or "chem" for chemical class
#' @param find_synonym find all synonyms or just return same synonym (T/F)
#' @param NameOrIds whether input is "names" or "ids" (default is "ids")
#' @return a list of rampIds for "path" or a dataframe of chemClass info
getRaMPInfoFromAnalytes<-function(analytes,
                                  NameOrIds = "ids",
                                  PathOrChem = "path",
                                  find_synonym = FALSE){
    if(PathOrChem == "path"){
        if(NameOrIds == "names"){
            synonym <- rampFindSynonymFromSynonym(synonym=analytes,
                                                  find_synonym=find_synonym)

            synonym <- data.frame(synonym)
            colnames(synonym)[1]="commonName"
            synonym$commonName <- tolower(synonym$commonName)
            print(dim(synonym))
            if(nrow(synonym)==0) {
                stop("Could not find any matches to the analytes entered.  If pasting, please make sure the names are delimited by end of line (not analyte per line)\nand that you are selecting 'names', not 'ids'");
            }
            return(synonym)
        } else if (NameOrIds == "ids"){
            sourceramp <- rampFindSourceRampId(sourceId=analytes)
            if (nrow(sourceramp)==0) {
                warning("Make sure you are actually inputting ids and not names (you have NameOrIds set to 'ids'. If you are, then no ids were matched in the RaMP database.")
		return(NULL)
            } else {
	            return(sourceramp)
   	    }
        } else {
            stop("Make sure NameOrIds is set to 'names' or 'ids'")
        }
    }else if(PathOrChem == "chem"){
        if(NameOrIds == "names"){
            stop("Please do not use common names when searching for chemical structures, as a single name often refers to many structures")
        }else if (NameOrIds == "ids"){
            chem_info<-rampFindClassInfoFromSourceId(sourceIds=analytes)
            return(chem_info)
        }
    }else{
        stop("'PathOrChem' must be set to 'path' or 'chem'")
    }
}

##' Return all analytes that map to unique pathways in a pathwaydf
##' @param inputdf internal df with pathwayramp ids
##' @return dataframe of all analytes that map to the input pathways
##' @author Andrew Christopher Patt
buildFrequencyTables<-function(inputdf){
    ## Get pathway ids that contain the user analytes
    pid <- unique(inputdf$pathwayRampId);
    list_pid <- sapply(pid,shQuote)
    list_pid <- paste(list_pid,collapse = ",")

    ## Retrieve compound ids associated with background pathways and count
    query <- paste0("select * from analytehaspathway where pathwayRampId in (",
                     list_pid,")")

    con <- connectToRaMP()

    input_RampIds <- DBI::dbGetQuery(con,query)
    DBI::dbDisconnect(con)

    return(input_RampIds)
}

##' Separate input ids into lists based on database of origin
##' @param input_RampIds list of analyte ramp IDs
##' @return pathway lists separated by source db
##' @author Andrew Christopher Patt
segregateDataBySource<-function(input_RampIds){
    # data frames for metabolites with pathwayRampID, Freq based  on Source(kegg, reactome, wiki)

    input_RampId_C <- input_RampIds[grep("RAMP_C", input_RampIds$rampId), ]
    unique_input_RampId_C <- unique(input_RampId_C[,c("rampId", "pathwayRampId")])
    unique_pathwayRampId_source <- unique(input_RampId_C[,c("pathwayRampId", "pathwaySource")])

    freq_unique_input_RampId_C <- as.data.frame(table(unique_input_RampId_C[,"pathwayRampId"]))

    names(freq_unique_input_RampId_C)[1] = 'pathwayRampId'
    merge_Pathwayfreq_source <- merge(freq_unique_input_RampId_C, unique_pathwayRampId_source, by="pathwayRampId")

    # subset metabolite data based on source -  kegg, reactome, wiki

    input_kegg_metab <- subset(merge_Pathwayfreq_source, merge_Pathwayfreq_source$pathwaySource == "kegg")
    input_reactome_metab <- subset(merge_Pathwayfreq_source, merge_Pathwayfreq_source$pathwaySource == "reactome")
    input_wiki_metab <- subset(merge_Pathwayfreq_source, merge_Pathwayfreq_source$pathwaySource == "wiki")

# data frames for Genes with pathawayRampID, Freq based  on Source(kegg, reactome, wiki, hmdb)

    input_RampId_G <- input_RampIds[grep("RAMP_G", input_RampIds$rampId), ]
    unique_input_RampId_G <- unique(input_RampId_G[,c("rampId", "pathwayRampId")])
    unique_pathwayG_source <- unique(input_RampId_G[,c("pathwayRampId", "pathwaySource")])

    freq_unique_input_RampId_G <- as.data.frame(table(unique_input_RampId_G[,"pathwayRampId"]))

    names(freq_unique_input_RampId_G)[1] = 'pathwayRampId'
    merge_PathwayG_source <- merge(freq_unique_input_RampId_G, unique_pathwayG_source, by="pathwayRampId")

   # subset gene data based on source -  kegg, reactome, wiki

    input_kegg_gene <- subset(merge_PathwayG_source, merge_PathwayG_source$pathwaySource == "kegg")
    input_reactome_gene <- subset(merge_PathwayG_source, merge_PathwayG_source$pathwaySource == "reactome")
    input_wiki_gene <- subset(merge_PathwayG_source, merge_PathwayG_source$pathwaySource == "wiki")
    return(
        list(list(input_kegg_metab, input_reactome_metab, input_wiki_metab),
             list(input_kegg_gene, input_reactome_gene, input_wiki_gene)))
}
