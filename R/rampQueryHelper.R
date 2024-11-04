#' Find all synonym from a given metabolite's name
#' This function is used to filter out some super common synonyms like glyceride
#' Now, this function only format the user input, so the user vector, dataframe,
#' and entire string separated by comma are working.
#' @param synonym name to search for
#' @param full bool if return whole data.frame
#' @param returnRampIds bool to return ramp Ids with output
#' (there are some common synonyms that will mess up whole searching)
#' @param db a RaMP database object
#' @return a data frame that contains synonym in the first column rampId in the second column
#' @noRd
rampFindSynonymFromSynonym <- function( synonym, full = FALSE, returnRampIds = FALSE, db = RaMP()){
  if(is.character(synonym)){
    if(grepl("\n",synonym)[1]){
      list_metabolite <- strsplit(synonym,"\n")
      list_metabolite <- unlist(list_metabolite)
    ## } else if(grepl(",",synonym)[1]){
    ##   list_metabolite <- strsplit(synonym,",")
    ##   list_metabolite <- unlist(list_metabolite)
    } else {
      list_metabolite <- synonym
    }
  } else if(is.data.frame(synonym)){
    list_metabolite <- unlist(synonym)
  } else{
    message("Wrong Format of argument")
    return(NULL)
  }

  list_metabolite <- unique(list_metabolite)

  df1 <- db@api$getSynonymsForSynonym(list_metabolite)

  if(returnRampIds || nrow(df1) < 1) {
      return(df1)
  } else {
      rampid <- df1$rampId
      df2 <- db@api$getSynonymInfoForRampIDs(rampIds=rampid)

      df2 <- merge(df1,df2)
      if(full){
          return(df2)
      }
      synonym <- df2$Synonym
      return(synonym)
  }
}
#' Find all source from given list of RaMP Ids
#' @param rampId could be a data frame return by rampFindSynonymFromSynonym
#' containing all information related to synonym. Or can be a list of
#' rampId
#' @param full return whole searching result or not (TRUE/FALSE)
#' @param db a RaMP database object
#' @return a data frame that has all source Id in the column or the source table that has metabolites entry
#' @noRd
rampFindSourceFromId <- function(rampId = "",full = TRUE, db = RaMP()){
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
  df <- db@api$getSourceInfoForRampIDs(rampIds = list_id)

  if(full){
    return(df)
  } else{
    return(df[,1])
  }
}


#' Fast search given a list of metabolites source Id
#' @param sourceId a vector of synonym that need to be searched
#' @param db a RaMP database object
#' @return a list contains all metabolits as name and pathway inside.
#' @noRd
rampFastPathFromSource<- function( sourceId, db = RaMP()){
  # progress<- shiny::Progress$new()
  # progress$set(message = "Querying databases ...",value = 0)
  now <- proc.time()
  # on.exit(dbDisconnect(con))

  list_metabolite <- unique(sourceId)

  df1 <- db@api$getAllSourceInfoForSourceIDs(sourceIds = list_metabolite)

  colnames(df1)[1] <-"sourceId2"
  #return(df1)
  rampid <- df1$rampId
  df2 <- db@api$getAllPathwaysForRampIDs(rampIds = rampid)

  #return(df2)
  id_list <- unique(df2$pathwayRampId)
  df3 <- db@api$getPathwayInfoForRampIDs(pathwayRampIds = id_list)

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
#' @param db a RaMP database object
#' @return data.frame that has sourceId and rampId and source as columns
#' @noRd
rampFindSourceRampId <- function( sourceId, db = RaMP()){

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
  df <- db@api$getSourceInfoFromSourceIDs(sourceIds = list_metabolite)

  return(df)
}

#' Check for id prefixes
#' By default, the function checks on whether there are at least 90% of input user analytes with
#' appropriate RaMP-supported prefixes.
#' @param idList list of ids (typically input by user)
#' @param percCutoff percent cutof used to throw a warning that many ids are missing
#' @noRd
checkIdPrefixes <- function(idList, percCutoff=0.9) {
  idCount <- length(idList)
  prefixCount <- 0
  for(id in idList) {
    if(grepl(":",id, fixed = TRUE)) {
      prefixCount <- prefixCount + 1
    }
  }
  if(prefixCount/idCount < percCutoff) {
    print(paste("RaMP expects ids to be prefixed with the source database." + (idCount-prefixCount) + " of " + idCount + " ids lack prefixes.\n", sep=""))
    #warnObj <- Warning(warn, call=TRUE, immediate=TRUE)
    print("Common metabolite prefixes: CAS:, chebi:, chemspider:, hmdb:, kegg:, LIPIDMAPS:, pubchem:")
    print("Examples: kegg:C02712, hmdb:HMDB04824, CAS:2566-39-4. The input list may contain a variety of id types.")
  }
}


#' queryReport is a helper function to report on the number of query items that were found and missed, and the list of missed query values.
#' @param queryList is a list object that contains all user input query values
#' @param foundList is a list object of all user input query values that were retrieved during the query
#' @returns returns a list object with three return values, 'query_list_size', 'found_list_size', 'missed_query_elements'
#' The 'size' values are integers for the size of the input query and the number of input query values found.
#' The missed_query_elements is a list containing the subset of query values that are not found during the query.
#' @noRd
queryReport <- function(queryList, foundList) {
  querySummary = list()
  querySummary[["query_list_size"]] <- length(unique(queryList))
  querySummary[["found_list_size"]] <- length(unique(foundList))
  querySummary[["missed_query_elements"]] <- setdiff(queryList, foundList)
  querySummary
}

#' Utility method to return metabolite counts found in compound class categories
#' based on an input data compound class data object from the getChemClass function
#' The returned counts for each class category are for both the metabolite id query list
#' and for the larger full or user-defined population of metabolite ids.
#' This method is used in the exported chemicalClassEnrichment function
#' @param classData Data object returned from a call to getChemClass
#' This input contains lists of chemical classes that pertain to a query list of metabolites and pertaining to
#' metabolites in a larger metabolite population.
#' @param inferIdMapping if FALSE, the survey only reports on class annotations made directly on the input ids.
#' If inferIdMapping is set to TRUE, the ids are cross-referenced or mapped to related ids that contain metabolite class annotations.
#' The default is FALSE. Note that this utility method is typcally used within chemical class enrichment and is passed the value of this parameter.
#' @returns a list object with two keys, 'mets' and 'pop' that each has a table of metabolite or population
#' chemical classes and metabolite counts per class. This supports the chemicalClassEnrichment function.
#' @noRd
getTotalFoundInCategories <- function(classData, inferIdMapping=FALSE) {
  counts <- list()

  print("check total summary")

  # met list values for each class category
  if(inferIdMapping) {
    metsData2 <- classData$met_classes
    metsData2 <- metsData2[,c('ramp_id','class_level_name', 'class_name', 'directIdClassHits')]
    metsData2 <- unique(metsData2)
    # need to sum direct hits on class levels by directIdClassHits
    metsClassLevelInfo = stats::aggregate(metsData2$directIdClassHits, by=list(metsData2$class_level_name), sum)
    colnames(metsClassLevelInfo) <- c("Var1","Freq")
    counts[["mets"]] <- metsClassLevelInfo

  } else {
    counts[["mets"]] <- data.frame(table(classData$met_classes$class_level_name))
  }
  # pop list values for each class category
  # two routes depending on whether we have a user provide population or the all-DB population
  if(!is.null(classData$pop_classes)) {
    # if we have a population classes object use 'table' to grab the tally

    if(inferIdMapping) {
      popData2 <- classData$pop_classes
      popData2 <- popData2[,c('ramp_id','class_level_name', 'class_name', 'directIdClassHits')]
      popData2 <- unique(popData2)

      popStats <- stats::aggregate(popData2$directIdClassHits, by=list(popData2$class_level_name),sum)
    } else {
      popStats <- data.frame(table(classData$pop_classes$class_level_name))
    }
    colnames(popStats) <- c("Var1", "Freq")
    counts[["pop"]] <- popStats
  } else {

    # merge the results for different class categories
    cntSum <- classData$count_summary[[1]]
    if(length(classData$count_summary) >1) {
      for(m in 2:length(classData$count_summary)) {
        cntSum <- rbind(cntSum,classData$count_summary[[m]])
      }
    }

    print("getting population totals")
    # use aggregate to get the count tally for each class category (instead of 'table')
    if(inferIdMapping) {
      cSum <- stats::aggregate(cntSum$pop_count, by=list(cntSum$class_level), FUN=sum)
    } else {
      cSum <- stats::aggregate(cntSum$pop_count, by=list(cntSum$class_level), FUN=sum)
    }
    colnames(cSum) <- c("Var1", "Freq")
    #cNames <- cSum$Group.1
    #cSum <- data.frame(t(cSum$x))
    #colnames(cSum) <- cNames
    counts[["pop"]] <- cSum
  }
  return(counts)
}


#' Helper function to return a Benjamini-Hochberg (BH) corrected p-values. The method supports
#' chemical class enrichment p-value corrections.
#' @param resultMat is a dataframe containing p-values in a column named 'p-value'
#' @returns the input dataFrame with BH corrected p-values with column name 'adjP_BH'
#' @noRd
bhCorrect <- function(resultMat) {
  resultMat <- resultMat[order(resultMat$`p-value`),]
  bhPvals <- stats::p.adjust(resultMat$`p-value`, method = "BH")
  resultMat$adjP_BH <- bhPvals
  return(resultMat)
}

#' Get class info for an input of metabolite source Ids
#' @importFrom rlang .data
#' @param sourceIds a vector of analytes (genes or metabolites) that need to be searched
#' @param db a RaMP database object
#' @return a dataframe of chemClass info
#' @noRd
rampFindClassInfoFromSourceId<-function(sourceIds, db = RaMP()){
    sourceIds <- unique(sourceIds)
    checkIdPrefixes(idList = sourceIds)
    idsToCheck <- Filter(function(x) !grepl("hmdb|chebi|LIPIDMAPS", x), sourceIds)

    potentialMultiMappings <- db@api$getAllSourceInfoForSourceIDs(sourceIds = idsToCheck)
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

    metsData <- db@api$getChemicalClassFromSourceIDs(sourceIds = sourceIds)

    metsData <- dplyr::filter(metsData, .data$sourceId %in% sourceIds)

    return(metsData)
}

#' Internal function for extracting annotations, used by pathway and chemical enrichment test functions
#' @param analytes a vector of analytes (genes or metabolites) that need to be searched
#' @param PathOrChem return "path" information for pathways or "chem" for chemical class
#' @param namesOrIds whether input is "names" or "ids" (default is "ids")
#' @param db a RaMP database object
#' @return a list of rampIds for "path" or a dataframe of chemClass info
#' @noRd
getRaMPInfoFromAnalytes<-function( analytes,
                                  namesOrIds = "ids",
                                  PathOrChem = "path", db = RaMP()){
    if(PathOrChem == "path"){
        if(namesOrIds == "names"){
            synonym <- rampFindSynonymFromSynonym(synonym=analytes,
                                                  returnRampIds=FALSE)

            colnames(synonym)[1]="commonName"
            synonym$commonName <- tolower(synonym$commonName)
            print(dim(synonym))
            if(nrow(synonym)==0) {
                stop("Could not find any matches to the analytes entered.  If pasting, please make sure the names are delimited by end of line (not analyte per line)\nand that you are selecting 'names', not 'ids'");
            }
            return(synonym)
        } else if (namesOrIds == "ids"){
            sourceramp <- rampFindSourceRampId(db = db, sourceId=analytes)
            if (nrow(sourceramp)==0) {
                warning("Make sure you are actually inputting ids and not names (you have namesOrIds set to 'ids'. If you are, then no ids were matched in the RaMP database.")
		return(NULL)
            } else {
	            return(sourceramp)
   	    }
        } else {
            stop("Make sure namesOrIds is set to 'names' or 'ids'")
        }
    }else if(PathOrChem == "chem"){
        if(namesOrIds == "names"){
            stop("Please do not use common names when searching for chemical structures, as a single name often refers to many structures")
        }else if (namesOrIds == "ids"){
            chem_info<-rampFindClassInfoFromSourceId(sourceIds=analytes)
            return(chem_info)
        }
    }else{
        stop("'PathOrChem' must be set to 'path' or 'chem'")
    }
}

##' Return all analytes that map to unique pathways in a pathwaydf
##' @param inputdf internal df with pathwayramp ids
##' @param pathwayDefinitions If "RaMP" (default), use pathway definitions within RaMP-DB. Else, supply path to gmx file containing custom pathway definitions. GMX files are a tab-separated format that contain one analyte set per column, with the name of the set in the first row, and constituent analytes in subsequent rows
##' @param analyteType "genes" or "metabolites"
##' @param db a RaMP database object
##' @return dataframe of all analytes that map to the input pathways
##' @author Andrew Christopher Patt
#' @importFrom utils head
#' @noRd
buildFrequencyTables<-function( inputdf, pathwayDefinitions="RaMP", analyteType, db = RaMP()) {

  if(pathwayDefinitions == "RaMP") {
	print("Now in buildFrequencyTables function")
	print(mode(inputdf))
	print(class(inputdf))
	head(inputdf$pathwayRampId)

    ## Get pathway ids that contain the user analytes
    pid <- unique(inputdf$pathwayRampId);
    inputRampIds <- db@api$getAllRampIDsForAllPathwayRampIDs(pathwayRampIds = pid)

    return(inputRampIds)
  } else {
    tryCatch(
      {
        if (analyteType == "metabolites") {
          pathwayDefinitions <- readxl::read_excel(pathwayDefinitions, sheet = 1)
        } else if (analyteType == "genes") {
          pathwayDefinitions <- readxl::read_excel(pathwayDefinitions, sheet = 2)
        }
      },
      error = function(e) {
        print("Pathway file could not be found or is improperly formatted. Please supply path to GMX file for custom pathway definitions")
      }
    )
    inputRampIds <- data.frame(rampId=character(),
                                pathwayRampId=character())
    pid <- unique(inputdf$pathwayRampId)

    for(i in pid){
      temp <- data.frame(rampId = pathwayDefinitions[,i][which(!is.na(pathwayDefinitions[,i])),],
                         pathwayRampId = i)
      colnames(temp) <- c("rampId","pathwayRampId")
      inputRampIds <- rbind(inputRampIds,temp)
    }
    inputRampIds$pathwaySource = "custom"
    if(analyteType == "metabolites"){
      inputRampIds$rampId = paste0("RAMP_C",inputRampIds$rampId)
    }else if(analyteType == "genes"){
      inputRampIds$rampId = paste0("RAMP_G",inputRampIds$rampId)
    }
    return(inputRampIds)
  }
}

##' Separate input ids into lists based on database of origin
##' @param inputRampIds list of analyte ramp IDs
##' @return pathway lists separated by source db
##' @author Andrew Christopher Patt
#' @noRd
segregateDataBySource<-function(inputRampIds){
  # data frames for metabolites with pathwayRampID, Freq based  on Source(kegg, reactome, wiki)
  input_RampId_C <- inputRampIds[grep("RAMP_C", inputRampIds$rampId), ]
  unique_input_RampId_C <- unique(input_RampId_C[,c("rampId", "pathwayRampId")])
  unique_pathwayRampId_source <- unique(input_RampId_C[,c("pathwayRampId", "pathwaySource")])

  freq_unique_input_RampId_C <- as.data.frame(table(unique_input_RampId_C[,"pathwayRampId"]))

  names(freq_unique_input_RampId_C)[1] = 'pathwayRampId'
  merge_Pathwayfreq_source_C <- merge(freq_unique_input_RampId_C,
                                    unique_pathwayRampId_source, by="pathwayRampId")

  # subset metabolite data based on source

  ## input_kegg_metab <- subset(merge_Pathwayfreq_source,
  ##                            merge_Pathwayfreq_source$pathwaySource == "kegg")
  ## input_reactome_metab <- subset(merge_Pathwayfreq_source,
  ##                                merge_Pathwayfreq_source$pathwaySource == "reactome")
  ## input_wiki_metab <- subset(merge_Pathwayfreq_source,
  ##                            merge_Pathwayfreq_source$pathwaySource == "wiki")
  ## input_custom_metab <- subset(merge_Pathwayfreq_source,
  ##                              merge_Pathwayfreq_source$pathwaySource == "custom")

  input_metab <- lapply(unique(merge_Pathwayfreq_source_C$pathwaySource),
                        function(x){
    return(subset(merge_Pathwayfreq_source_C,merge_Pathwayfreq_source_C$pathwaySource == x))
  })
  names(input_metab) <- unique(merge_Pathwayfreq_source_C$pathwaySource)

  # data frames for Genes with pathawayRampID, Freq based  on Source(kegg, reactome, wiki, hmdb)

  input_RampId_G <- inputRampIds[grep("RAMP_G", inputRampIds$rampId), ]
  unique_input_RampId_G <- unique(input_RampId_G[,c("rampId", "pathwayRampId")])
  unique_pathwayG_source <- unique(input_RampId_G[,c("pathwayRampId", "pathwaySource")])

  freq_unique_input_RampId_G <- as.data.frame(table(unique_input_RampId_G[,"pathwayRampId"]))

  names(freq_unique_input_RampId_G)[1] = 'pathwayRampId'
  merge_PathwayG_source <- merge(freq_unique_input_RampId_G, unique_pathwayG_source, by="pathwayRampId")

  # subset gene data based on source -  kegg, reactome, wiki

  ## input_kegg_gene <- subset(merge_PathwayG_source, merge_PathwayG_source$pathwaySource == "kegg")
  ## input_reactome_gene <- subset(merge_PathwayG_source, merge_PathwayG_source$pathwaySource == "reactome")
  ## input_wiki_gene <- subset(merge_PathwayG_source, merge_PathwayG_source$pathwaySource == "wiki")
  ## input_custom_gene <- subset(merge_PathwayG_source, merge_PathwayG_source$pathwaySource == "custom")
  input_gene <- lapply(unique(merge_PathwayG_source$pathwaySource), function(x){
    return(subset(merge_PathwayG_source,merge_PathwayG_source$pathwaySource == x))
  })
  names(input_gene) <- unique(merge_PathwayG_source$pathwaySource)
  out = list(input_metab,input_gene)
  names(out) <- c("metab","gene")
  return(out)
}

##' Return list of duplicate Wikipathway IDs from Reactome. This may be unnecessary in the future
##' @return List of duplicate Wikipathway IDs from Reactome.
##' @param db a RaMP database object
##' @author Andrew Patt
#' @noRd
find_duplicate_pathways <- function(db = RaMP()){
  .Deprecated("findDuplicatePathways")
  return(findDuplicatePathways(db = db))
}

##' Return list of duplicate Wikipathway IDs from Reactome. This may be unnecessary in the future
##' @return List of duplicate Wikipathway IDs from Reactome.
##' @param db a RaMP database object
##' @author John Braisted
#' @noRd
findDuplicatePathways <- function(db = RaMP()) {
  reactomePIDs <- db@api$getRampIdsForPathways(pathwayType = 'reactome')
  df2 <- db@api$getExactMatchingPathways()

  dupReturnList <- list(nrow(df2))
  # preference for reactome over wiki or kegg
  for(r in 1:nrow(df2)) {
    if(df2[r,1] %in% reactomePIDs[,1]) {
      dupReturnList[[r]] <- df2[r,2]
    } else if(df2[r,2] %in% reactomePIDs[,1]) {
      dupReturnList[[r]] <- df2[r,1]
    } else {
      dupReturnList[[r]] <- df2[r,2]
    }
  }

  return(unlist(dupReturnList))
}



#' Filter pathways by p-value cutoff for display and clustering
#' @param enrichResults The data frame generated by runFisherTest
#' @param pValType Specifies which p-value to use as the filter threshold.
#' Permitted values are 'pval' and 'fdr' for chemical class and pathway enrichment.
#' Pathway enrichment also includes an optional 'holm' value for holm p-value corrections. Default is 'fdr'.
#' @param pValCutoff return pathways where pValType p-values are < pValCutoff
#' @return list:[[1]]Dataframe with pathway enrichment results, only significant pathways
#' [[2]]analyte type
#' @examples
#' \dontrun{
#' analyteList <- c("MDM2", "TP53", "glutamate", "creatinine")
#'
#' fisher.results <- runCombinedFisherTest(analytes = analyteList, namesOrIds = 'names')
#' filtered.fisher.results <- FilterFishersResults(fisher.results, pValType='fdr', pValCutoff = 0.10)
#' }
#' @export
filterEnrichResults <- function(enrichResults, pValType = 'fdr', pValCutoff = 0.1) {

  print("Filtering Fisher Results...")

  result_type <- enrichResults$result_type

  if (pValType == "fdr")
  {
    pvalToFilter <- "Pval_FDR"
  }

  if (pValType == "holm")
  {
    pvalToFilter <- "Pval_Holm"
  }

  if (pValType == "pval")
  {
    if (result_type == 'pathway_enrichment' | result_type == 'reactionClass_enrichment')
    {
      if (enrichResults$analyteType == 'both')
      {
        pvalToFilter <- "Pval_combined"
      } else if (enrichResults$analyteType == 'genes') {
        pvalToFilter <- "Pval_Gene"
      } else if (enrichResults$analyteType == 'metabolites' | enrichResults$analyteType == 'chebi') {
        pvalToFilter <- "Pval_Metab"
      } else if (enrichResults$analyteType == 'uniprot') {
        pvalToFilter <- "Pval_Prot"
      }
    } else if (result_type == 'ontology_enrichment' | result_type == 'chemical_class_enrichment')
    {
      pvalToFilter <- "Pval"
    }
  }

  for (i in 1:length(enrichResults))
  {
    if (is(enrichResults[[i]], 'data.frame'))
    {
      resultDf <- enrichResults[[i]]
      resultDf <- subset(resultDf, resultDf[[pvalToFilter]] <= pValCutoff)
      enrichResults[[i]] <- resultDf
    }
  }

  return(enrichResults)
}


#'getChemicalClassRampIdsConn is a helper function that takes a list of metabolite ids, a list of 'population' metabolite ids
#' and a MariaDB Connection object. The method returns metabolite class information for the metabolite list and a population of all ramp metabolites.
#' @importFrom rlang .data
#' @param mets a list object of prefixed metabolite ids of interest
#' @param pop a list object of prefixed metabolite ids, representing a larger population of metabolites from which the mets were selected.
#' @param inferIdMapping if FALSE, the survey only reports on class annotations made directly on the input ids.
#' If inferIdMapping is set to TRUE, the ids are cross-referenced or mapped to related ids that contain metabolite class annotations.
#' The default is TRUE.
#' @returns a list object containing three objects 'count_summary', 'met_classes' and 'met_query_report'.
#' The count_summary is a dataframe containing metabolite classes and number of metabolites in each class.
#' The met_classes is a detailed listing of compound classes associated with each input metabolite
#' The met_query_report indicates the number of input metabolites, how many were found in the DB and the list of metabolites not found in RaMP DB.
#' @param db a RaMP database object
#' @noRd
getChemicalClassRampIdsConn <- function( mets, pop, inferIdMapping=TRUE, db = RaMP()) {

  mets <- unique(mets)

  checkIdPrefixes(idList = mets)

  pop <- unique(pop)

  checkIdPrefixes(idList = pop)

  result <- list()

  metsData <- db@api$getClassesForAnalytes(analytes = mets, inferIdMapping = inferIdMapping, includeAnalyteName = TRUE)

  # need to filter for our specific source ids
  metsData <- dplyr::filter(metsData, .data$sourceId %in% mets)

  # get query summary
  metQueryReport <- queryReport(queryList = mets, foundList = metsData$sourceId)

  if(inferIdMapping) {
    # if inferring mapping through ramp ids, the count has to be reduced to only counting source ids from the metabolite_class table
    metsData2 <- unique(metsData[,c("ramp_id","class_level_name","class_name","directIdClassHits")])
    metsCountData <- stats::aggregate(metsData2$directIdClassHits, by=list(metsData2$class_level_name, metsData2$class_name), sum)
    colnames(metsCountData) <- c("class_level", "class_name", "freq")
  } else {
    metsCountData <- data.frame(table(metsData$class_level_name,metsData$class_name))
    colnames(metsCountData) <- c("class_level", "class_name", "freq")
    metsCountData <- metsCountData[metsCountData$freq != 0,]
    metsCountData <- metsCountData[order(-metsCountData$freq),]
  }

  print("...finished metabolite list query...")

  popData <- db@api$getClassesForAnalytes(analytes = pop, inferIdMapping = inferIdMapping, includeAnalyteName = FALSE)

  popData <- dplyr::filter(popData, .data$sourceId %in% pop)

  # get query summary
  popQueryReport <- queryReport(queryList = pop, foundList = popData$sourceId)


  if(inferIdMapping) {
    # if inferring mapping through ramp ids, the count has to be reduced to only counting source ids from the metabolite_class table
    popData2 <- unique(popData[,c("ramp_id","class_level_name","class_name","directIdClassHits")])
    popCountData <- stats::aggregate(popData2$directIdClassHits, by=list(popData2$class_level_name, popData2$class_name), sum)
  } else {
    popCountData <- data.frame(table(popData$class_level_name,popData$class_name))
  }

  # popCountData <- data.frame(table(popData$class_level_name, popData$class_name))
  colnames(popCountData) <- c("class_level", "class_name", "freq")
  popCountData <- popCountData[popCountData$freq != 0,]

  print("...finished population list query...")
  print("...collating data...")

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


#'getChemicalClassRampIdsFullPopConn2 is a helper function that takes a list of metabolite ids and a MariaDB Connection object
#'and returns metabolite class information for the metabolite list and a population of all ramp metabolites.
#' @importFrom rlang .data
#' @param mets a list object of prefixed metabolite ids of interest
#' @param inferIdMapping if FALSE, the survey only reports on class annotations made directly on the input ids.
#' If inferIdMapping is set to TRUE, the ids are cross-referenced or mapped to related ids that contain metabolite class annotations.
#' The default is TRUE.
#' @param db a RaMP database object
#' @returns a list object containing three objects 'count_summary', 'met_classes' and 'met_query_report'.
#' The count_summary is a dataframe containing metabolite classes and number of metabolites in each class.
#' The met_classes is a detailed listing of compound classes associated with each input metabolite
#' The met_query_report indicates the number of input metabolites, how many were found in the DB and the list of metabolites not found in RaMP DB.
#' @noRd
getChemicalClassRampIdsFullPopConn <- function( mets, inferIdMapping=TRUE, db = RaMP()) {

  mets <- unique(mets)

  checkIdPrefixes(idList = mets)

  result <- list()

  metsData <- db@api$getClassesForAnalytes(analytes = mets, inferIdMapping = inferIdMapping, includeAnalyteName = TRUE)

  # need to filter for our specific source ids
  metsData <- dplyr::filter(metsData, .data$sourceId %in% mets)

  # get query summary
  metQueryReport <- queryReport(queryList = mets, foundList = metsData$sourceId)

  emptyMetsResult = FALSE

  if(nrow(metsData) > 0) {

    if(inferIdMapping) {
      # if inferring mapping through ramp ids, the count has to be reduced to only counting source ids from the metabolite_class table
      metsData2 <- unique(metsData[,c("ramp_id","class_level_name","class_name","directIdClassHits")])
      metsCountData <- stats::aggregate(metsData2$directIdClassHits, by=list(metsData2$class_level_name, metsData2$class_name), sum)
    } else {
      metsCountData <- data.frame(table(metsData$class_level_name,metsData$class_name))
    }

    colnames(metsCountData) <- c("class_level", "class_name", "freq")
    metsCountData <- metsCountData[metsCountData$freq != 0,]
    metsCountData <- metsCountData[order(-metsCountData$freq),]

    print("...finished metabolite list query...")

  } else {

    emptyMetsResult = TRUE

    print("...finished metabolite list query, Warning: NO query term matches in RaMP DB...")
    # build and empty result for the mets data
    metsCountData <- data.frame(matrix(ncol=3, nrow=0))
    colnames(metsCountData) <- c("class_level", "class_name", "freq")
  }

  # get full population counts for all classes
  popCountData <- db@api$getMetaboliteCountsForClasses()

  colnames(popCountData) <- c("class_level", "class_name", "freq")
  popCountData <- popCountData[popCountData$freq != 0,]

  print("...finished DB population query...")
  print("...collating data...")

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

    # handle NAs more gracefully
    subTable[is.na(subTable)] <- 0

    # append result for the class category
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


#' listToQueryString utility method to convert an id list to a comma separate string, with single quoted values.
#'
#' @param ids list of ids (can be names or ids)
#'
#' @return comma separated list of single quoted ids or names
#' @noRd
listToQueryString <- function(ids) {
  isStr <- paste0("'", paste0(ids, collapse = "','"), "'", sep="")
  return (isStr)
}


#' filterPathwaysByAnalyteCount utility method filtered a dataframe based on the number of analytes associated with rampPathwayIds contained in the dataframe.
#' Like fisher exact code, this one retains pathways with analyte count >= minPathwaySize, and having analyte_count < max_path_size
#'
#' @param pathwayDataframe a dataframe containing at least one column that contains rampPathwayIds
#' @param pathwayRampIdColName the column name containing the rampPathwayIds
#' @param minPathwaySize the minimum number of pathway members (genes and metabolites) to include the pathway in the output (default = 5)
#' @param maxPathwaySize the maximum number of pathway members (genes and metaboltes) to include the pathway in the output (default = 150)
#' @param db a RaMP database object
#' @noRd
filterPathwaysByAnalyteCount <- function( pathwayDataframe, pathwayRampIdColName = 'pathwayRampId', minPathwaySize = 5, maxPathwaySize = 150, db = RaMP()) {
  pwIds <- unlist(pathwayDataframe[[pathwayRampIdColName]])

  res <- db@api$getAnalyteCountsForPathways(pathwayRampIds = pwIds)
  res <- res[res$analyte_count >= minPathwaySize & res$analyte_count < maxPathwaySize,]
  keeperPW <- unlist(res$pathwayRampId)
  pathwayDataframe <- pathwayDataframe[pathwayDataframe[[pathwayRampIdColName]] %in% keeperPW, ]
  return(pathwayDataframe)
}


#' Creates the input dataframe for the sunburst plot created in 'plotReactionClasses'
#'
#' @param reactionClassesResults output of getReactionClassesForAnalytes()
#' @importFrom grDevices adjustcolor
#' @noRd
buildReactionClassesSunburstDataframe <- function(reactionClassesResults) {

  #create empty table for sunburst information
  sunburst_ontology_reactionclass <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(sunburst_ontology_reactionclass) <- c("ids", "labels", "parents")

  #dataframe of EC Level 1 information
  level1 <- data.frame(
    "ids" = reactionClassesResults$class_ec_level_1$ecNumber,
    "labels" = paste((reactionClassesResults$class_ec_level_1$rxnClass),
                     paste("EC Number:", reactionClassesResults$class_ec_level_1$ecNumber),
                     sep = "\n"),
    "parents" = "",
    "hovertemplate" = paste((reactionClassesResults$class_ec_level_1$rxnClass),
                            paste("EC Number:", reactionClassesResults$class_ec_level_1$ecNumber),
                            paste(
                              reactionClassesResults$class_ec_level_1$metCount,
                              "input metabolites out of",
                              reactionClassesResults$class_ec_level_1$totalMetsInRxnClass ,
                              "total"
                            ),
                            paste(
                              reactionClassesResults$class_ec_level_1$proteinCount,
                              "input proteins out of",
                              reactionClassesResults$class_ec_level_1$totalProteinsInRxnClass ,
                              "total"
                            ),
                            paste(
                              reactionClassesResults$class_ec_level_1$reactionCount,
                              "reactions hit out of",
                              reactionClassesResults$class_ec_level_1$totalRxnsInClass ,
                              "total"
                            ),
                            sep = "\n"))


  #dataframe of EC Level 2 information
  EC_number_split <- unlist(strsplit(reactionClassesResults$class_ec_level_2$ecNumber, split = "\\."))

  level2 <- data.frame(
    "ids" = paste0(
      paste0(EC_number_split[seq(1, length(EC_number_split), 4)], ".-.-.-"), #Level 1 of EC_number
      "-", #split notation
      reactionClassesResults$class_ec_level_2$ecNumber), #Level 2 EC_number
    "labels" = paste(
      reactionClassesResults$class_ec_level_2$rxnClass,
      paste("EC Number:", reactionClassesResults$class_ec_level_2$ecNumber),
      sep = "\n"),
    "parents" = paste0(EC_number_split[seq(1, length(EC_number_split), 4)], ".-.-.-"), #Level 1 of EC_number
    "hovertemplate" = paste(
      reactionClassesResults$class_ec_level_2$rxnClass,
      paste("EC Number:", reactionClassesResults$class_ec_level_2$ecNumber),
      paste(
        reactionClassesResults$class_ec_level_2$metCount,
        "input metabolites out of",
        reactionClassesResults$class_ec_level_2$totalMetsInRxnClass ,
        "total"),
      paste(
        reactionClassesResults$class_ec_level_2$proteinCount,
        "input proteins out of",
        reactionClassesResults$class_ec_level_2$totalProteinsInRxnClass ,
        "total"),
      paste(
        reactionClassesResults$class_ec_level_2$reactionCount,
        "reactions hit out of",
        reactionClassesResults$class_ec_level_2$totalRxnsInClass ,
        "total"),
      sep = "\n")
  )

  EC_number_split <-unlist(strsplit(reactionClassesResults$class_ec_level_3$ecNumber, split = "\\."))

  level3 <- data.frame(
    "ids" = paste0(
      paste0(EC_number_split[seq(1, length(EC_number_split), 4)],
             ".",
             EC_number_split[seq(2, length(EC_number_split), 4)], ".-.-"), #Level 2 of EC_number
      "-",
      reactionClassesResults$class_ec_level_3$ecNumber), #Level 3 of EC_number
    "labels" = paste(
      reactionClassesResults$class_ec_level_3$rxnClass,
      paste("EC Number:", reactionClassesResults$class_ec_level_3$ecNumber),
      sep = "\n"),
    "parents" = paste0(
      paste0(EC_number_split[seq(1, length(EC_number_split), 4)], ".-.-.-"), #Level 1 of EC_number
      "-",
      paste0(EC_number_split[seq(1, length(EC_number_split), 4)], ".", EC_number_split[seq(2, length(EC_number_split), 4)], ".-.-")), #Level 2 of EC_number
    "hovertemplate" = paste(
      reactionClassesResults$class_ec_level_3$rxnClass,
      paste("EC Number:", reactionClassesResults$class_ec_level_3$ecNumber),
      paste(
        reactionClassesResults$class_ec_level_3$metCount,
        "input metabolites out of",
        reactionClassesResults$class_ec_level_3$totalMetsInRxnClass ,
        "total"),
      paste(
        reactionClassesResults$class_ec_level_3$proteinCount,
        "input proteins out of",
        reactionClassesResults$class_ec_level_3$totalProteinsInRxnClass ,
        "total"),
      paste(
        reactionClassesResults$class_ec_level_3$reactionCount,
        "reactions hit out of",
        reactionClassesResults$class_ec_level_3$totalRxnsInClass ,
        "total"),
      sep = "\n")
  )

  EC_number_split <- unlist(strsplit(reactionClassesResults$class_ec_level_4$ecNumber, split = "\\."))

  level4 <- data.frame(
    "ids" = paste0(
      paste0(
        EC_number_split[seq(1, length(EC_number_split), 4)],
        ".",
        EC_number_split[seq(2, length(EC_number_split), 4)],
        ".",
        EC_number_split[seq(3, length(EC_number_split), 4)],
        ".-"), #Level 3 of EC_number
      "-",
      reactionClassesResults$class_ec_level_4$ecNumber), #Level 4 of EC_number
    "labels" = paste(
      reactionClassesResults$class_ec_level_4$rxnClass,
      paste(reactionClassesResults$class_ec_level_4$ecNumber),
      paste(
        "EC Number:",
        "Metabolite Count:",
        reactionClassesResults$class_ec_level_4$metCount),
      sep = "\n"),
    "parents" = paste0(
      paste0(
        EC_number_split[seq(1, length(EC_number_split), 4)],
        ".",
        EC_number_split[seq(2, length(EC_number_split), 4)],
        ".-.-"), #Level 2 of EC_number
      "-",
      paste0(
        EC_number_split[seq(1, length(EC_number_split), 4)],
        ".",
        EC_number_split[seq(2, length(EC_number_split), 4)],
        ".",
        EC_number_split[seq(3, length(EC_number_split), 4)],
        ".-")), #Level 3 of EC_number
    "hovertemplate" = paste(
      reactionClassesResults$class_ec_level_4$rxnClass,
      paste("EC Number:",reactionClassesResults$class_ec_level_4$ecNumber),
      paste(
        reactionClassesResults$class_ec_level_4$metCount,
        "input metabolites out of",
        reactionClassesResults$class_ec_level_4$totalMetsInRxnClass ,
        "total"),
      paste(
        reactionClassesResults$class_ec_level_4$proteinCount,
        "input proteins out of",
        reactionClassesResults$class_ec_level_4$totalProteinsInRxnClass ,
        "total"),
      paste(
        reactionClassesResults$class_ec_level_4$reactionCount,
        "reactions hit out of",
        reactionClassesResults$class_ec_level_4$totalRxnsInClass ,
        "total"),
      sep = "\n")
  )

  sunburst_ontology_reactionclass <- rbind(level1, level2, level3, level4)

  colors_sunburst <- Polychrome::sortByLuminance(Polychrome::light.colors())[c(6, 7, 8, 9, 10, 11, 12)]

  for (i in 1:nrow(sunburst_ontology_reactionclass))
  {
    num_of_dashes <-
      length(which(
        strsplit(sunburst_ontology_reactionclass$labels[i], split = "\\.")[[1]] == "-"
      ))
    if (num_of_dashes == 3)
    {
      sunburst_ontology_reactionclass$color[i] <- adjustcolor(colors_sunburst[as.numeric(strsplit(sunburst_ontology_reactionclass$ids[i], split = "\\.")[[1]][1])], alpha.f = 0.8)
    }
    if (num_of_dashes == 2)
    {
      sunburst_ontology_reactionclass$color[i] <- adjustcolor(colors_sunburst[as.numeric(strsplit(sunburst_ontology_reactionclass$ids[i], split = "\\.")[[1]][1])], alpha.f = 0.6)
    }
    if (num_of_dashes == 1)
    {
      sunburst_ontology_reactionclass$color[i] <- adjustcolor(colors_sunburst[as.numeric(strsplit(sunburst_ontology_reactionclass$ids[i], split = "\\.")[[1]][1])], alpha.f = 0.4)
    }
    if (num_of_dashes == 0)
    {
      sunburst_ontology_reactionclass$color[i] <- adjustcolor(colors_sunburst[as.numeric(strsplit(sunburst_ontology_reactionclass$ids[i], split = "\\.")[[1]][1])], alpha.f = 0.2)
    }
  }
  return(sunburst_ontology_reactionclass)

}


#' Creates the input dataframe for the upset plot created in 'plotAnalyteOverlapPerRxnLevel'
#'
#' @param reactionsResults output of getReactionClassesForAnalytes()
#' @param includeCofactorMets whether or not to include metabolite cofactors (TRUE/FALSE)
#' @noRd
buildAnalyteOverlapPerRxnLevelUpsetDataframe <- function(reactionsResults, includeCofactorMets = FALSE) {
  if (includeCofactorMets == FALSE)
  {
    reactionsResults$met2rxn <- reactionsResults$met2rxn %>% dplyr::filter(.data$isCofactor == 0)
  }
  if(nrow(reactionsResults$met2rxn)>0)
  {
    met2rxn_EC <- reactionsResults$met2rxn %>% dplyr::filter(!dplyr::if_any("ecNumber", is.na))
    if(nrow(met2rxn_EC)>0)
    {
      EC_number_split_met <- unlist(strsplit(met2rxn_EC$ecNumber,split="\\."))
      input2reactions_mets <- cbind(
        c(met2rxn_EC$metSourceId),
        c(met2rxn_EC$ecNumber),
        c(paste0(EC_number_split_met[seq(1, length(EC_number_split_met), 4)]))
      )


    }

    met2rxn_NoEC <- reactionsResults$met2rxn %>% dplyr::filter(dplyr::if_any("ecNumber", is.na))
  }
  if(nrow(reactionsResults$prot2rxn)>0)
  {
    prot2rxn_EC <- reactionsResults$prot2rxn %>% dplyr::filter(!dplyr::if_any("ecNumber", is.na))
    if(nrow(prot2rxn_EC)>0)
    {
      EC_number_split_prot <- unlist(strsplit(prot2rxn_EC$ecNumber,split="\\."))
      input2reactions_prot <- cbind(
        c(prot2rxn_EC$uniprot),
        c(prot2rxn_EC$ecNumber),
        c(paste0(EC_number_split_prot[seq(1, length(EC_number_split_prot), 4)]))
      )
    }

    prot2rxn_NoEC <- reactionsResults$prot2rxn %>% dplyr::filter(dplyr::if_any("ecNumber", is.na))
  }

  if(exists("input2reactions_mets") && exists("input2reactions_prot"))
  {
    input2reactions <- as.data.frame(
      rbind(input2reactions_mets, input2reactions_prot))
  } else if (exists("input2reactions_mets"))
  {
    input2reactions <- as.data.frame(input2reactions_mets)
  } else if (exists("input2reactions_prot"))
  {
    input2reactions <- as.data.frame(input2reactions_prot)
  }

  input2reactions_list <- split(input2reactions$V1, input2reactions$V3)

  if (length(input2reactions_list) != 7)
  {
    seq <- 1:7
    missing_ecNum <- as.character(seq[!seq %in% as.numeric(names(input2reactions_list))])

    for(i in 1:length(missing_ecNum))
    {
      input2reactions_list[[missing_ecNum[i]]] <- vector()
    }

  }

  for (i in 1:7)
  {
    if(names(input2reactions_list)[i] == "1")
    {
      names(input2reactions_list)[i] = "Oxidoreductases: 1.-.-.- "
    }
    else if(names(input2reactions_list)[i] == "2")
    {
      names(input2reactions_list)[i] = "Transferases: 2.-.-.-"
    }
    else if(names(input2reactions_list)[i] == "3")
    {
      names(input2reactions_list)[i] = "Hydrolases: 3.-.-.-"
    }
    else if(names(input2reactions_list)[i] == "4")
    {
      names(input2reactions_list)[i] = "Lyases: 4.-.-.-"
    }
    else if(names(input2reactions_list)[i] == "5")
    {
      names(input2reactions_list)[i] = "Isomerases: 5.-.-.-"
    }
    else if(names(input2reactions_list)[i] == "6")
    {
      names(input2reactions_list)[i] = "Ligases: 6.-.-.-"
    }
    else if(names(input2reactions_list)[i] == "7")
    {
      names(input2reactions_list)[i] = "Translocases: 7.-.-.-"
    }
  }

  if(nrow(met2rxn_NoEC) > 0 && nrow(prot2rxn_NoEC) >0)
  {
    NoEC <-c(met2rxn_NoEC$metSourceId, prot2rxn_NoEC$uniprot)
  } else if (nrow(met2rxn_NoEC) > 0)
  {
    NoEC <- met2rxn_NoEC$metSourceId
  } else if (nrow(prot2rxn_NoEC) >0)
  {
    NoEC <- prot2rxn_NoEC$uniprot
  }

  input2reactions_list$"Non-Enzymatic" <- NoEC

  for (i in 1:length(input2reactions_list))
  {
    if (length(input2reactions_list[[i]]) == 0)
    {
      next
    } else
    {input2reactions_list[[i]] <- unlist(strsplit(input2reactions_list[[i]], "[|]"))}
  }

  return(input2reactions_list)
}



#' Creates the input dataframe for the interactive plot created in 'plotChemicalClass'
#'
#' @param chemicalClassResults output of getChemClass()
#' @noRd
buildChemicalClassDataframe <- function(chemicalClassResults) {

  chemicalClassResults_split <- split(chemicalClassResults$met_classes, chemicalClassResults$met_classes$source)

  sunburst_ontology_chemicalClass <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(sunburst_ontology_chemicalClass) <- c("ids", "labels", "parents")

  sunburst_ontology_chemicalClass <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(sunburst_ontology_chemicalClass) <- c("ids", "labels", "parents")

  if (length(chemicalClassResults_split) == 2)
  {
    hmdb_levels <- split(chemicalClassResults_split$hmdb, chemicalClassResults_split$hmdb$class_level_name)
    lipidmaps_levels <- split(chemicalClassResults_split$lipidmaps, chemicalClassResults_split$lipidmaps$class_level_name)

    missing_level2 <- which(is.na(merge(hmdb_levels$ClassyFire_super_class, hmdb_levels$ClassyFire_class, by = 1, all = TRUE)$class_name.y)==TRUE)

    hmdb_collate <- merge(hmdb_levels$ClassyFire_super_class, hmdb_levels$ClassyFire_class, by = 1)

    sunburst_ontology_chemicalClass <- rbind(sunburst_ontology_chemicalClass, data.frame("ids" = unique(hmdb_levels$ClassyFire_super_class$class_name), "labels" = unique(hmdb_levels$ClassyFire_super_class$class_name), "parents"= ""))

    if("Lipids and lipid-like molecules" %in% sunburst_ontology_chemicalClass$ids == FALSE)
    {
      sunburst_ontology_chemicalClass <- rbind(sunburst_ontology_chemicalClass, data.frame("ids" = "Lipids and lipid-like molecules", "labels" = "Lipids and lipid-like molecules", "parents"= ""))
    }

    sunburst_ontology_chemicalClass <- rbind(sunburst_ontology_chemicalClass, unique(data.frame("ids" = paste0(hmdb_collate$class_name.x,"-", hmdb_collate$class_name.y), "labels" = hmdb_collate$class_name.y, "parents"= hmdb_collate$class_name.x)))

    if(length(missing_level2)>0)
    {
      sunburst_ontology_chemicalClass <- rbind(sunburst_ontology_chemicalClass, unique(data.frame("ids" = paste0(hmdb_levels$ClassyFire_super_class[missing_level2,]$class_name,"-", hmdb_levels$ClassyFire_super_class[missing_level2,]$common_names), "labels" = paste0(hmdb_levels$ClassyFire_super_class[missing_level2,]$common_names, '\n', hmdb_levels$ClassyFire_super_class[missing_level2,][,1]), "parents"= hmdb_levels$ClassyFire_super_class[missing_level2,]$class_name)))
    }

    missing_level3 <- which(is.na(merge(hmdb_levels$ClassyFire_class, hmdb_levels$ClassyFire_sub_class, by = 1, all = TRUE)$class_name.y)==TRUE)

    if(length(missing_level3)>0)
    {
      sunburst_ontology_chemicalClass <- rbind(sunburst_ontology_chemicalClass, unique(data.frame("ids" = paste0(hmdb_collate[missing_level3,]$class_name.y,"-", hmdb_collate[missing_level3,]$common_names.x), "labels" = paste0(hmdb_collate[missing_level3,]$common_names.x, '\n', hmdb_collate[missing_level3,][,1]), "parents"= paste0(hmdb_collate[missing_level3,]$class_name.x,"-", hmdb_collate[missing_level3,]$class_name.y))))
    }

    hmdb_collate <- merge(hmdb_collate, hmdb_levels$ClassyFire_sub_class, by = 1)

    sunburst_ontology_chemicalClass <- rbind(sunburst_ontology_chemicalClass, unique(data.frame("ids" = paste0(hmdb_collate$class_name.y,"-", hmdb_collate$class_name), "labels" = hmdb_collate$class_name, "parents"= paste0(hmdb_collate$class_name.x,"-", hmdb_collate$class_name.y))))


    sunburst_ontology_chemicalClass <- rbind(sunburst_ontology_chemicalClass, unique(data.frame("ids" = paste0(hmdb_collate$class_name,"-", hmdb_collate$common_names.x), "labels" = paste0(hmdb_collate$common_names.x, '\n', hmdb_collate[,1]), "parents"= paste0(hmdb_collate$class_name.y,"-", hmdb_collate$class_name))))

    lipidmaps_levels$LipidMaps_category$class_name <- gsub(" \\[.*","",lipidmaps_levels$LipidMaps_category$class_name)
    lipidmaps_levels$LipidMaps_main_class$class_name <- gsub(" \\[.*","",lipidmaps_levels$LipidMaps_main_class$class_name)
    lipidmaps_levels$LipidMaps_sub_class$class_name <- gsub(" \\[.*","",lipidmaps_levels$LipidMaps_sub_class$class_name)

    if (any(unique(paste0("Lipids and lipid-like molecules-", lipidmaps_levels$LipidMaps_category$class_name)) %in% sunburst_ontology_chemicalClass$ids == FALSE))
    {
      index <- which(unique(paste0("Lipids and lipid-like molecules-", lipidmaps_levels$LipidMaps_category$class_name)) %in% sunburst_ontology_chemicalClass$ids == FALSE)

      sunburst_ontology_chemicalClass <- rbind(sunburst_ontology_chemicalClass, unique(data.frame("ids" = paste0("Lipids and lipid-like molecules-", lipidmaps_levels$LipidMaps_category$class_name), "labels" = lipidmaps_levels$LipidMaps_category$class_name, "parents"= "Lipids and lipid-like molecules"))[index,])
    }

    lipidmaps_collate <- merge(lipidmaps_levels$LipidMaps_category, lipidmaps_levels$LipidMaps_main_class, by = 1)
    lipidmaps_collate <- merge(lipidmaps_collate, lipidmaps_levels$LipidMaps_sub_class, by = 1)

    if (any(unique(paste0(lipidmaps_collate$class_name.x, "-", lipidmaps_collate$class_name.y)) %in% sunburst_ontology_chemicalClass$ids == FALSE))
    {
      index <- which(unique(paste0(lipidmaps_collate$class_name.x, "-", lipidmaps_collate$class_name.y)) %in% sunburst_ontology_chemicalClass$ids == FALSE)

      sunburst_ontology_chemicalClass <- rbind(sunburst_ontology_chemicalClass, unique(data.frame("ids" = paste0(lipidmaps_collate$class_name.x, "-", lipidmaps_collate$class_name.y), "labels" = lipidmaps_collate$class_name.y, "parents"= paste0("Lipids and lipid-like molecules-", lipidmaps_collate$class_name.x)))[index,])
    }

    sunburst_ontology_chemicalClass <- rbind(sunburst_ontology_chemicalClass, unique(data.frame("ids" = paste0(lipidmaps_collate$class_name.y, "-", lipidmaps_collate$class_name), "labels" = lipidmaps_collate$class_name, "parents"= paste0(lipidmaps_collate$class_name.x, "-", lipidmaps_collate$class_name.y))))

    sunburst_ontology_chemicalClass <- rbind(sunburst_ontology_chemicalClass, unique(data.frame("ids" = paste0(lipidmaps_collate$class_name, "-", lipidmaps_collate$common_names.x), "labels" = paste0(lipidmaps_collate$common_names.x, '\n', lipidmaps_collate[,1]), "parents"= paste0(lipidmaps_collate$class_name.y, "-", lipidmaps_collate$class_name))))

  }
  if (length(chemicalClassResults_split) == 1)
  {
    if(names(chemicalClassResults_split[1]) == "hmdb")
    {
      hmdb_levels <- split(chemicalClassResults_split$hmdb, chemicalClassResults_split$hmdb$class_level_name)

      missing_level2 <- which(is.na(merge(hmdb_levels$ClassyFire_super_class, hmdb_levels$ClassyFire_class, by = 1, all = TRUE)$class_name.y)==TRUE)

      hmdb_collate <- merge(hmdb_levels$ClassyFire_super_class, hmdb_levels$ClassyFire_class, by = 1)

      sunburst_ontology_chemicalClass <- rbind(sunburst_ontology_chemicalClass, data.frame("ids" = unique(hmdb_levels$ClassyFire_super_class$class_name), "labels" = unique(hmdb_levels$ClassyFire_super_class$class_name), "parents"= ""))

      sunburst_ontology_chemicalClass <- rbind(sunburst_ontology_chemicalClass, unique(data.frame("ids" = paste0(hmdb_collate$class_name.x,"-", hmdb_collate$class_name.y), "labels" = hmdb_collate$class_name.y, "parents"= hmdb_collate$class_name.x)))

      if(length(missing_level2)>0)
      {
        sunburst_ontology_chemicalClass <- rbind(sunburst_ontology_chemicalClass, unique(data.frame("ids" = paste0(hmdb_levels$ClassyFire_super_class[missing_level2,]$class_name,"-", hmdb_levels$ClassyFire_super_class[missing_level2,]$common_names), "labels" = paste0(hmdb_levels$ClassyFire_super_class[missing_level2,]$common_names, '\n', hmdb_levels$ClassyFire_super_class[missing_level2,][,1]), "parents"= hmdb_levels$ClassyFire_super_class[missing_level2,]$class_name)))
      }

      missing_level3 <- which(is.na(merge(hmdb_levels$ClassyFire_class, hmdb_levels$ClassyFire_sub_class, by = 1, all = TRUE)$class_name.y)==TRUE)

      if(length(missing_level3)>0)
      {
        sunburst_ontology_chemicalClass <- rbind(sunburst_ontology_chemicalClass, unique(data.frame("ids" = paste0(hmdb_collate[missing_level3,]$class_name.y,"-", hmdb_collate[missing_level3,]$common_names.x), "labels" = paste0(hmdb_collate[missing_level3,]$common_names.x, '\n', hmdb_collate[missing_level3,][,1]), "parents"= paste0(hmdb_collate[missing_level3,]$class_name.x,"-", hmdb_collate[missing_level3,]$class_name.y))))
      }

      hmdb_collate <- merge(hmdb_collate, hmdb_levels$ClassyFire_sub_class, by = 1)

      sunburst_ontology_chemicalClass <- rbind(sunburst_ontology_chemicalClass, unique(data.frame("ids" = paste0(hmdb_collate$class_name.y,"-", hmdb_collate$class_name), "labels" = hmdb_collate$class_name, "parents"= paste0(hmdb_collate$class_name.x,"-", hmdb_collate$class_name.y))))


      sunburst_ontology_chemicalClass <- rbind(sunburst_ontology_chemicalClass, unique(data.frame("ids" = paste0(hmdb_collate$class_name,"-", hmdb_collate$common_names.x), "labels" = paste0(hmdb_collate$common_names.x, '\n', hmdb_collate[,1]), "parents"= paste0(hmdb_collate$class_name.y,"-", hmdb_collate$class_name))))
    }

    else if (names(chemicalClassResults_split[1]) == "lipidmaps")
    {
      lipidmaps_levels <- split(chemicalClassResults_split$lipidmaps, chemicalClassResults_split$lipidmaps$class_level_name)

      sunburst_ontology_chemicalClass <- rbind(sunburst_ontology_chemicalClass, unique(data.frame("ids" = lipidmaps_levels$LipidMaps_category$class_name, "labels" = lipidmaps_levels$LipidMaps_category$class_name, "parents"= "")))

      lipidmaps_collate <- merge(lipidmaps_levels$LipidMaps_category, lipidmaps_levels$LipidMaps_main_class, by = 1)
      lipidmaps_collate <- merge(lipidmaps_collate, lipidmaps_levels$LipidMaps_sub_class, by = 1)

      sunburst_ontology_chemicalClass <- rbind(sunburst_ontology_chemicalClass, unique(data.frame("ids" = paste0(lipidmaps_collate$class_name.x, "-", lipidmaps_collate$class_name.y), "labels" = lipidmaps_collate$class_name.y, "parents"= lipidmaps_collate$class_name.x)))

      sunburst_ontology_chemicalClass <- rbind(sunburst_ontology_chemicalClass, unique(data.frame("ids" = paste0(lipidmaps_collate$class_name.y, "-", lipidmaps_collate$class_name), "labels" = lipidmaps_collate$class_name, "parents"= paste0(lipidmaps_collate$class_name.x, "-", lipidmaps_collate$class_name.y))))

      sunburst_ontology_chemicalClass <- rbind(sunburst_ontology_chemicalClass, unique(data.frame("ids" = paste0(lipidmaps_collate$class_name, "-", lipidmaps_collate$common_names.x), "labels" = paste0(lipidmaps_collate$common_names.x, '\n', lipidmaps_collate[,1]), "parents"= paste0(lipidmaps_collate$class_name.y, "-", lipidmaps_collate$class_name))))

    }
  }

  rownames(sunburst_ontology_chemicalClass) <- c(1:nrow(sunburst_ontology_chemicalClass))

  return(sunburst_ontology_chemicalClass)

}


#' Creates the input dataframe for the interactive plot created in 'plotCataNetwork'
#'
#' @param rampFastCataResults output of getChemClass()
#' @noRd
buildCataNetworkDataframe <- function(rampFastCataResults)
{
  #Set Edges
  myedges = rampFastCataResults[,c("query_relation", "input_common_name","rxn_partner_common_name", "Source")]
  colnames(myedges)[2:3] <-c("from","to")

  myedges$color <- ifelse(myedges$Source=="HMDB", "#ca1f7b", ifelse(myedges$Source=="Rhea", "#cc5500", "#008080"))
  myedges$highlight <- myedges$color

  #Set nodes
  mynodes=c(unique(myedges$from),unique(myedges$to))
  mycol=c(rep("black",length(unique(myedges$from))),
          rep("#d2e5f6",length(unique(myedges$to))))
  mynames <- mynodes

  mynodes <- data.frame(color=mycol,id=mynames,label=mynames)

  duplicates <- mynodes[mynodes[,3] %in% mynodes[duplicated(mynodes[3]),3],]

  if(nrow(duplicates)>0)
  {
    mynodes <- mynodes[-c(as.numeric(rownames(duplicates))),]

    duplicates <- split(duplicates, duplicates$id)

    for (i in 1:length(duplicates))
    {
      if(length(which(grepl("black", duplicates[[i]]$color)))>0)
      {
        single <- unique(duplicates[[i]][grepl("black", duplicates[[i]]$color),])
      } else
      {
        single <- duplicates[[i]][1,]
      }

      mynodes <- rbind(mynodes, single)
    }
  }

  for (i in 1:nrow(mynodes))
  {
    if (mynodes[i,1] == "black")
    {
      if (myedges[which(myedges$from == mynodes[i,2])[1], 1]=="met2gene" | myedges[which(myedges$from == mynodes[i,2])[1], 1]=="met2protein")
      {
        mynodes$shape[i] <- "dot"
      } else {
        mynodes$shape[i] <- "square"
      }
    } else if (mynodes[i,1] == "#d2e5f6")
    {
      if (myedges[which(myedges$to == mynodes[i,2])[1], 1]=="met2gene" | myedges[which(myedges$to == mynodes[i,2])[1], 1]=="met2protein")
      {
        mynodes$shape[i] <- "square"
      } else {
        mynodes$shape[i] <- "dot"
      }
    }
  }

  return(list("mynodes" = mynodes, "myedges" = myedges))
}
