#' Find all synonym from a given metabolite's name
#' This function is used to filter out some super common synonyms like glyceride
#' Now, this function only format the user input, so the user vector, dataframe,
#' and entire string separated by comma are working.
#' @param synonym name to search for
#' @param full bool if return whole data.frame
#' @param return_rampIds bool to return ramp Ids with output
#' (there are some common synonyms that will mess up whole searching)
#' @param db a RaMP database object
#' @return a data frame that contains synonym in the first column rampId in the second column
rampFindSynonymFromSynonym <- function( synonym,full = FALSE,
	return_rampIds = FALSE, db = RaMP()){
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
  list_metabolite <- sapply(list_metabolite,shQuote)
  list_metabolite <- paste(list_metabolite,collapse = ",")

  query <- paste0("select Synonym as origins, rampId from analytesynonym where Synonym in (",
                  list_metabolite,
                  ");")

  if(get("is_sqlite", pkg.globals)) {
    query <- paste0("select Synonym as origins, rampId from analytesynonym where Synonym COLLATE NOCASE in (",
                    list_metabolite,
                    ");")
  }

  df1 <- RaMP::runQuery(query, db)

  if(return_rampIds || nrow(df1) < 1) {
      return(df1)
  } else {
      rampid <- df1$rampId
      rampid <- sapply(rampid,shQuote)
      rampid <- paste(rampid,collapse = ",")
      query <- paste0("select * from analytesynonym where rampId in (",rampid,");")

      df2 <- RaMP::runQuery(query, db)

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
  list_id <- sapply(list_id,shQuote)
  list_id <- paste(list_id,collapse = ",")

  query <- paste0("select * from source where rampId in (",list_id,");")

  df <- RaMP::runQuery(query, db)

  if(full){
    return(df)
  } else{
    return(df[,1])
  }
}


#' Fast search given a list of metabolites source Id
#' @param sourceid a vector of synonym that need to be searched
#' @param find_synonym bool if find all synonyms or just return same synonym
#' @param db a RaMP database object
#' @return a list contains all metabolits as name and pathway inside.
rampFastPathFromSource<- function( sourceid, find_synonym = FALSE, db = RaMP()){
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

  df1 <- RaMP::runQuery(query1, db)

  colnames(df1)[1] <-"sourceId2"
  #return(df1)
  rampid <- df1$rampId
  rampid <- sapply(rampid,shQuote)
  rampid <- paste(rampid,collapse = ",")
  query2 <- paste0("select * from analytehaspathway where
                   rampId in (",rampid,");")

  df2 <- RaMP::runQuery(query2, db)

  #return(df2)
  id_list <- unique(df2$pathwayRampId)
  id_list <- sapply(id_list,shQuote)
  id_list <- paste(id_list,collapse = ",")
  print(id_list)
  query3 <- paste0("select * from pathway where pathwayRampId in (",
                   id_list,");")

  df3 <- RaMP::runQuery(query3, db)

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
  list_metabolite <- sapply(list_metabolite,shQuote)
  list_metabolite <- paste(list_metabolite,collapse = ",")

  query <- paste0("select sourceId,IDtype as analytesource, rampId from source where sourceId in (",list_metabolite,");")
  df <- RaMP::runQuery(query, db)

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
#' @param inferIdMapping if FALSE, the survey only reports on class annotations made directly on the input ids.
#' If inferIdMapping is set to TRUE, the ids are cross-referenced or mapped to related ids that contain metabolite class annotations.
#' The default is FALSE. Note that this utility method is typcally used within chemical class enrichment and is passed the value of this parameter.
#' @returns a list object with two keys, 'mets' and 'pop' that each has a table of metabolite or population
#' chemical classes and metabolite counts per class. This supports the chemicalClassEnrichment function.
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
  # two routes depending on wether we have a user provide population or the all-DB population
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
bhCorrect <- function(resultMat) {
  resultMat <- resultMat[order(resultMat$`p-value`),]
  bhPvals <- stats::p.adjust(resultMat$`p-value`, method = "BH")
  resultMat$adjP_BH <- bhPvals
  return(resultMat)
}

#' Get class info for an input of metabolite source Ids
#' @param sourceIds a vector of analytes (genes or metabolites) that need to be searched
#' @param db a RaMP database object
#' @return a dataframe of chemClass info
rampFindClassInfoFromSourceId<-function(sourceIds, db = RaMP()){
    sourceIds <- unique(sourceIds)
    checkIdPrefixes(sourceIds)
    idsToCheck <- sapply(sourceIds,function(x){
        if(!grepl("hmdb|chebi|LIPIDMAPS",x)){
            return(x)
        }
    })
    idsToCheck <- paste(idsToCheck, collapse = "','")
    idsToCheck <- paste("'" ,idsToCheck, "'", sep = "")

    sql <- paste("select * from source where sourceId in (",idsToCheck,")")

    potentialMultiMappings <- RaMP::runQuery(sql, db)

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

    sql <- paste("select distinct a.ramp_id, b.sourceId, a.class_level_name, a.class_name, a.source from metabolite_class a, source b
          where b.rampId = a.ramp_id and b.sourceId in (",metStr,")")

    metsData <- RaMP::runQuery(sql, db)

    metsData <- subset(metsData, "sourceId" %in% sourceIds)

    return(metsData)
}

#' Internal function for extracting annotations, used by pathway and chemical enrichment test functions
#' @param analytes a vector of analytes (genes or metabolites) that need to be searched
#' @param PathOrChem return "path" information for pathways or "chem" for chemical class
#' @param namesOrIds whether input is "names" or "ids" (default is "ids")
#' @param db a RaMP database object
#' @return a list of rampIds for "path" or a dataframe of chemClass info
getRaMPInfoFromAnalytes<-function( analytes,
                                  namesOrIds = "ids",
                                  PathOrChem = "path", db = RaMP()){
    if(PathOrChem == "path"){
        if(namesOrIds == "names"){
            synonym <- rampFindSynonymFromSynonym(synonym=analytes,
                                                  return_rampIds=FALSE)

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
##' @param pathway_definitions If "RaMP" (default), use pathway definitions within RaMP-DB. Else, supply path to gmx file containing custom pathway definitions. GMX files are a tab-separated format that contain one analyte set per column, with the name of the set in the first row, and constituent analytes in subsequent rows
##' @param analyte_type "genes" or "metabolites"
##' @param db a RaMP database object
##' @return dataframe of all analytes that map to the input pathways
##' @author Andrew Christopher Patt
#' @importFrom utils head
buildFrequencyTables<-function( inputdf, pathway_definitions="RaMP", analyte_type, db = RaMP()) {

  if(pathway_definitions == "RaMP") {
	print("Now in buildFrequencyTables function")
	print(mode(inputdf))
	print(class(inputdf))
	head(inputdf$pathwayRampId)

    ## Get pathway ids that contain the user analytes
    pid <- unique(inputdf$pathwayRampId);
    list_pid <- sapply(pid,shQuote)
    list_pid <- paste(list_pid,collapse = ",")

    ## Retrieve compound ids associated with background pathways and count
    query <- paste0("select * from analytehaspathway where pathwayRampId in (",
                    list_pid,")")

    input_RampIds <- RaMP::runQuery(query, db)

    return(input_RampIds)
  } else {
    tryCatch(
      {
        if (analyte_type == "metabolites") {
          pathway_definitions <- readxl::read_excel(pathways, sheet = 1)
        } else if (analyte_type == "genes") {
          pathway_definitions <- readxl::read_excel(pathways, sheet = 2)
        }
      },
      error = function(e) {
        print("Pathway file could not be found or is improperly formatted. Please supply path to GMX file for custom pathway definitions")
      }
    )
    input_RampIds <- data.frame(rampId=character(),
                                pathwayRampId=character())
    pid <- unique(inputdf$pathwayRampId)

    for(i in pid){
      temp <- data.frame(rampId = pathway_definitions[,i][which(!is.na(pathway_definitions[,i])),],
                         pathwayRampId = i)
      colnames(temp) <- c("rampId","pathwayRampId")
      input_RampIds <- rbind(input_RampIds,temp)
    }
    input_RampIds$pathwaySource = "custom"
    if(analyte_type == "metabolites"){
      input_RampIds$rampId = paste0("RAMP_C",input_RampIds$rampId)
    }else if(analyte_type == "genes"){
      input_RampIds$rampId = paste0("RAMP_G",input_RampIds$rampId)
    }
    return(input_RampIds)
  }
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

  input_RampId_G <- input_RampIds[grep("RAMP_G", input_RampIds$rampId), ]
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
find_duplicate_pathways <- function(db = RaMP()){

  .Deprecated("findDuplicatPathways")

  pathway_overlap = analyte_result
  duplicate_pairs = data.frame(Pathway1=character(),Pathway2=character())
  for(i in 1:ncol(pathway_overlap)){
    duplicates <- which(pathway_overlap[,i]==1)
    duplicates <- duplicates[-which(duplicates==i)]
    if(length(duplicates!=0)){
      duplicate_pairs <- rbind(duplicate_pairs,
                               data.frame(Pathway1=colnames(pathway_overlap)[duplicates],
                                          Pathway2=colnames(pathway_overlap)[i]))
    }
  }
  query <- "select * from analytehaspathway where pathwaySource != 'hmdb';"

  allpids <- RaMP::runQuery(query, db)

  duplicate_pathways <- apply(duplicate_pairs, 1, function(x){
    path1 <- x[1]
    path2 <- x[2]
    path1_source <- unique(allpids[which(allpids$pathwayRampId==path1),"pathwaySource"])
    path2_source <- unique(allpids[which(allpids$pathwayRampId==path2),"pathwaySource"])
    return(data.frame(path1=path1, path2=path2, path1_source=path1_source, path2_source=path2_source))
  })
  duplicate_pathways = do.call(rbind,duplicate_pathways)
  # If one duplicate is wiki, return reactome. Else, return pathway two
  duplicate_pathways <- apply(duplicate_pathways,1,function(x){
    if(x[3]=="wiki" & x[4]=="reactome"){
      return(x[1])
    }else if(x[3]=="reactome" & x[4]=="wiki"){
      return(x[2])
    }else{
      return(x[2])
    }
  })
  names(duplicate_pathways)=NULL
  return(duplicate_pathways)
}

##' Return list of duplicate Wikipathway IDs from Reactome. This may be unnecessary in the future
##' @return List of duplicate Wikipathway IDs from Reactome.
##' @param db a RaMP database object
##' @author John Braisted
findDuplicatePathways <- function(db = RaMP()) {

  query <- "select pathwayRampId from pathway where type = 'reactome';"

  reactomePIDs <- RaMP::runQuery(query, db)

  ar <- db@dbSummaryObjCache$analyte_result
  diag(ar) <- 0.0
  ar[ar != 1.0] <- 0.0
  colHits <- colnames(ar)[colSums(ar) >= 1.0]
  rowHits <- colnames(ar)[rowSums(ar) >= 1.0]
  ar2 <- ar[rowHits, colHits]
  n = 0

  for(r in rownames(ar2)) {
    colHits <- colnames(ar2)[ar2[r,]==1.0]
    rowHits <- rep(r, length(colHits))
    df <- data.frame(colHits)
    df <- cbind(df, rowHits)
    if(n == 0) {
      df2 <- df
    } else {
      df2 <- rbind(df2, df)
    }
    n = n + 1
  }

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
#' @param fishers_df The data frame generated by runFisherTest
#' @param pval_type Specifies which p-value to use as the filter threshold.
#' Permitted values are 'pval' and 'fdr' for chemical class and pathway enrichment.
#' Pathway enrichment also includes an optional 'holm' value for holm p-value corrections. Default is 'fdr'.
#' @param pval_cutoff return pathways where pval_type p-values are < pval_cutoff
#' @return list:[[1]]Dataframe with pathway enrichment results, only significant pathways
#' [[2]]analyte type
#' @examples
#' \dontrun{
#' analyteList <- c("MDM2", "TP53", "glutamate", "creatinine")
#'
#' fisher.results <- runCombinedFisherTest(analytes = analyteList, namesOrIds = 'names')
#' filtered.fisher.results <- FilterFishersResults(fisher.results, pval_type='fdr', pval_cutoff = 0.10)
#' }
#' @export
FilterFishersResults <- function(fishers_df, pval_type = 'fdr', pval_cutoff = 0.1) {

  print("Filtering Fisher Results...")

  # Check to see whether the output is from ORA performed on genes and metabolites
  # or genes or metabolites
  result_type <- fishers_df$result_type

  if(result_type == 'pathway_enrichment') {
    analyte_type <- fishers_df$analyte_type
    fishers_df <- fishers_df$fishresults

    print("Fisher Result Type: Pathway Enrichment")

    if (analyte_type != 'both') {
      if (pval_type == 'holm') {
        return(list(fishresults = fishers_df[which(fishers_df[, "Pval_Holm"] <=
                                                     pval_cutoff), ], analyte_type = analyte_type))
      } else if (pval_type == 'fdr') {
        return(list(fishresults = fishers_df[which(fishers_df[, "Pval_FDR"] <=
                                                     pval_cutoff), ], analyte_type = analyte_type))
      } else if (pval_type == 'pval') {
        return(list(fishresults = fishers_df[which(fishers_df[, "Pval"] <=
                                                     pval_cutoff), ], analyte_type = analyte_type))
      } else {
        warning(paste0("The pval_type parameter should be one of three values, 'fdr', 'holm' or 'pval', entered value pval_type= ", pval_type))
        return(NULL)
      }
    } else { # ORA was performed on both genes and metabolites:
      if (pval_type == 'holm') {
        return(list(fishresults = fishers_df[which(fishers_df[, "Pval_combined_Holm"] <=
                                                     pval_cutoff), ], analyte_type = analyte_type))
      } else if (pval_type == 'fdr') {
        return(list(fishresults = fishers_df[which(fishers_df[, "Pval_combined_FDR"] <=
                                                     pval_cutoff), ], analyte_type = analyte_type))
      } else if (pval_type == 'pval') {
        return(list(fishresults = fishers_df[which(fishers_df[, "Pval_combined"] <=
                                                     pval_cutoff), ], analyte_type = analyte_type))
      } else {
        warning(paste0("The pval_type parameter should be one of three values, 'fdr', 'holm' or 'pval', entered value pval_type= ", pval_type))
        return(NULL)
      }
    }
  } else if(result_type == "chemical_class_enrichment") {

    print("Fisher Result Type: Chemical Class Enrichstrment")

    if(pval_type == 'pval') {
      criteriaCol <- 'p-value'
    } else if (pval_type == 'fdr') {
      criteriaCol <- 'adjP_BH'
    } else {
      warning(paste0("The pval_type parameter should be one of three values, 'fdr' or 'pval' for chemical class enrichment, entered value pval_type= ", pval_type))
      return(NULL)
    }

    for(result in names(fishers_df)) {

      #if(class(fishers_df[[result]]) == 'data.frame') {
       if(methods::is(fishers_df[[result]], 'data.frame')) {
        print(result)
        resultDf <- fishers_df[[result]]
        resultDf <- subset(resultDf, resultDf[[criteriaCol]] <= pval_cutoff)
        fishers_df[[result]] <- resultDf
      }
    }
    return(fishers_df)
  }
}


#'chemicalClassSurveyRampIdsConn is a helper function that takes a list of metabolite ids, a list of 'population' metabolite ids
#' and a MariaDB Connection object. The method returns metabolite class information for the metabolite list and a population of all ramp metabolites.
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
chemicalClassSurveyRampIdsConn <- function( mets, pop, inferIdMapping=TRUE, db = RaMP()) {

  mets <- unique(mets)

  checkIdPrefixes(mets)

  pop <- unique(pop)

  checkIdPrefixes(pop)

  result <- list()

  # first handle metabolites of interest
  metStr <- paste(mets, collapse = "','")
  metStr <- paste("'" ,metStr, "'", sep = "")


  isSQLite <- .is_sqlite(db)

  # if inferring ID mapping, the query goes through the source table to map input id to ramp id, then map to related ids having chem class annotations
  # if not ID mapping, then the match is directly on the input source ids. HMDB and LipidMaps IDs are supported directly, May 2023.
  if(inferIdMapping) {
    sql <- paste("select distinct a.ramp_id, b.sourceId, group_concat(distinct b.commonName order by b.commonName asc separator '; ') as common_names,
                   a.class_level_name, a.class_name, a.source as source, count(distinct(a.class_source_id)) as directIdClassHits from metabolite_class a, source b
                   where b.rampId = a.ramp_id and b.sourceId in (",metStr,")
                   group by a.ramp_Id, b.sourceId, a.class_level_name, a.class_name, a.source")

    if(isSQLite) {
      sql <- paste("select distinct a.ramp_id, b.sourceId, group_concat(distinct b.commonName COLLATE NOCASE) as common_names,
                  a.class_level_name, a.class_name, a.source as source, count(distinct(a.class_source_id)) as directIdClassHits from metabolite_class a, source b
                  where b.rampId = a.ramp_id and b.sourceId in (",metStr,")
                  group by a.ramp_Id, b.sourceId, a.class_level_name, a.class_name, a.source")
    }
  } else {
    sql = paste("select distinct c.ramp_id, c.class_source_id, group_concat(distinct s.commonName order by s.commonName asc separator '; ') as common_names,
                 c.class_level_name, c.class_name, c.source as source, count(distinct(c.class_source_id)) as directIdClassHits
                 from metabolite_class c, source s
                 where c.class_source_id in (",metStr,") and s.sourceId = c.class_source_id
                 group by c.class_source_id, c.class_level_name, c.class_name, c.source, c.ramp_id")

    if(isSQLite) {
      sql = paste("select distinct c.ramp_id, c.class_source_id, group_concat(distinct s.commonName COLLATE NOCASE) as common_names,
                 c.class_level_name, c.class_name, c.source as source, count(distinct(c.class_source_id)) as directIdClassHits
                 from metabolite_class c, source s
                 where c.class_source_id in (",metStr,") and s.sourceId = c.class_source_id
                 group by c.class_source_id, c.class_level_name, c.class_name, c.source, c.ramp_id")
    }

  }

  metsData <- RaMP::runQuery(sql, db)

  # get query summary
  metQueryReport <- queryReport(mets, metsData$sourceId)

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

  # Population info
  popStr <- paste(pop, collapse = "','")
  popStr <- paste("'" ,popStr, "'", sep = "")

  # a similar query on population ids, id mapping matches on mapped source ids, no id mapping matches input ids directly on annotated ids
  if(inferIdMapping) {
    sql <- paste("select distinct a.ramp_id, b.sourceId, a.class_level_name, a.class_name, a.source,
                  count(distinct(a.class_source_id)) as directIdClassHits
                  from metabolite_class a, source b
                  where b.rampId = a.ramp_id and b.sourceId in (",popStr,")
                 group by a.ramp_Id, b.sourceId, a.class_level_name, a.class_name, a.source")
  } else {
    sql <- paste("select distinct c.ramp_id, c.class_source_id, c.class_level_name, c.class_name, c.source,
                 count(distinct(c.class_source_id)) as directIdClassHits
                 from metabolite_class c
                 where c.class_source_id in (",popStr,")
                 group by c.class_source_id, c.class_level_name, c.class_name, c.source, c.ramp_id")
  }

  # ("select distinct c.ramp_id, c.class_source_id, group_concat(distinct s.commonName order by s.commonName asc separator '; ') as common_names,
  # c.class_level_name, c.class_name, c.source as source, count(distinct(c.class_source_id)) as directIdClassHits
  # from metabolite_class c, source s
  # where c.class_source_id in (",metStr,") and s.sourceId = c.class_source_id
  # group by c.class_source_id, c.class_level_name, c.class_name")

  popData <- RaMP::runQuery(sql, db)

  if(inferIdMapping) {
    popData <- subset(popData, "sourceId" %in% pop)
  } else {
    popData <- subset(popData, "class_source_id" %in% pop)
  }

  #need to filter for our source ids
  # popData <- subset(popData, "sourceId" %in% pop)

  # get query summary
  popQueryReport <- queryReport(pop, popData$sourceId)


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


#'chemicalClassSurveyRampIdsFullPopConn2 is a helper function that takes a list of metabolite ids and a MariaDB Connection object
#'and returns metabolite class information for the metabolite list and a population of all ramp metabolites.
#' @param mets a list object of prefixed metabolite ids of interest
#' @param inferIdMapping if FALSE, the survey only reports on class annotations made directly on the input ids.
#' If inferIdMapping is set to TRUE, the ids are cross-referenced or mapped to related ids that contain metabolite class annotations.
#' The default is TRUE.
#' @param db a RaMP database object
#' @returns a list object containing three objects 'count_summary', 'met_classes' and 'met_query_report'.
#' The count_summary is a dataframe containing metabolite classes and number of metabolites in each class.
#' The met_classes is a detailed listing of compound classes associated with each input metabolite
#' The met_query_report indicates the number of input metabolites, how many were found in the DB and the list of metabolites not found in RaMP DB.
chemicalClassSurveyRampIdsFullPopConn <- function( mets, inferIdMapping=TRUE, db = RaMP()) {

  mets <- unique(mets)

  checkIdPrefixes(mets)

  result <- list()

  # first handle metabolites of interest
  metStr <- paste(mets, collapse = "','")
  metStr <- paste("'" ,metStr, "'", sep = "")

  isSQLite = .is_sqlite(db)

  # Id mapping matches on source ids mapped via ramp ids in the source table. No id mapping matches on input ids directly.
  if(inferIdMapping) {
    sql <- paste("select distinct a.ramp_id, b.sourceId, group_concat(distinct b.commonName order by b.commonName asc separator '; ') as common_names,
     a.class_level_name, a.class_name, a.source as source, count(distinct(a.class_source_id)) as directIdClassHits from metabolite_class a, source b
          where b.rampId = a.ramp_id and b.sourceId in (",metStr,")
               group by a.ramp_Id, b.sourceId, a.class_level_name, a.class_name, a.source")
    if(isSQLite) {
      sql <- paste("select distinct a.ramp_id, b.sourceId, group_concat(distinct(b.commonName) COLLATE NOCASE) as common_names,
          a.class_level_name, a.class_name, a.source as source, count(distinct(a.class_source_id)) as directIdClassHits from metabolite_class a, source b
          where b.rampId = a.ramp_id and b.sourceId in (",metStr,")
          group by a.ramp_Id, b.sourceId, a.class_level_name, a.class_name, a.source")
    }
  } else {
    sql = paste("select distinct c.ramp_id, c.class_source_id, group_concat(distinct s.commonName order by s.commonName asc separator '; ') as common_names,
               c.class_level_name, c.class_name, c.source, count(distinct(c.class_source_id)) as directIdClassHits from metabolite_class c, source s
               where c.class_source_id in (",metStr,") and s.sourceId = c.class_source_id group by c.class_source_id, c.class_level_name, c.class_name")

    if(isSQLite) {
      sql = paste("select distinct c.ramp_id, c.class_source_id, group_concat(distinct(s.commonName) COLLATE NOCASE) as common_names,
                  c.class_level_name, c.class_name, c.source, count(distinct(c.class_source_id)) as directIdClassHits from metabolite_class c, source s
                  where c.class_source_id in (",metStr,") and s.sourceId = c.class_source_id group by c.class_source_id, c.class_level_name, c.class_name")
    }
  }

  metsData <- RaMP::runQuery(sql, db)

  # get query summary
  metQueryReport <- queryReport(mets, metsData$sourceId)

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
  sql <- paste("select class_level_name, class_name, count(1) as pop_hits from metabolite_class
                 group by class_level_name, class_name")

  popCountData <- RaMP::runQuery(sql, db)

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
#' @param analytes list of analytes (can be names or ids)
#'
#' @return comma separated list of single quoted analyte ids or names
#'
listToQueryString <- function(analytes) {
  analyteStr <- paste0("'", paste0(analytes, collapse = "','"), "'", sep="")
  return (analyteStr)
}


#' filterPathwaysByAnalyteCount utility method filtered a dataframe based on the number of analytes associated with rampPathwayIds contained in the dataframe.
#' Like fisher exact code, this one retains pathways with analyte count >= minPathwaySize, and having analyte_count < max_path_size
#'
#' @param pathway_dataframe a dataframe containing at least one column that contains rampPathwayIds
#' @param pathway_ramp_id_col_name the column name containing the rampPathwayIds
#' @param minPathwaySize the minimum number of pathway members (genes and metabolites) to include the pathway in the output (default = 5)
#' @param maxPathwaySize the maximum number of pathway memnbers (genes and metaboltes) to include the pathway in the output (default = 150)
#' @param db a RaMP databse object
filterPathwaysByAnalyteCount <- function( pathway_dataframe, pathway_ramp_id_col_name = 'pathwayRampId', minPathwaySize = 5, maxPathwaySize = 150, db = RaMP()) {
  pwIds <- unlist(pathway_dataframe[[pathway_ramp_id_col_name]])
  pwIdsStr <- listToQueryString(pwIds)

  sql <- paste0("select pathwayRampId, count(distinct(rampId)) as analyte_count from analytehaspathway where pathwayRampId in (", pwIdsStr,") group by pathwayRampId")

  res <- RaMP::runQuery(sql, db=db)
  res <- res[res$analyte_count >= minPathwaySize & res$analyte_count < maxPathwaySize,]
  keeperPW <- unlist(res$pathwayRampId)
  pathway_dataframe <- pathway_dataframe[pathway_dataframe[[pathway_ramp_id_col_name]] %in% keeperPW, ]
  return(pathway_dataframe)
}


#' Creates the input dataframe for the sunburst plot created in 'plotReactionClasses'
#'
#' @param reactionClassesResults output of getReactionClassesForAnalytes()
#' @importFrom grDevices adjustcolor
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

buildAnalyteOverlapPerRxnLevelUpsetDataframe <- function(reactionsResults, includeCofactorMets = FALSE) {

  if(nrow(reactionsResults$met2rxn)>0)
  {
    reactionsResults$met2rxn <- reactionsResults$met2rxn %>% dplyr::filter(!dplyr::if_any(.data$ecNumber, is.na))
    if(nrow(reactionsResults$met2rxn)>0)
    {
      EC_number_split_met <- unlist(strsplit(reactionsResults$met2rxn$ecNumber,split="\\."))
      input2reactions_mets <- cbind(
        c(reactionsResults$met2rxn$metSourceId),
        c(reactionsResults$met2rxn$ecNumber),
        c(paste0(EC_number_split_met[seq(1, length(EC_number_split_met), 4)]))
      )
    }
  }
  if (includeCofactorMets == FALSE)
  {
    reactionsResults$met2rxn <- reactionsResults$met2rxn %>% dplyr::filter(.data$isCofactor == 0)
  }
  if(nrow(reactionsResults$prot2rxn)>0)
  {
    reactionsResults$prot2rxn <- reactionsResults$prot2rxn %>% dplyr::filter(!dplyr::if_any(.data$ecNumber, is.na))
    if(nrow(reactionsResults$prot2rxn)>0)
    {
      EC_number_split_prot <- unlist(strsplit(reactionsResults$prot2rxn$ecNumber,split="\\."))
      input2reactions_prot <- cbind(
        c(reactionsResults$prot2rxn$uniprot),
        c(reactionsResults$prot2rxn$ecNumber),
        c(paste0(EC_number_split_prot[seq(1, length(EC_number_split_prot), 4)]))
      )
    }
  }

  if(exists("input2reactions_mets") && exists("input2reactions_prot"))
  {
    input2reactions <- as.data.frame(
      rbind(input2reactions_mets, input2reactions_prot))
  }
  else if (exists("input2reactions_mets"))
  {
    input2reactions <- as.data.frame(input2reactions_mets)
  }
  else if (exists("input2reactions_prot"))
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
      input2reactions_list <- c(input2reactions_list, i = list(NULL))
      names(input2reactions_list)[(length(input2reactions_list))] <- missing_ecNum[i]
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

  return(input2reactions_list)
}



#' Creates the input dataframe for the interactive plot created in 'plotChemicalClassSurvery'
#'
#' @param chemicalClassSurveryResults output of getReactionClassesForAnalytes()

buildChemicalClassSurveryDataframe <- function(chemicalClassSurveryResults) {

  chemicalClassSurveryResults_split <- split(chemicalClassSurveryResults$met_classes, chemicalClassSurveryResults$met_classes$source)

  sunburst_ontology_chemicalClass <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(sunburst_ontology_chemicalClass) <- c("ids", "labels", "parents")

  sunburst_ontology_chemicalClass <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(sunburst_ontology_chemicalClass) <- c("ids", "labels", "parents")

  if (length(chemicalClassSurveryResults_split) == 2)
  {
    hmdb_levels <- split(chemicalClassSurveryResults_split$hmdb, chemicalClassSurveryResults_split$hmdb$class_level_name)
    lipidmaps_levels <- split(chemicalClassSurveryResults_split$lipidmaps, chemicalClassSurveryResults_split$lipidmaps$class_level_name)

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
  if (length(chemicalClassSurveryResults_split) == 1)
  {
    if(names(chemicalClassSurveryResults_split[1]) == "hmdb")
    {
      hmdb_levels <- split(chemicalClassSurveryResults_split$hmdb, chemicalClassSurveryResults_split$hmdb$class_level_name)

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

    else if (names(chemicalClassSurveryResults_split[1]) == "lipidmaps")
    {
      lipidmaps_levels <- split(chemicalClassSurveryResults_split$lipidmaps, chemicalClassSurveryResults_split$lipidmaps$class_level_name)

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
