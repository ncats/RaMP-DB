# queries to retrieve and analyze chemical classes

#' Returns chemical class information comparing a metabolite subset to a larger metabolite population.
#'
#' @param mets a list object of source prepended metaboite ids, representing a metabolite set of interest
#' @param background an optional list of source prepended metaboite ids to be used as the background reference of
#' metabolites for enrichment. The background can be either a list of ids, a file name containing the id list,
#' one id per column (no file header row) or a specificed biospecimen type (available biospecimen types: "Blood",
#' "Adipose tissue", "Heart", "Urine", "Brain", "Liver", "Kidney","Saliva", or "Feces").
#' @param background_type one of 'database' (all analytes in the RaMP Database), 'list' (a list of input ids),
#' or 'file' in which case the background parameter will be a file path, or 'biospecimen' where the specified background parameter is
#' a RaMP HMDB metabolite ontology term (see background parameter, above, for the most common biospecimen background values).
#' @param includeRaMPids include internal RaMP identifiers (default is "FALSE")
#' @return Returns chemcial class information data including class count tallies and comparisons between metabolites of interest and the metabolite population,
#' metabolite mappings to classes, and query summary report indicating the number of input metabolites that were resolved and listing those metabolite ids
#' that are not found in the database.
#'
#' The returned object (return_obj below) contains three or four main result areas. Use str(return_obj) to see the structure described here.
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
#' # The background population can be a separate ID list (preferred) or all database entries
#' # (skip pop parameter).
#' pkg.globals <- setConnectionToRaMP(dbname="ramp2",username="root",conpass="",host = "localhost")
#' metClassResult <- chemicalClassSurvey(mets = metList)
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
chemicalClassSurvey <- function(mets, background = "database", background_type="database", includeRaMPids = FALSE){

  
  conn <- connectToRaMP()
  print("Starting Chemical Class Survey")

  if(background_type == "file") {
    bkgrnd <- utils::read.table(background, header=F)[,1]

    filteredMets <- mets[mets %in% bkgrnd]
    print(paste0("Number of input query ids: ",length(mets)))
    print(paste0("Number of input query ids found in supplied biospecimen background: ",length(filteredMets)))
    print(paste0("Excluded input query ids that are NOT found in supplied biospecimen background: ",(length(mets)-length(filteredMets))))

    if((length(filteredMets)) < 1) {
      warning("All input query mets were not found in the specified background ID list.")
      return(list())
    }

    mets <- filteredMets

  } else if(background_type == "list") {
    bkgrnd = background

    filteredMets <- mets[mets %in% bkgrnd]
    print(paste0("Number of input query ids: ",length(mets)))
    print(paste0("Number of input query ids found in supplied biospecimen background: ",length(filteredMets)))
    print(paste0("Excluded input query ids that are NOT found in supplied biospecimen background: ",(length(mets)-length(filteredMets))))

    if((length(filteredMets)) < 1) {
      warning("All input query mets were not found in the specified background ID list.")
      return(list())
    }

    mets <- filteredMets

  } else if(background_type == "database") {
    # use the full database as background
    bkgrnd = 'database'
  } else if (background_type == "biospecimen") {

    print(paste0("Biospecimen background specified: ", background))

    query <- paste0("select s.rampId, s.sourceId from source s, analytehasontology ao, ontology o
    where o.commonName in ('", background, "') and o.rampOntologyId=ao.rampOntologyId and s.rampId = ao.rampCompoundId")
    con <- connectToRaMP()
    bg <- RMariaDB::dbGetQuery(con, query)
    RMariaDB::dbDisconnect(con)

    # cases check if bg is empty (suggest to query for biospecimen types in ramp)
    if(is.null(bg) || nrow(bg) == 0) {
      warning("The input biospecimen type is not represented in the RaMP Database. Use RaMP::getOntologies() to see available terms.")
      return(list())
    }
    # check that all input query ids are inside the biospecimen type list
    # make sure that query ids are culled for ids outside the specified biospeciment background
    filteredMets <- mets[mets %in% bg$sourceId]
    print(paste0("Number of input query ids: ",length(mets)))
    print(paste0("Number of input query ids found in supplied biospecimen background: ",length(filteredMets)))
    print(paste0("Excluded input query ids that are NOT found in supplied biospecimen background: ",(length(mets)-length(filteredMets))))

    # check that we still have some query metabolites
    if((length(filteredMets)) < 1) {
      warning(paste0("All input query mets were not found in the specified biospecimen background: ",background))
      return(list())
    }

    # try to make a set of source ids assoicated with unique rampids
    bg <- bg[!duplicated(bg$rampId),]

    # use unique list of source ids for the biospecimen type
    # note that rampId counts will be used downstream for statistics, reducing redundancy
    bkgrnd <- unique(unlist(bg$sourceId))

    print(paste0("Number of metabolites in biospecimen background: ", length(bkgrnd)))

    # set mets to the mets contained within the background
    mets <- filteredMets
  }
  else {
    stop("background_type was not specified correctly. Please specify one of the following options: 'database', 'file', 'list', or 'biospecimen'.")
  }

  if(background_type == "database"){
      res <- chemicalClassSurveyRampIdsFullPopConn(mets, conn)
  } else {
    res <- chemicalClassSurveyRampIdsConn(mets, bkgrnd, conn)
  }
  RMariaDB::dbDisconnect(conn)

  print("Finished Chemical Class Survey")

  # It's possible that the chemical class survey result is empty, no hits for metabolites
  # or all metabolites are not in the background
  if(is.null(res) || length(res) == 0) {
    warning("Empty Chemical Class Result")
    return(list())
  }

  if(includeRaMPids){
    return(res)
  }else{
    if(!is.null(res$met_classes)) {
      res$met_classes<-res$met_classes %>% cleanup
      res$met_classes<-res$met_classes %>% cleanup
    }
    if(!is.null(res$pop_classes)) {
      res$pop_classes<-res$pop_classes %>% cleanup
      res$pop_classes<-res$pop_classes %>% cleanup
    }
    return(res)
  }
}


#' Returns chemical class information comparing a metabolite subset to a metabolite population, including Fisher Exact Test
#' enrichment p-values and FDR values.
#'
#' @param mets a vector of source prepended metabolite ids
#' @param background an optional list of source prepended metaboite ids to be used as the background reference of
#' metabolites for enrichment. The background can be either a list of ids, a file name containing the id list,
#' one id per column (no file header row) or a specificed biospecimen type (available biospecimen types: "Blood",
#' "Adipose tissue", "Heart", "Urine", "Brain", "Liver", "Kidney","Saliva", or "Feces").
#' @param background_type one of 'database' (all analytes in the RaMP Database), 'list' (a list of input ids),
#' or 'file' in which case the background parameter will be a file path, or 'biospecimen' where the specified background parameter is
#' a RaMP HMDB metabolite ontology term (see background parameter, above. for the most common biospecimen background values).
#' @return a list of dataframes, each holding chemical classs enrichment statistics for specific chemical classification systems,
#' such as HMDB Classyfire class categories and LIPIDMAPS class categories.  The results list chemical classes, metabolite hits counts,
#' Fisher Exact p-values and Benjamini-Hochberg corrected p-values (FDR estimates)
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

#' enrichedClassStats <- chemicalClassEnrichment(mets = metList)
#'}
#' @export
chemicalClassEnrichment <- function(mets, background = "database", background_type = "database") {
  print("Starting Chemical Class Enrichment")

  classData <- chemicalClassSurvey(mets = mets,
                                   background = background,
                                   background_type = background_type,
                                   includeRaMPids = TRUE)

  # Quit if empty survey result
  if(is.null(classData) || length(classData) < 1) {
    warning("Empty survey result. Please see other warnings.")
    return(list())
  }

  enrichmentStat <- list()

  totalCountInfo <- getTotalFoundInCategories(classData)

  breakPoint <- TRUE

  for (category in names(classData$count_summary)) {

    categoryData <- classData$count_summary[[category]]

    # only run if we have data for that class category
    # ANNNND we have that category represented in $mets
    if(nrow(categoryData) > 0 && category %in% names(totalCountInfo$mets)) {
      resultRow <- 1

      category <- category
      totMetCnt <- totalCountInfo$mets[[category]]
      totPopCnt <- as.integer(totalCountInfo$pop[[category]])
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

          row = list(as.character(category), as.character(className), contingencyMat[1,1], contingencyMat[1,1] + contingencyMat[2,1],
                     totMetCnt, totPopCnt, p)

          resultMat[resultRow, ] <- row

          resultRow <- resultRow + 1
        }
      }
      resultMat
      resultMat <- bhCorrect(resultMat)
      enrichmentStat[[as.character(category)]] <- resultMat
    } else {
      # just append an empty data frame place holder for that category type.
      resultMat <- data.frame(matrix(ncol=7, nrow=0))
      colnames(resultMat) <- c("category", "class_name", "met_hits", "pop_hits",
                               "met_size", "pop_size", "p-value")
      enrichmentStat[[as.character(category)]] <- resultMat
    }
  }
  enrichmentStat[['result_type']] <- 'chemical_class_enrichment'
  print("Finished Chemical Class Enrichment")
  return(enrichmentStat)
}

