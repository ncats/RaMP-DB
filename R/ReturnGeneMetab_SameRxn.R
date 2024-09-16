#' Retrieves analytes that involved in same reaction as input metabolite
#'
#' @param analytes a vector of analytes that need to be searched
#' @param namesOrIds whether input is "names" or "ids" (default is "ids")
#' @param db a RaMP database object
#' @return a list of two dataframes containing query results from HMDB and Rhea. If the input is a metabolite, the function will output
#' gene transcript common names and source IDs that are known to catalyze
#' reactions in the same pathway as that metabolite. Conversely, if the input
#' is a gene, the function will return the common name and source id of metabolites
#' known to be catalyzed directly or indirectly by that gene. Input ids and common names will be returned.
#' If no input ids or names are found in the database, the return value will be an empty data frame, 0 rows.
#'
#' @examples
#' \dontrun{
#' inputs.of.interest <- c("kegg:C00186" , "hmdb:HMDB0000148", "kegg:C00780",
#'  "hmdb:HMDB0000064", "ensembl:ENSG00000115850", "uniprot:Q99259")
#'
#' new.transcripts <- rampFastCata( analytes = inputs.of.interest, db = rampDB )
#' }
#' @export
rampFastCata <- function( analytes="none", namesOrIds="ids", db = RaMP() ) {

  rampId <- pathwayRampId <- c()
  if(length(analytes)==1){
    if(analytes=="none"){
      stop("Please provide input analytes")}}

  if (!(namesOrIds %in% c('names','ids'))){
    stop('Please specify search by "names" or "ids"')
  }

  now <- proc.time()
  if(is.character(analytes)){
    if(grepl("\n",analytes)[1]){
      list_metabolite <- strsplit(analytes,"\n")
      list_metabolite <- unlist(list_metabolite)
    } else if(grepl(",",analytes)[1]){
      list_metabolite <- strsplit(analytes,"\n")
      list_metabolite <- unlist(list_metabolite)
    } else {
      list_metabolite <- analytes
    }
  } else if(is.data.frame(analytes)){
    list_metabolite <- unlist(analytes)
  } else {stop("The input 'analytes' is not a recognized format. Please check input.")}

  list_metabolite <- unique(list_metabolite)
  list_metabolite <- sapply(list_metabolite,shQuote)
  list_metabolite <- paste(list_metabolite,collapse = ",")

  data_access <- DataAccessObject$new(db = db)

  if(namesOrIds == 'ids') {

    print("Analyte ID-based reaction partner query.")

    df1 <- data_access$getRxnPartnersFromMetIDs(metaboliteIDs = list_metabolite)
    print("Building metabolite to gene relations.")

    print(paste0("Number of met2gene relations: ",(nrow(df1))))

    df2 <- data_access$getRxnPartnersFromGeneIDs(geneIDs = list_metabolite)
    print("Building gene to metabolite relations.")

  } else {

    # working on 'names' query
    # note that we now bring in the synonyms table

    print("Analyte name-based reaction partner query.")
    print("Building metabolite to gene relations.")
    df1 <- data_access$getRxnPartnersFromMetNames(metaboliteNames = list_metabolite)

    print(paste0("Number of met2gene relations: ",(nrow(df1))))
    df2 <- data_access$getRxnPartnersFromGeneNames(geneNames = list_metabolite)

    print("Building gene to metabolite relations.")
    print(paste0("Number of gene2met relations: ",(nrow(df2))))
  }

  if(!is.null(df1) && nrow(df1) > 0) {
    df1$query_relation <- 'met2gene'
    result <- df1
    if(!is.null(df2) && nrow(df2) > 0) {
      df2$query_relation <- 'gene2met'
      result <- rbind(result, df2)
    }
  } else {
    if(!is.null(df2) && nrow(df2) > 0) {
      df2$query_relation <- 'gene2met'
      result <- df2
    } else {
      # default handling of empty result
      # empty df1 requires use of tibble/tidyr add_column
      df1 <- tibble::add_column(df1, 'query_relation'=NA)
      result <- df1
    }
  }

  # remove rampId column
  result <- subset(result, select=-c(rampId))
  # move relation first
  result <- result[,c(ncol(result), 1:(ncol(result)-1))]

  print(paste0("Total Relation Count: ", (nrow(result))))

  rheaResult <- RaMP::getRheaAnalyteReactionAssociations(db=db, analytes=analytes, includeRheaRxnDetails = F, humanProtein = T)

  resultList <- list()
  resultList[["HMDB_Analyte_Associations"]] <- result
  resultList[["Rhea_Analyte_Associations"]] <- rheaResult



  if(nrow(resultList$HMDB_Analyte_Associations) >0)
  {
    resultList$HMDB_Analyte_Associations$rxn_partner_ids <- NA
    resultList$HMDB_Analyte_Associations$Source <- "HMDB"
  }

  if(nrow(resultList$Rhea_Analyte_Associations) >0)
  {
    colnames(resultList$Rhea_Analyte_Associations)[3] <- "input_common_name"
    resultList$Rhea_Analyte_Associations$Source <- "Rhea"
  }

  if (nrow(resultList$HMDB_Analyte_Associations) >0 && nrow(resultList$Rhea_Analyte_Associations) >0)
  {
    resultDF <- rbind(resultList$HMDB_Analyte_Associations, resultList$Rhea_Analyte_Associations)
    resultDF[which(do.call(paste0, resultDF[,3:4]) %in% do.call(paste0, resultDF[duplicated(resultDF[3:4]),3:4])),]$Source <- "Both"
    duplicates <- subset(resultDF, resultDF$Source=='Both')
    resultDF <- subset(resultDF, resultDF$Source!='Both')

    resultDF <- rbind(resultDF, duplicates[!is.na(duplicates$rxn_partner_ids),])

  } else if (nrow(resultList$HMDB_Analyte_Associations) >0)
  {
    resultDF <- resultList$HMDB_Analyte_Associations
  } else if (nrow(resultList$Rhea_Analyte_Associations) >0)
  {
    resultDF <- resultList$Rhea_Analyte_Associations
  }


  return(resultDF)
}

