# queries to retrieve and analyze chemical properties

#' Returns chemical properties given a metabolite list
#'
#' @param mets a list object of source prepended metabolite ids, representing a metabolite set of interest
#' @param propertyList an optional list of specific properties to extract.  Options include 'all' (default),  'smiles', 'inchi_key', 'inchi_key_prefix', 'inchi', 'mw', 'monoisotop_mass', 'formula', 'common_name'.
#' If a props list is not supplied, all property fields will be returned.
#' @return Returns chemical property information for the list of input metabolites and a query report reporting on the number of metabolite ids that were matched and the list of un-matched input ids.
#'
#' The returned object (return_obj below) contains two results. Use str(return_obj) to see the structure described here.
#'
#' \strong{return_obj$chem_props} chemical properties for all matched input ids.
#'
#' \strong{return_obj$query_report} this reports on the query list size, the number of input ids
#' that were found in the database, and a list of metabolite ids that were not found in the database.
#'@examples
#'\dontrun{
#' # metabolite list of interest
#' metabolites.of.interest = c("pubchem:64969",
#'                              "chebi:16958",
#'                              "chemspider:20549",
#'                              "kegg:C05598",
#'                              "chemspider:388809",
#'                              "pubchem:53861142",
#'                              "hmdb:HMDB0001138",
#'                              "hmdb:HMDB0029412")
#'
#' chemical.classes <- getChemicalProperties(mets = metabolites.of.interest, db = rampDB)
#'}
#' @export
getChemicalProperties <- function(mets, propertyList = 'all', db = RaMP() ){

  message("Starting Chemical Property Query")

  mets <- unique(mets)
  checkIdPrefixes(mets)
  result <- list()

  # first handle metabolites of interest
  metStr <- paste(mets, collapse = "','")
  metStr <- paste("'" ,metStr, "'", sep = "")

  if(length(grep("all",propertyList))==1) {
    sql <- paste0("select * from chem_props where chem_source_id in (",metStr,")")
  } else {
      propList <- buildPropertyList(db, propertyList);
      if(startsWith(propList, "Error")) {
        message(propList)
        return(NULL)
      }
      sql <- paste("select",propList,"from chem_props",
                   "where chem_source_id in (",metStr,")")
  }

  metsData <- RaMP:::runQuery(sql, db)
  foundMets <- unique(metsData$chem_source_id)

  result[['chem_props']] <- metsData

  queryNotes <- queryReport(mets, foundMets)

  result[['query_report']] <- queryNotes

  result[['chem_props']] <- result[['chem_props']] %>% cleanup

  if(!is.null(result)) message("Finished Chemical Property Query")

  # if(typeof(result) == 'list') {
  #   result <- data.frame(result)
  # }

  return(result)
}

# Internal function to validate property list
# @param propList an optional list of specific properties to extract.  Options include 'all' (default),  'iso_smiles', 'inchi_key', 'inchi_key_prefix', 'inchi', 'mw', 'monoisotop_mass', 'formula', 'common_name'.
buildPropertyList <- function(db = RaMP(), propList) {

  # validate that all properties are valid
  #  validProperties <- c('smiles', 'inchi_key', 'inchi_key_prefix', 'inchi', 'mw', 'monoisotop_mass', 'formula', 'common_name')

  if(.is_sqlite(db)) {
    sql = 'pragma table_info(chem_props)'
    ramptypes <- RaMP:::runQuery(sql, db)
    ramptypes <- unlist(ramptypes$name)
  } else {
    sql = 'describe chem_props'
    ramptypes <- RaMP:::runQuery(sql, db)
    ramptypes <- unlist(ramptypes$Field)
  }




  validProperties <- setdiff(ramptypes, c("ramp_id", "chem_data_source", "chem_source_id"))


  haveInvalidProps = FALSE


  invalidProps = c()
  for(prop in propList) {
    if (!(prop %in% validProperties)) {
      haveInvalidProps = TRUE
      invalidProps <- c(invalidProps, prop)
    }
  }
  if(haveInvalidProps) {
    propword <- 'Error: A requested chemical property name is invalid.'
    if(length(invalidProps) > 1) propword <- 'Error: Some requested chemical property names are invalid.'
    errorMsg <- paste0(utils::str(length(invalidProps)), propword, "\n")
    #errorMsg <- paste0(errorMsg, "Valid property names: 'smiles', 'inchi_key', 'inchi_key_prefix', 'inchi', 'mw', 'monoisotop_mass', 'formula', 'common_name'\n")
    errorMsg <- paste0(errorMsg, "Valid property names:", paste(validProperties,collapse=", "))
    errorMsg <- paste0(errorMsg, "Invalid input property name list:")
    for(badName in invalidProps) {
      errorMsg <- paste(errorMsg, badName)
    }
    return(errorMsg)
  }

  propStr <- paste(propList, collapse = ",")

  #  propStr <- gsub("smiles", "iso_smiles", propStr)
  #  propStr <- gsub("formula", "mol_formula", propStr)

  propStr <- paste("chem_source_id, ramp_id,",propStr)

  return(propStr)
}
