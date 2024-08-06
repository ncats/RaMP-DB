#' Function that query database to find ontology information based on
#' the given list of analytes
#' @param analytes a vector of analytes or a analytes delimited by new line character
#' @param namesOrIds specify the type of given data
#' @param db a RaMP database object
#' @return dataframe that contains searched ontology from given analytes
#'
#' @examples
#' \dontrun{
#' analytes.of.interest <- c("chebi:15422", "hmdb:HMDB0000064", "hmdb:HMDB0000148", "wikidata:Q426660")
#'
#' new.ontologies <- RaMP::getOntoFromMeta(db = rampDB, analytes = analytes.of.interest)
#' }
#' @export
getOntoFromMeta <- function(analytes, namesOrIds = "ids", db = RaMP()) {
  if (!(namesOrIds %in% c("ids", "names"))) {
    stop("Specifiy the type of given data to 'ids' or 'names'")
  }

  now <- proc.time()
  if (is.character(analytes)) {
    if (grepl("\n", analytes)[1]) {
      list_metabolite <- strsplit(analytes, "\n")
      list_metabolite <- unlist(list_metabolite)
    } else if (grepl(",", analytes)[1]) {
      list_metabolite <- strsplit(analytes, "\n")
      list_metabolite <- unlist(list_metabolite)
    } else {
      list_metabolite <- analytes
    }
  } else if (is.data.frame(analytes)) {
    list_metabolite <- unlist(analytes)
  }
  list_metabolite <- unique(list_metabolite)
  list_metabolite <- sapply(list_metabolite, shQuote)
  list_metabolite <- paste(list_metabolite, collapse = ",")

  if (namesOrIds == "ids") {
    sql <- paste0("select * from source where sourceId in (", list_metabolite, ");")
  } else if (namesOrIds == "names") {
    sql <- paste0("select * from source where rampId in (select * from (select rampId from analytesynonym where Synonym in (", list_metabolite, ")) as subquery);")
    cat(file = stderr(), "query sql in Package call with -- ", sql, "\n")
  }

  df <- RaMP::runQuery(sql, db)

  if (nrow(df) == 0) {
    message("This source id
            does not exist in the source table")
    return(NULL)
  }

  rampid <- unique(df$rampId)
  rampid <- sapply(rampid, shQuote)
  rampid <- paste(rampid, collapse = ",")

  sql <- paste0(
    "select * from analytehasontology where rampCompoundId in (",
    rampid, ");"
  )

  df2 <- RaMP::runQuery(sql, db)

  if (nrow(df2) == 0) {
    message("No searching result because these metabolites are not linked to ontology")
    return(NULL)
  }

  rampontoid <- unique(df2$rampOntologyId)
  rampontoid <- sapply(rampontoid, shQuote)
  rampontoid <- paste(rampontoid, collapse = ",")
  sql <- paste0(
    "select * from ontology where rampOntologyId in (",
    rampontoid, ");"
  )

  df3 <- RaMP::runQuery(sql, db)

  mdf <- unique(merge(df3, df2, all.x = T))
  mdf <- unique(merge(mdf, df,
    all.x = T, by.x = "rampCompoundId",
    by.y = "rampId"
  ))
  colnames(mdf)[colnames(mdf) == "commonName.x"] <- "Ontology"
  colnames(mdf)[colnames(mdf) == "commonName.y"] <- "Metabolites"

  mdf <- mdf[c("Metabolites", "sourceId", "IDtype", "Ontology", "HMDBOntologyType")]

  # need to make unique list (JB, 12/1/21)
  mdf <- unique(mdf)

  # colnames(mdf) <- c('Metabolites','Ontology','Ontology_Type')
  return(mdf)
}


#' function that query database to find analytes in given ontologies
#' @param ontology a vector of ontology or ontologies delimited by new line character
#' @param db a RaMP database object
#' @return dataframe that  contains searched analytes from given ontology
#' @examples
#' \dontrun{
#' ontologies.of.interest <- c("Colon", "Liver", "Lung")
#'
#' new.metabolites <- RaMP::getMetaFromOnto(db = rampDB, ontology = ontologies.of.interest)
#' }
#' @importFrom rlang .data
#' @export
getMetaFromOnto <- function(ontology, db = RaMP()) {

  print("Retreiving Metabolites for input ontology terms.")
  now <- proc.time()
  if(is.character(ontology)){
    if(grepl("\n",ontology)[1]){
      list_ontology <- strsplit(ontology,"\n")
      list_ontology <- unlist(list_ontology)
    } else if(grepl(",",ontology)[1]){
      list_ontology <- strsplit(ontology,",")
      list_ontology <- unlist(list_ontology)
    } else {
      list_ontology <- ontology
    }
  } else if(is.data.frame(ontology)){
    list_ontology <- unlist(ontology)
  }

  list_ontology <- unique(list_ontology)

  allontos <- getOntologies(db = db)
  matched_ontos <- unlist(lapply(list_ontology,
                                 function(x) grep(paste0("^",x,"$"), allontos$commonName)))

  # Only proceed if df has anything returned
  if (length(matched_ontos) > 0) {

    print(paste0("Found ",length(matched_ontos)," ontology term matches."))

    ontologyList <- paste0("'", paste(list_ontology, collapse="','"),"'")

    sql = paste0("select rampId,
          group_concat(distinct s.sourceId order by s.sourceId separator '; ') as source_ids,
          group_concat(distinct s.commonName order by s.commonName separator '; ') as common_names, o.commonName, o.HMDBOntologyType
          from source s, analytehasontology ao, ontology o where ao.rampOntologyId in (
          select distinct rampOntologyId from ontology where commonName in (",ontologyList,"))
          and o.rampOntologyId = ao.rampOntologyId and s.rampId = ao.rampCompoundId
          group by o.commonName, s.rampId, o.HMDBOntologyType")

    if(RaMP:::.is_sqlite(db)) {
      sql = paste0("select rampId,
          group_concat(distinct s.sourceId COLLATE NOCASE) as source_ids,
          group_concat(distinct s.commonName COLLATE NOCASE) as common_names, o.commonName, o.HMDBOntologyType
          from source s, analytehasontology ao, ontology o where ao.rampOntologyId in (
          select distinct rampOntologyId from ontology where commonName in (",ontologyList,"))
          and o.rampOntologyId = ao.rampOntologyId and s.rampId = ao.rampCompoundId
          group by o.commonName, s.rampId, o.HMDBOntologyType")
    }

    mdf_final <- RaMP::runQuery(sql, db)

    mdf_final <- unique(mdf_final)
    mdf_final <- mdf_final[,c(4,5,3,2)]
    colnames(mdf_final) <- c('ontologyTerm', 'ontologyCategory', 'metNames', 'metIds')

    print(paste0("Found ",nrow(mdf_final)," metabolites associated with the input ontology terms."))
    print("Finished getting metabolies from ontology terms.")

    return(mdf_final)
  } else {
    warning("The input ontology terms were not found in RaMP.\nRun the getOntologies() function to see available ontology terms.")
    return(NA)
  }
}
