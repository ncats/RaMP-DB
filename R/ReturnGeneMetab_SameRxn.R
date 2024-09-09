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

  isSQLite = .is_sqlite(x = db)

  if(namesOrIds == 'ids') {

    print("Analyte ID-based reaction partner query.")

    # Ugghhh... SQLite specific query changes...

    metQuery <- paste0("select c.sourceId as input_analyte, group_concat(distinct c.commonName order by c.commonName asc separator '; ') as input_common_names,
  group_concat(distinct g.commonName order by g.commonName asc separator '; ') as rxn_partner_common_name,
  group_concat(distinct g.sourceId order by g.sourceId asc separator '; ') as rxn_partner_ids,
  g.rampId from catalyzed r
  join source g on r.rampGeneId = g.rampId
  join source c on r.rampCompoundId = c.rampId
  where c.sourceId in (",list_metabolite,") group by g.rampId, c.sourceId")

    if(isSQLite) {
    metQuery <- paste0("select c.sourceId as input_analyte, group_concat(distinct c.commonName COLLATE NOCASE) as input_common_names,
  group_concat(distinct g.commonName COLLATE NOCASE) as rxn_partner_common_name,
  group_concat(distinct g.sourceId COLLATE NOCASE) as rxn_partner_ids,
  g.rampId from catalyzed r
  join source g on r.rampGeneId = g.rampId
  join source c on r.rampCompoundId = c.rampId
  where c.sourceId in (",list_metabolite,") group by g.rampId, c.sourceId")
    }

    print("Building metabolite to gene relations.")

    df1 <- runQuery(sql = metQuery, db = db)

    print(paste0("Number of met2gene relations: ",(nrow(df1))))


    geneQuery <- paste0("select g.sourceId as input_analyte, group_concat(distinct g.commonName order by g.commonName asc separator '; ') as input_common_names,
  group_concat(distinct c.commonName order by c.commonName asc separator '; ') as rxn_partner_common_name,
  group_concat(distinct c.sourceId order by c.sourceId asc separator '; ') as rxn_partner_ids,
  c.rampId from catalyzed r
  join source g on r.rampGeneId = g.rampId
  join source c on r.rampCompoundId = c.rampId
  where g.sourceId in (", list_metabolite,") group by c.rampId, g.sourceId")

    if(isSQLite) {
    geneQuery <- paste0("select g.sourceId as input_analyte, group_concat(distinct g.commonName COLLATE NOCASE) as input_common_names,
  group_concat(distinct c.commonName COLLATE NOCASE) as rxn_partner_common_name,
  group_concat(distinct c.sourceId COLLATE NOCASE) as rxn_partner_ids,
  c.rampId from catalyzed r
  join source g on r.rampGeneId = g.rampId
  join source c on r.rampCompoundId = c.rampId
  where g.sourceId in (", list_metabolite,") group by c.rampId, g.sourceId")
    }
    print("Building gene to metabolite relations.")

    df2 <- runQuery(sql = geneQuery, db = db)

  } else {

    # working on 'names' query
    # note that we now bring in the synonyms table

    print("Analyte name-based reaction partner query.")

    metQuery <- paste0("select s.Synonym as input_analyte, group_concat(distinct c.commonName order by c.commonName asc separator '; ') as input_common_names,
  group_concat(distinct g.commonName order by g.commonName asc separator '; ') as rxn_partner_common_name,
  group_concat(distinct g.sourceId order by g.sourceId asc separator '; ') as rxn_partner_ids,
  g.rampId from catalyzed r
  join source g on r.rampGeneId = g.rampId
  join source c on r.rampCompoundId = c.rampId
  join analytesynonym s on s.rampId = r.rampCompoundId
  where s.Synonym in (",list_metabolite,") group by g.rampId, s.Synonym")

    if(isSQLite) {
      metQuery <- paste0("select s.Synonym as input_analyte, group_concat(distinct c.commonName COLLATE NOCASE) as input_common_names,
  group_concat(distinct g.commonName COLLATE NOCASE) as rxn_partner_common_name,
  group_concat(distinct g.sourceId COLLATE NOCASE) as rxn_partner_ids,
  g.rampId from catalyzed r
  join source g on r.rampGeneId = g.rampId
  join source c on r.rampCompoundId = c.rampId
  join analytesynonym s on s.rampId = r.rampCompoundId
  where s.Synonym in (",list_metabolite,") group by g.rampId, s.Synonym")
    }

    print("Building metabolite to gene relations.")

    df1 <- runQuery(sql = metQuery, db = db)

    print(paste0("Number of met2gene relations: ",(nrow(df1))))

    geneQuery <- paste0("select s.Synonym as input_analyte, group_concat(distinct g.commonName order by g.commonName asc separator '; ') as input_common_names,
  group_concat(distinct c.commonName order by c.commonName asc separator '; ') as rxn_partner_common_name,
  group_concat(distinct c.sourceId order by c.sourceId asc separator '; ') as rxn_partner_ids,
  c.rampId from catalyzed r
  join source g on r.rampGeneId = g.rampId
  join source c on r.rampCompoundId = c.rampId
  join analytesynonym s on s.rampId = r.rampGeneId
  where s.Synonym in (", list_metabolite,") group by c.rampId, s.Synonym")

    if(isSQLite) {
      geneQuery <- paste0("select s.Synonym as input_analyte, group_concat(distinct g.commonName COLLATE NOCASE) as input_common_names,
  group_concat(distinct c.commonName COLLATE NOCASE) as rxn_partner_common_name,
  group_concat(distinct c.sourceId COLLATE NOCASE) as rxn_partner_ids,
  c.rampId from catalyzed r
  join source g on r.rampGeneId = g.rampId
  join source c on r.rampCompoundId = c.rampId
  join analytesynonym s on s.rampId = r.rampGeneId
  where s.Synonym in (", list_metabolite,") group by c.rampId, s.Synonym")
    }
    print("Building gene to metabolite relations.")

    df2 <- runQuery(sql = geneQuery, db = db)

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

  return(resultList)
}

