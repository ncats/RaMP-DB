#' Retrieves analytes that involved in same reaction as input metabolite
#'
#' @param analytes a vector of analytes that need to be searched
#' @param NamesOrIds whether input is "names" or "ids" (default is "ids")
#' @return a dataframe containing query results. If the input is a metabolite, the function will output
#' gene transcript common names and source IDs that are known to catalyze
#' reactions in the same pathway as that metabolite. Conversely, if the input
#' is a gene, the function will return the common name and source id of metabolites
#' known to be catalyzed directly or indirectly by that gene. Input ids and common names will be returned.
#' If no input ids or names are found in the database, the return value will be an empty data frame, 0 rows.
#'
#' @examples
#' \dontrun{
#' pkg.globals <- setConnectionToRaMP(dbname="ramp2",username="root",conpass="",host = "localhost")
#' rampFastCata(analytes="creatine",NameOrIds="names")
#' }
#' @export
rampFastCata <- function(db = RaMP(), analytes="none", NamesOrIds="ids") {

  rampId <- pathwayRampId <- c()
  if(length(analytes)==1){
    if(analytes=="none"){
      stop("Please provide input analytes")}}
  if (!(NamesOrIds %in% c('names','ids')))
    stop('Please specify search by "names" or "ids"')

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

  isSQLite = .is_sqlite(db)

  if(NamesOrIds == 'ids') {

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

    df1 <- RaMP::runQuery(metQuery, db)

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

    df2 <- RaMP::runQuery(geneQuery, db)

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

    df1 <- RaMP::runQuery(metQuery, db)

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

    df2 <- RaMP::runQuery(geneQuery, db)

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

  rheaResult <- RaMP::getRheaAnalyteReactionAssociations(db=db, analytes=analytes, includeRheaRxnDetails = F, humanProteins = T)

  resultList <- list()
  resultList[["HMDB_Analyte_Associations"]] <- result
  resultList[["Rhea_Analyte_Associations"]] <- rheaResult

  return(resultList)
}






#' Retrieves analytes that involved in same reaction as input metabolite
#'
#' @param analytes a vector of analytes that need to be searched
#' @param NamesOrIds whether input is "names" or "ids" (default is "ids")
#' @return a dataframe containing query results. If the input is a metabolite, the function will output
#' gene transcript common names and source IDs that are known to catalyze
#' reactions in the same pathway as that metabolite. Conversely, if the input
#' is a gene, the function will return the common name and source id of metabolites
#' known to be catalyzed directly or indirectly by that gene.
#'
#' @examples
#' \dontrun{
#' pkg.globals <- setConnectionToRaMP(dbname="ramp2",username="root",conpass="",host = "localhost")
#' rampFastCata(analytes="creatine",NamesOrIds="names")
#' }
rampFastCataOriginal <- function(db = RaMP(), analytes="none", NamesOrIds="ids") {
    if(length(analytes)==1){
        if(analytes=="none"){
    stop("Please provide input analytes")}}
  if (!(NamesOrIds %in% c('names','ids')))
    stop('Please specify search by "names" or "ids"')

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

  #  print(list_metabolite)

  # Retrieve RaMP analyte ids
  if (NamesOrIds == 'names'){
    #    query1 <- paste0("select Synonym as analyte1,rampId,geneOrCompound as type1 from analytesynonym where Synonym in (",list_metabolite,");")
    query1 <- paste0("select rampId,geneOrCompound as type1,Synonym as InputAnalyte from analytesynonym where Synonym in (",list_metabolite,");")
  } else if (NamesOrIds == 'ids'){
    #    query1 <- paste0('select rampId,geneOrCompound as type1,commonName as InputMetabolite from analytesynonym where rampId in (select rampId from source where sourceId in (',list_metabolite,'));')
    query1 <- paste0('select rampId,geneOrCompound as type1,commonName as InputMetabolite from source where sourceId in (',list_metabolite,');')
  }

  # Retrieves Name, RaMPID and type (gene or compound) for input
  df1 <- RaMP::runQuery(query1, db)

  #print(df1$rampId)
  df_c <- df_g <- NULL
  mdf_c <- mdf_g <- NULL
  # Process metabolite ids
  mdf_cfin2 <- mdf_gfin2 <- c()
  if(length(grep("RAMP_C",df1$rampId))!=0){
    df_c <- df1[grep("RAMP_C",df1$rampId),]
    print("Get Compound ...")
    c_id <- unique(df_c$rampId)
    if(length(c_id) == 0){
      message("Input metabolites do not have catalyzation information")
      mdf_cfin2 <- NULL #return(NULL)
    } else {
      print(length(c_id))
      c_id <- sapply(c_id,shQuote)
      c_id <- paste(c_id,collapse = ",")

      # Retrieve rampid of genes that are in same reaction
      query_c <- paste0("select rampCompoundId as rampId,rampGeneId as rampId2 from catalyzed where rampCompoundId in (",c_id,");")
      print("Geting gene Id from Compound Id ...")

      df_c2 <- RaMP::runQuery(query_c, db)

      if(nrow(df_c2) == 0){
        message("No genes found in same reaction as input metabolite")
        mdf_cfin2 <- NULL
      } else {
        print("Getting names from gene Id ...")
        analyte2_list <- unique(df_c2$rampId2)
        analyte2_list <- sapply(analyte2_list,shQuote)
        analyte2_list <- paste(analyte2_list,collapse = ",")
        # Get names for metabolite ids
        query2 <- paste0("select * from source
             		where rampId in (",analyte2_list,");")

        df_c3 <- RaMP::runQuery(query2, db)

        if(nrow(df_c3) == 0){
          message("Cannot retrieve names for those metabolites")
          mdf_cfin2 <- NULL
        } else {
          # Now merge it all:
          mdc_c <- merge(df_c,df_c2)
          colnames(df_c3)[which(colnames(df_c3)=="rampId")]="rampId2"
          mdf_cfin <- merge(mdc_c,df_c3)
          #print(colnames(mdf_cfin))
          if (NamesOrIds == 'names'){
            mdf_cfin <- mdf_cfin[,c("InputAnalyte","sourceId","IDtype","commonName")]
          } else if (NamesOrIds == 'ids'){
            mdf_cfin <- mdf_cfin[,c("InputMetabolite","sourceId","IDtype","commonName")]
          }
          colnames(mdf_cfin) <- c("Input_Metabolite","Gene_sourceId","Gene_IDtype",
                                  "Gene_CommonName")

          # Collapse source ids:
          mdf_cfin$temp <- paste(mdf_cfin[,"Input_Metabolite"],mdf_cfin[,"Gene_CommonName"])
          tempout <- data.frame(Input_Metabolite=NA,Gene_CommonName=NA,Gene_sourceIds=NA)
          mdf_cfin2=c()
          for (i in unique(mdf_cfin$temp)) {
            temp <- mdf_cfin[which(mdf_cfin$temp==i),]
            tempout$Input_Metabolite=temp[1,"Input_Metabolite"]
            tempout$Gene_sourceIds <-
              paste(temp$Gene_sourceId,
                    collapse="; ")
            tempout$Gene_CommonName=temp[1,"Gene_CommonName"]
            mdf_cfin2 <- rbind(mdf_cfin2,tempout)
          }
          colnames(mdf_cfin2) <- c("Input_Analyte","Input_CatalyzedBy_CommonName",
                                   "Input_CatalyzedBy_SourceIds")
        } # end else couldn't retrieve names for metabolites
      } # end else couldn't find metabolite ids
    } # no catalyzation information
  } # end if compound ids found

  # Do analagous for genes
  if(length(grep("RAMP_G",df1$rampId))!=0){
    print("Also find gene inside")
    df_g <- df1[grep("RAMP_G",df1$rampId),]
    print("Get gene ...")
    g_id <- df_g$rampId
    g_id <- sapply(unique(g_id),shQuote)
    g_id <- paste(g_id,collapse = ",")

    if(length(g_id) == 0){
      message("No IDs found for input genes")
      mdf_gfin2 <- NULL #return(NULL)
    } else {
      # Get rampID for genes and catalyzed metabolites
      query_g <- paste0("select * from catalyzed where rampGeneId in (",g_id,");")

      df_g2 <- RaMP::runQuery(con,query_g, db)

      if(nrow(df_g2) == 0){
        message("Could not find metabolites in same reaction as input genes")
        mdf_gfin2=c()
      } else {
        analyte2_list <- df_g2$rampCompoundId
        analyte2_list <- sapply(analyte2_list,shQuote)
        analyte2_list <- paste(analyte2_list,collapse = ",")

        # Get names for metabolite IDs
        query2 <- paste0("select * from source where rampId in (",analyte2_list,");")

        df_g3 <-RaMP::runQuery(query2, db)

        if(nrow(df_g3) == 0){
          message("Cannot retrieve names for those genes")
          mdf_gfin2 <- NULL
        } else {
          # Now merge it all:
          mdc_g <- merge(df_g,df_g2)
          colnames(df_g3)[which(colnames(df_g3)=="rampId")]="rampId2"
          mdf_gfin <- merge(mdc_g,df_g3)
          colnames(mdf_gfin)[colnames(mdf_gfin) == 'InputMetabolite'] = 'InputAnalyte'
          print(colnames(mdf_gfin))
          mdf_gfin <- mdf_gfin[,c("InputAnalyte","sourceId","IDtype","commonName")]
          colnames(mdf_gfin) <- c("Input_Gene","Gene_sourceId","Gene_IDtype",
                                  "Gene_CommonName")
          mdf_gfin <- mdf_gfin[!duplicated(mdf_gfin),]

          # Collapse source ids:
          mdf_gfin$temp <- paste(mdf_gfin[,"Input_Gene"],mdf_gfin[,"Gene_CommonName"])
          tempout <- data.frame(Input_Analyte=NA,Input_CatalyzedBy_CommonName=NA,
                                Input_CatalyzedBy_SourceIds=NA)
          mdf_gfin2=c()
          for (i in unique(mdf_gfin$temp)) {
            temp <- mdf_gfin[which(mdf_gfin$temp==i),]
            tempout$Input_Analyte=temp[1,"Input_Gene"]
            tempout$Input_CatalyzedBy_CommonName=temp[1,"Gene_CommonName"]
            tempout$Input_CatalyzedBy_SourceIds <-
              paste(temp$Gene_sourceId,collapse="; ")
            mdf_gfin2 <- rbind(mdf_gfin2,tempout)
          }
        } # else couldn't retrieve names for those genes
      } # couldn't find metabolites catalyzed by input genes
    } # end else couldn't find ids for input gene names
  } # end gene
  mdf <- rbind(mdf_cfin2,mdf_gfin2)
  print("Done ...")
  print("timing ...")
  print(proc.time() - now)
  #  if(!is.null(mdf)) {
  #  	colnames(mdf) <- c("Input_Analyte","Input_CatalyzedBy_CommonName",
  #		"Input_CatalyzedBy_SourceIds")
                                        #   }
  return(mdf)
}
