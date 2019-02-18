#' processData function generates pathway RampId frequency (gene or metabolite) based on pathway source (hmdb,kegg,reactome,wiki)
#' @param conpass password for database access (string)
#' @param dbname name of the mysql database (default is "ramp")
#' @param username username for database access (default is "root")
#' @param host host name for database access (default is "localhost")
#' @return R object (FT_data.Rdata) with dataframes (hmdb_metab,hmdb_gene,kegg_gene,kegg_metab,reactome_gene,reactome_metab,wiki_gene,wiki_metab)

processData <- function(conpass=NULL,
                        dbname="ramp",username="root",
                        host = "localhost") {

  if(is.null(conpass)) {
    stop("Please define the password for the mysql connection")
  }

  # get all rows form analytehaspathway
  query <- "select * from analytehaspathway"
  con <- DBI::dbConnect(RMySQL::MySQL(), user = username,
                        password = conpass,
                        dbname = dbname,
                        host = host)

  allRampIds <- DBI::dbGetQuery(con,query)


  if(is.null(allRampIds)) {

    stop("Data doesn't exist")

  } else {


    # length(unique(allRampIds$pathwayRampId))
    # total unique pathwayRampIDs - 51,526

    # data frames for metabolites with pathawayRampID, Freq based  on Source(kegg, reactome, wiki, hmdb)

    allRampId_C <- allRampIds[grep("RAMP_C", allRampIds$rampId), ]
    unique_allRampId_C <- unique(allRampId_C[,c("rampId", "pathwayRampId")])
    unique_pathwayRampId_source <- unique(allRampId_C[,c("pathwayRampId", "pathwaySource")])

    # length(unique_pathwayRampId_source$pathwayRampId) - 36,039

    freq_unique_allRampId_C <- as.data.frame(table(unique_allRampId_C[,"pathwayRampId"]))

    # length of total RAMP_C pathwayRampIDs - 36,039

    names(freq_unique_allRampId_C)[1] = 'pathwayRampId'
    merge_Pathwayfreq_source <- merge(freq_unique_allRampId_C, unique_pathwayRampId_source, by="pathwayRampId")

    # subset data based on source -  kegg, reactome, wiki, hmdb

    kegg_metab <- subset(merge_Pathwayfreq_source, merge_Pathwayfreq_source$pathwaySource == "kegg")
    #length(kegg_metab$pathwayRampId) - 264
    reactome_metab <- subset(merge_Pathwayfreq_source, merge_Pathwayfreq_source$pathwaySource == "reactome")
    #length(reactome_metab$pathwayRampId) - 1866
    wiki_metab <- subset(merge_Pathwayfreq_source, merge_Pathwayfreq_source$pathwaySource == "wiki")
    #length(wiki_metab$pathwayRampId) - 213
    hmdb_metab <- subset(merge_Pathwayfreq_source, merge_Pathwayfreq_source$pathwaySource == "hmdb")
    #length(hmdb_metab$pathwayRampId) - 33,696



    # data frames for Genes with pathawayRampID, Freq based  on Source(kegg, reactome, wiki, hmdb)

    allRampId_G <- allRampIds[grep("RAMP_G", allRampIds$rampId), ]
    unique_allRampId_G <- unique(allRampId_G[,c("rampId", "pathwayRampId")])
    unique_pathwayG_source <- unique(allRampId_G[,c("pathwayRampId", "pathwaySource")])

    # length(unique(unique_pathwayG_source$pathwayRampId)) - 51,461

    freq_unique_allRampId_G <- as.data.frame(table(unique_allRampId_G[,"pathwayRampId"]))

    # length of total RAMP_G pathwayRampIDs - 51,461

    names(freq_unique_allRampId_G)[1] = 'pathwayRampId'
    merge_PathwayG_source <- merge(freq_unique_allRampId_G, unique_pathwayG_source, by="pathwayRampId")

    # subset data based on source -  kegg, reactome, wiki, hmdb

    kegg_gene <- subset(merge_PathwayG_source, merge_PathwayG_source$pathwaySource == "kegg")
    #length(kegg_gene$pathwayRampId) - 323
    reactome_gene <- subset(merge_PathwayG_source, merge_PathwayG_source$pathwaySource == "reactome")
    #length(reactome_gene$pathwayRampId) - 2155
    wiki_gene <- subset(merge_PathwayG_source, merge_PathwayG_source$pathwaySource == "wiki")
    #length(wiki_gene$pathwayRampId) - 408
    hmdb_gene <- subset(merge_PathwayG_source, merge_PathwayG_source$pathwaySource == "hmdb")
    #length(hmdb_gene$pathwayRampId) - 48575

    save(kegg_metab,
         reactome_metab,
         wiki_metab,
         hmdb_metab,
         kegg_gene,
         reactome_gene,
         wiki_gene,
         hmdb_gene,
         file = paste0(system.file(package = "RaMP"),"/extdata/FT_data.Rdata"))
  }
}

#' Load pathway overlap matrices for find_clusters function
#' @return A list of pathway overlap matrices and freqTables for clustering and fisher test

loadSysData <- function() {

  load(system.file(package = "RaMP",... = "extdata/genes_overlap_matrix.RData"))
  load(system.file(package = "RaMP",... = "extdata/metabolites_overlap_matrix.RData"))
  load(system.file(package = "RaMP",... = "extdata/analytes_overlap_matrix.RData"))
  load(system.file(package = "RaMP",... = "extdata/FT_data.Rdata"))
  return(list(
         genes_result = genes_result,
         metabolites_result = metabolites_result,
         analyte_result = analyte_result,
         hmdb_gene = hmdb_gene,
         hmdb_metab = hmdb_metab,
         kegg_gene = kegg_gene,
         kegg_metab = kegg_metab,
         reactome_gene = reactome_gene,
         reactome_metab = reactome_metab,
         wiki_gene = wiki_gene,
         wiki_metab = wiki_metab)
    )

}

#' sysdataObject() function generates sysdata.rda object where dataframes or matrices cannot be accessed by user,
#' but made availble for internal functions (runFisherTest() & findCluster())
#' @return sysdata.rda object in R directory with 3 Overlap Matrices (genes_result, metabolites_result, analyte_result) & dataframes generated from processData function.

sysdataObject <- function() {
  #load frequency tables and overlap matrices (genes_result, metabolites_result, analyte_result) using loadSysData() function
  sysDat <- loadSysData()
  stopifnot(is.matrix(sysDat$genes_result),
            is.matrix(sysDat$metabolites_result),
            is.matrix(sysDat$analyte_result),
            is.data.frame(sysDat$hmdb_gene),
            is.data.frame(sysDat$hmdb_metab),
            is.data.frame(sysDat$kegg_gene),
            is.data.frame(sysDat$kegg_metab),
            is.data.frame(sysDat$reactome_gene),
            is.data.frame(sysDat$reactome_metab),
            is.data.frame(sysDat$wiki_gene),
            is.data.frame(sysDat$wiki_metab))

  # uncomment devtools::use_data function to create sysdata.Rda object
  # devtools::use_data(sysDat$genes_result,
  #                    sysDat$metabolites_result,
  #                    sysDat$analyte_result,
  #                    sysDat$hmdb_gene,
  #                    sysDat$hmdb_metab,
  #                    sysDat$kegg_gene,
  #                    sysDat$kegg_metab,
  #                    sysDat$reactome_gene,
  #                    sysDat$reactome_metab,
  #                    sysDat$wiki_gene,
  #                    sysDat$wiki_metab,
  #                    overwrite = TRUE,
  #                    internal = TRUE)

}

#Steps to test dev package
#step1: remove RaMP package in R LibPath
#step2: $cd /RaMP-DB
#step3: open R console at RaMP-DB dir
#step4: call library(devtools) and devtools::install()
#step5: Access RaMP library using R lib path or call library(RaMP)
#access data: RaMP:::

