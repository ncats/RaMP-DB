#' Do fisher test for only one pathway from search result
#'
#' @param analytes a vector of analytes (genes or metabolites) that need to be searched
#' @param namesOrIds whether input is "names" or "ids" (default is "ids", must be the same for analytes and background)
#' @param total_genes number of genes analyzed in the experiment (e.g. background) (default is 20000, with assumption that analyte_type is "genes")
#' @param analyte_type "metabolites" or "genes" (default is "metabolites")
#' @param alternative alternative hypothesis test passed on to fisher.test().  Options are two.sided, greater, or less (default is "less")
#' @param minPathwaySize the minimum number of pathway members (genes and metabolites) to include the pathway in the output (default = 5)
#' @param maxPathwaySize the maximum number of pathway memnbers (genes and metaboltes) to include the pathway in the output (default = 150)
#' @param background_type type of background that is input by the user.  Opions are "database" if user wants all
#' analytes from the RaMP database will be used; "file", if user wants to input a file with a list of background
#' analytes; "list", if user wants to input a vector of analyte IDs; "biospecimen", if user wants to specify a
#' biospecimen type (e.g. blood, adipose tissue, etc.) and have those biospecimen-specific analytes used.  For genes,
#' only the "database" option is used.
#' @param background background to be used for Fisher's tests.  If parameter 'background_type="database"', this parameter
#' is ignored (default="database"); if parameter 'background_type= "file"', then 'background' should be a file name (with
#' directory); if 'background_type="list"', then 'background' should be a vector of RaMP IDs; if 'backgroud_type="biospecimen"'
#' then users should specify one of the following: "Blood", "Adipose tissue", "Heart", "Urine", "Brain", "Liver", "Kidney",
#' "Saliva", and "Feces"
#' @param pathway_definitions If "RaMP" (default), use pathway definitions within RaMP-DB. Else, supply path to gmx file containing custom pathway definitions. GMX files are a tab-separated format that contain one analyte set per column, with the name of the set in the first row, and constituent analytes in subsequent rows
#' @param include_smpdb Include pathways from smpdb/hmdb in analysis. Excluded by default since definitions are highly redundant
#' @param db a RaMP databse object
#' @return a dataframe with columns containing pathway ID, fisher's p value, user analytes in pathway, and total analytes in pathway

runFisherTest <- function(analytes,
                          total_genes = 20000,
                          namesOrIds = "ids",
                          analyte_type = "metabolites",
                          alternative = "less",
                          minPathwaySize = 5, maxPathwaySize = 150,
                          background_type = "database", background = "database",
                          pathway_definitions = "RaMP", include_smpdb = FALSE,
                          db = RaMP()) {
  if(analyte_type == "genes"){
    background_type = "database"
    print("Using database background for genes")
  }
  now <- proc.time()
  print("Fisher Testing ......")
  if (pathway_definitions != "RaMP") {
    pathwaydf <- getCustomPathwayFromAnalyte(analytes=analytes,
      pathway_definitions=pathway_definitions,
      analyte_type = analyte_type
      )
  } else {
    pathwaydf <- getPathwayFromAnalyte(db = db, analytes = analytes,
    includeRaMPids = TRUE,
      namesOrIds = namesOrIds,
      find_synonym = FALSE,
      include_smpdb = include_smpdb
    )
  }

  pathwayRampId <- rampId <- c()

  pathwayRampId <- rampId <- c()

  if (analyte_type == "metabolites") {
    pathwaydf <- pathwaydf[grep("RAMP_C_", pathwaydf$rampId), ]
  } else if (analyte_type == "genes") {
    pathwaydf <- pathwaydf[grep("RAMP_G_", pathwaydf$rampId), ]
  }


  # moved this check until we determine if we have analytes of a given type.
  if (nrow(pathwaydf) == 0) {
    return(NULL)
  }

#  if(class(background_type)=="list"){
   if(is(background_type, "list")){
    background = unlist(background)
  }

  if (pathway_definitions == "RaMP"){
    if (background_type == "list" & analyte_type == "metabolites") {
      backgrounddf <- getPathwayFromAnalyte(db = db, background,
                                            includeRaMPids = TRUE,
                                            namesOrIds = namesOrIds,
                                            include_smpdb = include_smpdb
                                            )
      print("Custom background specified, genes will be discarded")

    } else if (background_type=="file" & analyte_type == "metabolites") {
      userbkg <- utils::read.table(background, header=F)[,1]
      backgrounddf <- getPathwayFromAnalyte(db = db, analytes = userbkg,
                                            includeRaMPids = TRUE,
                                            namesOrIds = namesOrIds,
                                            include_smpdb = include_smpdb
                                            )
      print("Custom background specified, genes will be discarded")
    } else if (background_type == "biospecimen" & analyte_type == "metabolites") {
      biospecimen <- background
      if (biospecimen == "Adipose") {
        biospecimen <- "Adipose tissue"
      }

      # Get metabolites that belong to a specific biospecimen
      # query <- paste0(
      #   "SELECT analytehasontology.*, ontology.*, analytehaspathway.* from analytehasontology, ontology, analytehaspathway where ontology.commonName in ('",
      #   biospecimen,
      #   "') and ontology.rampOntologyId = analytehasontology.rampOntologyId and analytehasontology.rampCompoundId = analytehaspathway.rampId"
      # )

      # less data pull-back
      query <- paste0(
        "SELECT analytehaspathway.* from analytehasontology, ontology, analytehaspathway where ontology.commonName in ('",
        biospecimen,
        "') and analytehasontology.rampOntologyId = ontology.rampOntologyId and analytehasontology.rampCompoundId = analytehaspathway.rampId"
      )

      backgrounddf <- RaMP::runQuery(query, db)

      if (nrow(backgrounddf) == 0) {
        stop("Biospecimen background not found. Choices are 'Blood', 'Adipose', 'Heart', 'Urine', 'Brain', 'Liver', 'Kidney', 'Saliva', and 'Feces'")
      }

      # only keep the input metabolites (converted into pathwaydf in line above) that are in the biospecimen type specified
      pathwaydf <- with(pathwaydf, {
        pathwaydf %>%
          dplyr::filter(`rampId` %in% backgrounddf$rampId)
      })
      if (nrow(pathwaydf) == 0) {
        stop("There are no metabolites in your input that map to your selected biospecimen")
      }
    } else if (background_type == "database") {
    # do nothing, it's handled down below in if statements
    } else {
      stop("background_type was not specified correctly.  Please specify one of the following options: database, file, list, biospecimen")
    }
  }else{
    if (background_type == "list" & analyte_type == "metabolites") {
      backgrounddf <- getCustomPathwayFromAnalyte(analytes=background,
                                                  pathway_definitions=pathway_definitions,
                                                  analyte_type = analyte_type
                                                  )
      print("Custom background specified, genes will be discarded")
    } else if (background_type == "file" & analyte_type == "metabolites") {
      userbkg <- utils::read.table(background, header = F)[, 1]
      backgrounddf <- getCustomPathwayFromAnalyte(analytes=userbkg,
                                                  pathway_definitions=pathway_definitions,
                                                  analyte_type = analyte_type
                                                  )
      print("Custom background specified, genes will be discarded")
    }else if(analyte_type == "genes"){
      # do nothing, it's handled down below in if statements
    }else{
      stop("Only custom backgrounds are supported for custom pathway definitions. Please provide a 'list' or 'file' containing the analyte background")
    }
  }

  ## Check that all metabolites of interest are in the background
  if (background_type != "database") {
    if (length(setdiff(pathwaydf$rampId, backgrounddf$rampId) != 0)) {
      stop("All analytes in set of interest must also be in background")
    }
  }

  ## Initialize empty contingency table for later
  contingencyTb <- matrix(0, nrow = 2, ncol = 2)
  colnames(contingencyTb) <- c("In Pathway", "Not In Pathway")
  rownames(contingencyTb) <- c("All Metabolites", "User's Metabolites")
  ## Get pathway ids that contain the user analytes
  pid <- unique(pathwaydf$pathwayRampId)
  list_pid <- sapply(pid, shQuote)
  list_pid <- paste(list_pid, collapse = ",")

  # Get the total number of metabolites that are mapped to pathways in RaMP (that's the default background)
  # added conditional to not pull hmdb ids
  if(!include_smpdb){
    query <- "select distinct rampId, pathwaySource from analytehaspathway where pathwaySource != 'hmdb';"
  }else{
    query <- "select distinct rampId, pathwaySource from analytehaspathway;"
  }

  allids <- runQuery(sql = query, db = db)
  allids <- allids[!duplicated(allids), ]

  if ((analyte_type == "metabolites")) {
    totanalytes <- lapply(unique(allids$pathwaySource), function(x){
      return(allids %>% dplyr::filter(.data$pathwaySource==x) %>% nrow)
    })
    names(totanalytes) <- unique(allids$pathwaySource)
  } else if (analyte_type == "genes") {

    # for now we're using a fixed population size for genes
    # this can be enhanced to take a list of all measured genes
    # or use a subset of genes having pathway annotations within each source
    ## wiki_totanalytes <- react_totanalytes <- kegg_totanalytes <- total_genes
    totanalytes <- as.list(rep(total_genes,length(unique(allids$pathwaySource))))
    names(totanalytes) <- unique(allids$pathwaySource)
  } else {
    print("analyte_type must be 'metabolites' or 'genes'")
  }

  ## Input_RampIds is a table of all analytes included in pathways represented in the user set
  ## "User" refers to significant analytes
  input_RampIds <- buildFrequencyTables(inputdf = pathwaydf, pathway_definitions = pathway_definitions, analyte_type = analyte_type, db = db)
  if (is.null(input_RampIds)) {
    stop("Data doesn't exist")
  } else {
    segregated_id_list <- segregateDataBySource(input_RampIds = input_RampIds)
    if(analyte_type=="metabolites"){
      segregated_id_list <- segregated_id_list$metab
    }else{
      segregated_id_list <- segregated_id_list$gene
    }
  }

  # Loop through each pathway, build the contingency table, and calculate Fisher's Exact
  # test p-value
  pidCount <- 0
  pval <- oddsratio <- totinpath <- userinpath <- pidused <- c()
  for (i in pid) {
    ids_inpath <- pathwaydf[which(pathwaydf$pathwayRampId == i), "rampId"]

    pidCount <- pidCount + 1

    if (analyte_type == "metabolites") {
      # Check to make sure that this pathway does have metabolites
      if (length(grep("RAMP_C", ids_inpath)) == 0) {
        user_in_pathway <- 0
      } else {
        user_in_pathway <- length(unique(grep("RAMP_C", ids_inpath, value = TRUE)))
        if (background_type != "database") {
          ids_inpath_bg <- backgrounddf[which(backgrounddf$pathwayRampId == i), "rampId"]
          bg_in_pathway <- length(unique(grep("RAMP_C", ids_inpath_bg, value = TRUE)))
        }
      }

      if(pidCount == 1) {
        tot_user_analytes <- length(grep("RAMP_C", unique(pathwaydf$rampId)))
        if (background_type != "database") {
          tot_bg_analytes <- length(grep("RAMP_C", unique(backgrounddf$rampId)))
        }
      }

    } else { # if genes
      # Check to make sure that this pathway does have genes
      if (length(grep("RAMP_G", ids_inpath)) == 0) {
        user_in_pathway <- 0
      } else {
        user_in_pathway <- length(unique(grep("RAMP_G", ids_inpath, value = TRUE)))
      }

      if(pidCount == 1) {
        tot_user_analytes <- length(grep("RAMP_G", unique(pathwaydf$rampId)))
      }
    }
    pathway_index <- names(segregated_id_list)[which(sapply(segregated_id_list, function(x) i %in% x$pathwayRampId))]
    if(is.null(pathway_index) | is.na(pathway_index)){
      total_pathway_analytes = 0
    }else{
      total_pathway_analytes = totanalytes[[pathway_index]]
      tot_in_pathway <- segregated_id_list[[pathway_index]] %>% dplyr::filter(`pathwayRampId`==i) %>% dplyr::pull("Freq")
    }
    if (tot_in_pathway == 0 || user_in_pathway == 0) {
      pval <- c(pval, NA)
      oddsratio <- c(oddsratio, NA)
    } else {
      tot_out_pathway <- total_pathway_analytes - tot_in_pathway
      # fill the rest of the table out

      ## user_in_pathway <- length(unique(pathwaydf[which(pathwaydf$pathwayRampId==i),"rampId"]))
      if (background_type != "database") {
        bg_in_pathway <- length(unique(backgrounddf[which(backgrounddf$pathwayRampId == i), "rampId"]))
      }
      # EM - Corrected the following line that initially counted all input analytes without regard as to whether
      # whether they were genes or metabolites.
      # user_out_pathway <- length(unique(pathwaydf$rampId)) - user_in_pathway
      user_out_pathway <- tot_user_analytes - user_in_pathway

      if (background_type != "database") {
        bg_in_pathway <- length(unique(backgrounddf[which(backgrounddf$pathwayRampId == i), "rampId"]))
        bg_out_pathway <- tot_bg_analytes - bg_in_pathway
      }

      if (background_type == "database") {
        contingencyTb[1, 1] <- tot_in_pathway - user_in_pathway
      } else {
        contingencyTb[1, 1] <- bg_in_pathway
      }

      if (background_type == "database") {
        contingencyTb[1, 2] <- tot_out_pathway - user_out_pathway
      } else {
        contingencyTb[1, 2] <- bg_out_pathway
      }
      contingencyTb[2, 1] <- user_in_pathway
      contingencyTb[2, 2] <- user_out_pathway

      # Put the test into a try catch in case there's an issue, we'll have some details on the contingency matrix
      tryCatch(
        {
          result <- stats::fisher.test(contingencyTb, alternative = alternative)
        },
        error = function(e) {
          print(toString(e))
          print(i)
          print(contingencyTb)
          print(tot_in_pathway)
          print(tot_out_pathway)
          print(user_in_pathway)
          print(user_out_pathway)
          print(analyte_type)
          print(pathwaydf)
        }
      )
      pval <- c(pval, result$p.value)
      oddsratio <- c(oddsratio, result$estimate)
    } # End else tot_in_pathway is not zero

    userinpath <- c(userinpath, user_in_pathway)
    totinpath <- c(totinpath, tot_in_pathway)
    pidused <- c(pidused, i)
  } # end for loop

  print("")
  print(now - proc.time())
  print("")
  
  # Now run fisher's tests for all other pids (all pathways not covered in dataset)
  
  
  # only keep pathways that have >= minPathwaySize or < maxPathwaySize compounds
  keepers <- intersect(
    which(c(totinpath) >= minPathwaySize),
    which(c(totinpath) < maxPathwaySize)
  )
  
  # hist(totinpath,breaks=1000)
  print(paste0("Keeping ", length(keepers), " pathways"))
  # fdr <- stats::p.adjust(c(pval,pval2)[keepers],method="fdr")
  # holm <- stats::p.adjust(c(pval,pval2)[keepers],method="holm")
  print(paste0("Calculated p-values for ", length(pval), " pathways"))
  
  # format output (retrieve pathway name for each unique source id first
  out <- data.frame(
    pathwayRampId = pidused[keepers],
    Pval = pval[keepers], # FDR.Adjusted.Pval=fdr,
    # Holm.Adjusted.Pval=holm,
    Odds_Ratio = oddsratio[keepers],
    Num_In_Path = userinpath[keepers],
    Total_In_Path = totinpath[keepers]
  )
  
  # Remove duplicate pathways between wikipathways and reactome, only perfect overlaps
  # only make the dup list if it doesn't exist from a previous run in the session
  if( !exists('duplicate_pathways')) {
    duplicate_pathways <- findDuplicatePathways(db = db)
  }
  if (any(out$pathwayRampId %in% duplicate_pathways)) {
    out <- out[-which(out$pathwayRampId %in% duplicate_pathways), ]
  }

  out <- out[!duplicated(out), ]

  # for user is the output needed, based on what user input
  return(list(out, pathwaydf))
}

#' Do fisher test for only one pathway from search result
#' clicked on highchart
#' @param analytes a vector of analytes (genes or metabolites) that need to be searched
#' @param namesOrIds whether input is "names" or "ids" (default is "ids", must be the same for analytes and background)
#' @param total_genes number of genes analyzed in the experiment (e.g. background) (default is 20000, with assumption that analyte_type is "genes")
#' @param min_analyte if the number of analytes (gene or metabolite) in a pathway is
#' < min_analyte, do not report
#' @param alternative alternative hypothesis test passed on to fisher.test().  Options are two.sided, greater, or less (default is "less")
#' @param minPathwaySize the minimum number of pathway members (genes and metabolites) to include the pathway in the output (default = 5)
#' @param maxPathwaySize the maximum number of pathway memnbers (genes and metaboltes) to include the pathway in the output (default = 150)
#' @param includeRaMPids include internal RaMP identifiers (default is "FALSE")
#' @param background_type type of background that is input by the user.  Opions are "database" if user wants all
#' analytes from the RaMP database to be used as background; "file", if user wnats to input a file path with a list of background
#' analytes; "list", if user wants to input a vector of analyte IDs; "biospecimen", if user wants to specify a
#' biospecimen type (e.g. blood, adipose tissue, etc.) and have those biospecimen-specific analytes used.  For genes,
#' only the "database" option is used.
#' @param background background to be used for Fisher's tests.  If parameter 'background_type="database"', this parameter
#' is ignored (default="database"); if parameter 'background_type= "file"', then 'background' should be a file name (with
#' directory); if 'background_type="list"', then 'background' should be a vector of RaMP IDs; if 'backgroud_type="biospecimen"'
#' then users should specify one of the following: "Blood", "Adipose tissue", "Heart", "Urine", "Brain", "Liver", "Kidney",
#' "Saliva", and "Feces"
#' @param pathway_definitions If "RaMP" (default), use pathway definitions within RaMP-DB. Else, supply path to gmx file containing custom pathway definitions. GMX files are a tab-separated format that contain one analyte set per column, with the name of the set in the first row, and constituent analytes in subsequent rows. Please supply a .xls or .xlsx file. If supplying pathway definitions for genes and metabolites, ensure that metabolite definitions are on tab 1, and gene definitions are on tab2.
#' @param include_smpdb Include pathways from smpdb/hmdb in analysis. Excluded by default since definitions are highly redundant
#' @param db a RaMP databse object
#' @return a list containing two entries: [[1]] fishresults, a dataframe containing pathways with Fisher's p values
#' (raw and with FDR and Holm adjustment), number of user analytes in pathway, total number of analytes in pathway,
#' and pathway source ID/database. [[2]] analyte_type, a string specifying the type of analyte input into the function ("genes", "metabolites", or "both")
#' @examples
#' \dontrun{
#' fisher.results <- runCombinedFisherTest(analytes = c(
#' "hmdb:HMDB0000033",
#' "hmdb:HMDB0000052",
#' "hmdb:HMDB0000094",
#' "hmdb:HMDB0000161",
#' "hmdb:HMDB0000168",
#' "hmdb:HMDB0000191",
#' "hmdb:HMDB0000201",
#' "chemspider:10026",
#' "hmdb:HMDB0006059",
#' "Chemspider:6405",
#' "CAS:5657-19-2",
#' "hmdb:HMDB0002511",
#' "chemspider:20171375", "CAS:133-32-4",
#' "CAS:5746-90-7",
#' "CAS:477251-67-5",
#' "hmdb:HMDB0000695",
#' "chebi:15934",
#' "CAS:838-07-3",
#' "hmdb:HMDBP00789",
#' "hmdb:HMDBP00283",
#' "hmdb:HMDBP00284",
#' "hmdb:HMDBP00850"
#' ), db = rampDB )
#'
#' fisher.results <- runCombinedFisherTest(analytes = analyte.list, namesOrIds = "ids")
#' }
#' @export
runCombinedFisherTest <- function(
    analytes,
    namesOrIds = "ids",
    total_genes = 20000,
    min_analyte = 2,
    alternative = "less",
    minPathwaySize = 5,
    maxPathwaySize = 150,
    includeRaMPids = FALSE,
    background_type = "database",
    background = "database",
    pathway_definitions = "RaMP",
    include_smpdb = FALSE,
    db = RaMP()) {

  G <- M <- 0

  # Grab pathways that contain metabolites to run Fisher on metabolites
  # This will return all pathways that have at 8-120 metabolites/genes in them
  ## fishmetab <- pathwaydf[grep("RAMP_C_", pathwaydf$rampId), ]

  print("Running Fisher's tests on metabolites")

  outmetab <- runFisherTest(
    db = db,
    analytes = analytes,
    analyte_type = "metabolites",
    total_genes = total_genes,
    minPathwaySize = minPathwaySize,
    maxPathwaySize = maxPathwaySize,
    background_type = background_type,
    background = background,
    pathway_definitions = pathway_definitions,
    include_smpdb=include_smpdb,
    namesOrIds = namesOrIds
  )
  pathwaydf_metab <- outmetab[[2]]
  outmetab <- outmetab[[1]]
  if (!is.null(outmetab)) {
    M <- 1
  }


  ## Grab pathways that contain genes to run Fisher on genes
  ## fishgene <- pathwaydf[grep("RAMP_G_", pathwaydf$rampId), ]
  ## Genes are not evaluated if custom background is specified
  if (background_type == "database" & pathway_definitions == "RaMP") {
    print("Running Fisher's tests on genes")
    outgene <- runFisherTest(
      db = db,
      analytes = analytes,
      analyte_type = "genes",
      total_genes = total_genes,
      minPathwaySize = minPathwaySize,
      maxPathwaySize = maxPathwaySize,
      include_smpdb=include_smpdb,
      namesOrIds = namesOrIds
    )
    pathwaydf_gene <- outgene[[2]]
    outgene <- outgene[[1]]
  } else if (pathway_definitions != "RaMP") {
    outgene <- runFisherTest(
      db = db,
      analytes = analytes,
      analyte_type = "genes",
      total_genes = total_genes,
      minPathwaySize = minPathwaySize,
      maxPathwaySize = maxPathwaySize,
      background_type = background_type,
      background = background,
      pathway_definitions = pathway_definitions,
      include_smpdb=include_smpdb
    )
    pathwaydf_gene <- outgene[[2]]
    outgene <- outgene[[1]]
  } else {
    outgene <- NULL
    pathwaydf_gene <- NULL
  }
  # if no ids map to pathways, return an empty result.
  if ((is.null(pathwaydf_metab) || nrow(pathwaydf_metab) < 1) &&
    (is.null(pathwaydf_gene) || nrow(pathwaydf_gene) < 1)) {
    print("input analyte names did not map to RaMP pathways. Returning empty result.")
    return(list(fishresults = data.frame(), analyte_type = "none_empty_pathway_mapping", result_type = "pathway_enrichment"))
  }

  pathwaydf <- rbind(pathwaydf_metab, pathwaydf_gene)

  if (!is.null(outgene)) {
    G <- 1
  }
  if (is.null(outgene) & !is.null(outmetab)) {
    out <- outmetab
    fdr <- stats::p.adjust(out$Pval, method = "fdr")
    out <- cbind(out, fdr)
    colnames(out)[ncol(out)] <- "Pval_FDR"
    holm <- stats::p.adjust(out$Pval, method = "holm")
    out <- cbind(out, holm)
    colnames(out)[ncol(out)] <- "Pval_Holm"
    keepers <- which(out$Num_In_Path >= min_analyte)
    out2 <- merge(
      pathwaydf_metab[, c(
        "pathwayName", "pathwayRampId", "pathwayId",
        "pathwaySource"
      )],
      out[keepers, ],
      by = "pathwayRampId"
    )
  } else if (!is.null(outgene) && is.null(outmetab)) {
    out <- outgene
    fdr <- stats::p.adjust(out$Pval, method = "fdr")
    out <- cbind(out, fdr)
    colnames(out)[ncol(out)] <- "Pval_FDR"
    holm <- stats::p.adjust(out$Pval, method = "holm")
    out <- cbind(out, holm)
    colnames(out)[ncol(out)] <- "Pval_Holm"
    keepers <- which(out$Num_In_Path >= min_analyte)
    out2 <- merge(
      pathwaydf_gene[, c(
        "pathwayName", "pathwayRampId", "pathwayId",
        "pathwaySource"
      )],
      out[keepers, ],
      by = "pathwayRampId"
    )
  } else {
    # merge the results if both genes and metabolites were run
    G <- M <- 1
    allfish <- merge(outmetab, outgene,
      by = "pathwayRampId", all.x = T, all.y = T
    )
    colnames(allfish)[which(colnames(allfish) == "Pval.x")] <- "Pval_Metab"
    colnames(allfish)[which(colnames(allfish) == "Pval.y")] <- "Pval_Gene"
    colnames(allfish)[which(colnames(allfish) == "Odds_Ratio.x")] <- "Odds_Ratio_Metab"
    colnames(allfish)[which(colnames(allfish) == "Odds_Ratio.y")] <- "Odds_Ratio_Gene"
    colnames(allfish)[which(colnames(allfish) == "Total_In_Path.x")] <- "Total_In_Path_Metab"
    colnames(allfish)[which(colnames(allfish) == "Total_In_Path.y")] <- "Total_In_Path_Gene"
    colnames(allfish)[which(colnames(allfish) == "Num_In_Path.x")] <- "Num_In_Path_Metab"
    colnames(allfish)[which(colnames(allfish) == "Num_In_Path.y")] <- "Num_In_Path_Gene"

    # Calculate combined p-values for pathways that have both genes and metabolites
    gm <- intersect(which(!is.na(allfish$Pval_Metab)), which(!is.na(allfish$Pval_Gene)))
    combpval <- stats::pchisq(-2 * (log(allfish$Pval_Metab[gm]) + log(allfish$Pval_Gene[gm])),
      df = 2, lower.tail = FALSE
    )

    g <- which(is.na(allfish$Pval_Metab))
    gpval <- allfish$Pval_Gene[g]
    m <- which(is.na(allfish$Pval_Gene))
    mpval <- allfish$Pval_Metab[m]

    out <- rbind(allfish[gm, ], allfish[g, ], allfish[m, ])
    out <- cbind(out, c(combpval, gpval, mpval))
    colnames(out)[ncol(out)] <- "Pval_combined"
    fdr <- stats::p.adjust(out$Pval_combined, method = "fdr")
    out <- cbind(out, fdr)
    colnames(out)[ncol(out)] <- "Pval_combined_FDR"
    holm <- stats::p.adjust(out$Pval_combined, method = "holm")
    out <- cbind(out, holm)
    colnames(out)[ncol(out)] <- "Pval_combined_Holm"

    ## keepers <- intersect(
    ##   c(
    ##     which(out$Num_In_Path_Metab >= min_analyte),
    ##     which(is.na(out$Num_In_Path_Metab))
    ##   ),
    ##   c(
    ##     which(out$Num_In_Path_Gene >= min_analyte),
    ##     which(is.na(out$Num_In_Path_Gene))
    ##   )
    ## )

    keepers <- unique(
      c(
        which(out$Num_In_Path_Metab >= min_analyte),
        which(out$Num_In_Path_Gene >= min_analyte)
      )
    )

    # Now that p-values are calculated, only return pathways that are in the list
    # of pathways that contain user genes and metabolites
    ## pathwaydf <- getPathwayFromAnalyte(analytes,
    ##   includeRaMPids = TRUE,
    ##   namesOrIds = namesOrIds
    ##   )
    if(pathway_definitions!="RaMP"){
      pathwaydf$pathwayName = pathwaydf$pathwayRampId
    }
    out2 <- merge(
      pathwaydf[, c(
        "pathwayName", "pathwayRampId", "pathwayId",
        "pathwaySource"
      )],
      out[keepers, ],
      by = "pathwayRampId"
    )
  } # end merging when genes and metabolites were run
  out2 <- out2[!duplicated(out2), ]

  analyte_type <- c()
  if (G == 1 && M == 1) {
    analyte_type <- "both"
  } else if (G == 1 && M == 0) {
    analyte_type <- "genes"
  } else if (G == 0 && M == 1) {
    analyte_type <- "metabolites"
  }

  out2$analytes <- apply(out2, 1, function(x) {
    pathwayid <- x["pathwayRampId"]
    sigpathwaydf <- pathwaydf[which(pathwaydf$pathwayRampId == pathwayid), ]
    analytes <- sigpathwaydf[, "commonName"] %>%
      paste0(collapse = "|")  # KJK - changed semicolon to pipe, since lipids sometimes have ';' in the name
    return(analytes)
  })
  if (includeRaMPids) {
    return(list(fishresults = out2, analyte_type = analyte_type, result_type = "pathway_enrichment"))
  } else {
    return(list(fishresults = out2 %>% cleanup(), analyte_type = analyte_type, result_type = "pathway_enrichment"))
  }
}


#' Function that search analytes (gene or compounds)  or a list of analytes and
#' returns associated pathways
#'
#' @param analytes a vector of analytes (genes or metabolites) that need to be searched
#' @param find_synonym find all synonyms or just return same synonym (T/F)
#' @param namesOrIds whether input is "names" or "ids" (default is "ids")
#' @param includeRaMPids include internal RaMP identifiers (default is "FALSE")
#' @param include_smpdb Include pathways from smpdb/hmdb in analysis. Excluded by default since definitions are
#' highly redundant
#' @param minPathwaySize the minimum number of pathway members (genes and metabolites) to include the pathway in
#'  the output (default = 5)
#' @param maxPathwaySize the maximum number of pathway memnbers (genes and metaboltes) to include the pathway in
#' the output (default = 150)
#' @param db a RaMP database object
#' @importFrom rlang .data
#' @return a list contains all metabolites as name and pathway inside.
#' @examples
#' \dontrun{
#' getPathwayFromAnalyte(db = rampDB, namesOrIds="ids", c("ensembl:ENSG00000135679",
#'      "hmdb:HMDB0000064",
#'      "hmdb:HMDB0000148",
#'      "ensembl:ENSG00000141510"))
#'
#' getPathwayFromAnalyte(db = rampDB, namesOrIds = "names", c("Serotonin", "Bilirubin", "Urea"))
#' }
#' @export
getPathwayFromAnalyte <- function( analytes = "none",
                                  find_synonym = FALSE,
                                  namesOrIds = "ids",
                                  includeRaMPids = FALSE,
                                  include_smpdb = FALSE,
                                  minPathwaySize = 5,
                                  maxPathwaySize = 150,
                                  db = RaMP() ) {

  rampId <- pathwayRampId <- c()

  print("Starting getPathwayFromAnalyte()")
  if (is.null(analytes) || length(analytes) == 0) {
    warning("Input analyte list is NULL or empty. Aborting getPathwayFromAnalyte()")
    return(NULL)
  }

  if (!(namesOrIds %in% c("ids", "names"))) {
    warning(paste0(
      "namesOrIds must have a value in c('ids','names')\n",
      "Supplied namesOrIds falue = ('", namesOrIds, "')\nAborting getPathwayFromAnlyte()"
    ))
    return(NULL)
  }

  list_metabolite <- analytes
  list_metabolite <- sapply(list_metabolite, shQuote)
  list_metabolite <- paste(list_metabolite, collapse = ",")
  if (list_metabolite == "") {
    warning("Unable to retrieve metabolites")
    return(NULL)
  }

  isSQLite = .is_sqlite(x = db)

  if (namesOrIds == "ids") {

    print("Working on ID List...")

    sql <- paste0("select p.pathwayName, p.type as pathwaySource, p.sourceId as pathwayId, s.sourceId as inputId, group_concat(distinct s.commonName order by s.commonName separator '; ') as commonName, s.rampId, p.pathwayRampId from
                  source s,
                  analytehaspathway ap,
                  pathway p
                  where
                  s.sourceId in (", list_metabolite, ")
                  and
                  ap.rampId = s.rampId
                  and
                  p.pathwayRampId = ap.pathwayRampId
                  group by inputId, rampId, pathwayId, p.pathwayName, p.type, p.pathwayRampId
                  order by pathwayName asc")

    if(isSQLite) {
      sql <- paste0("select p.pathwayName, p.type as pathwaySource, p.sourceId as pathwayId, s.sourceId as inputId, group_concat(distinct s.commonName COLLATE NOCASE) as commonName, s.rampId, p.pathwayRampId from
                  source s,
                  analytehaspathway ap,
                  pathway p
                  where
                  s.sourceId in (", list_metabolite, ")
                  and
                  ap.rampId = s.rampId
                  and
                  p.pathwayRampId = ap.pathwayRampId
                  group by inputId, s.rampId, pathwayId, p.pathwayName, p.type, p.pathwayRampId
                  order by pathwayName asc")
    }

  } else {
    print("Working on analyte name list...")
    sql <- paste0(
      "select p.pathwayName, p.type as pathwaySource, p.sourceId as pathwayId, lower(asyn.Synonym) as commonName, group_concat(distinct s.sourceId order by s.sourceId separator '; ') as sourceIds, s.rampId, p.pathwayRampId
    from
    source s,
    analytesynonym asyn,
    analytehaspathway ap,
    pathway p
    where
    asyn.Synonym in (", list_metabolite, ")
    and
    s.rampId = asyn.rampId
    and
    ap.rampId = s.rampId
    and
    p.pathwayRampId = ap.pathwayRampId
    group by commonName, s.rampId, pathwayId, p.pathwayName, p.type, p.pathwayRampId
    order by pathwayName asc
  "
    )

    if(isSQLite) {
      sql <- paste0(
        "select p.pathwayName, p.type as pathwaySource, p.sourceId as pathwayId, lower(asyn.Synonym) as commonName, group_concat(distinct s.sourceId COLLATE NOCASE) as sourceIds, s.rampId, p.pathwayRampId
    from
    source s,
    analytesynonym asyn,
    analytehaspathway ap,
    pathway p
    where
    asyn.Synonym in (", list_metabolite, ")
    and
    s.rampId = asyn.rampId
    and
    ap.rampId = s.rampId
    and
    p.pathwayRampId = ap.pathwayRampId
    group by commonName, s.rampId, pathwayId, p.pathwayName, p.type, p.pathwayRampId
    order by pathwayName asc
  "
      )

    }

  }

  df2 <- runQuery(sql = sql, db = db)
  if(!include_smpdb){
    df2 <- df2 %>% dplyr::filter(.data$pathwaySource != "hmdb")
  }
  if (find_synonym && nrow(df2) > 0) {
    rampIds <- df2[, "rampId"]

    rampIds <- sapply(rampIds, shQuote)
    rampIds <- paste(rampIds, collapse = ",")

    sql <- paste0("select rampId as rampId, group_concat(distinct Synonym order by Synonym separator '; ')
     as synonyms from analytesynonym
     where rampId in (", rampIds, ") group by rampId")

    synonymsDf <- runQuery(sql = sql, db = db)

    if (nrow(synonymsDf) > 0) {
      df2 <- merge(df2, synonymsDf, by = "rampId")
    }
  }

  # filter by pathway size criteria

  df2 <- filterPathwaysByAnalyteCount(db = db, pathway_dataframe=df2, pathway_ramp_id_col_name = 'pathwayRampId', minPathwaySize = minPathwaySize, maxPathwaySize = maxPathwaySize)

  if (!includeRaMPids && nrow(df2) > 0) {
    df2 <- subset(df2, select = -c(rampId, pathwayRampId))
  }

  print("finished getPathwayFromAnalyte()")
  print(paste0("Found ", nrow(df2), " associated pathways."))

  return(df2)
}


#' Utility method supporting pathway analyses when file-based pathway lists are used rather than the ramp database.
#' @param analytes a vector of analytes (genes or metabolites) that need to be searched
#' @param pathway_definitions If "RaMP" (default), use pathway definitions within RaMP-DB. Else, supply path to gmx file containing custom pathway definitions. GMX files are a tab-separated format that contain one analyte set per column, with the name of the set in the first row, and constituent analytes in subsequent rows
#' @param analyte_type "genes" or "metabolites"
#' @return A pathwaydf compatible with runFisherTest
#' @importFrom rlang .data
#' @author Andrew Patt
getCustomPathwayFromAnalyte <- function(analytes, pathway_definitions, analyte_type) {
  print("Starting getCustomPathwayFromAnalyte()")
  tryCatch(
    {
      if (analyte_type == "metabolites") {
        pathway_definitions <- readxl::read_excel(pathway_definitions, sheet = 1)
      } else if (analyte_type == "genes") {
        pathway_definitions <- readxl::read_excel(pathway_definitions, sheet = 2)
      }
    },
    error = function(e) {
      print("Pathway file could not be found or is improperly formatted. Please supply path to GMX file for custom pathway definitions")
    }
  )

  pathwaydf <- data.frame("Analyte" = character(), "Pathway" = character())

  for (i in analytes) {
    for (j in 1:ncol(pathway_definitions)) {
      if (i %in% unlist(pathway_definitions[, j])) {
        pathwaydf <- rbind(pathwaydf, data.frame(
          "Analyte" = i,
          "Pathway" = colnames(pathway_definitions)[j]
        ))
      }
    }
  }

  print(paste0(
    "Found ", length(unique(pathwaydf$Analyte)),
    " out of ", length(unique(analytes)),
    " input analytes in pathway definitions file"
  ))

  pathwaydf$pathwayname <- pathwaydf$pathwayId <- pathwaydf$pathwayRampId <- pathwaydf$Pathway
  pathwaydf$inputId <- pathwaydf$commonName <- pathwaydf$Analyte
  if (analyte_type == "metabolites") {
    pathwaydf$rampId <- paste0("RAMP_C_", pathwaydf$inputId)
  } else if (analyte_type == "genes") {
    pathwaydf$rampId <- paste0("RAMP_G_", pathwaydf$inputId)
  }
  pathwaydf$pathwaySource <- "custom"
  pathwaydf <- subset(pathwaydf, select = -c(.data$Pathway, .data$Analyte))

  return(pathwaydf)
}



#' Perform fuzzy multiple linkage partitioning clustering on pathways
#' identified by Fisher's test
#'
#' @param fishers_df The full result object generated by runCombinedFisherTest
#' @param perc_analyte_overlap Minimum overlap for pathways to be considered similar
#' (Default = 0.5)
#' @param min_pathway_tocluster Minimum number of 'similar' pathways required to start
#' a cluster (medoid) (Default = 3)
#' @param perc_pathway_overlap Minimum overlap for clusters to merge (Default = 0.5)
#' @param db a RaMP database object
#'
#' @return list:[[1]] Pathway enrichment result with dataframe having a cluster assignment column added
#' [[2]] analyte type
#' [[3]] cluster assignment in the list form
#' @examples
#' \dontrun{
#' pathwaydf <- getPathwayFromAnalyte(c(
#' "ensembl:ENSG00000135679", "hmdb:HMDB0000064",
#' "hmdb:HMDB0000148", "ensembl:ENSG00000141510"
#' ))
#'
#' fisher.results <- runCombinedFisherTest(pathwaydf = pathwaydf)
#'
#' clustered.fisher.results <- findCluster(fisher.results)
#' }
#' @export
findCluster <- function(fishers_df, perc_analyte_overlap = 0.5,
                        min_pathway_tocluster = 2, perc_pathway_overlap = 0.5, db = RaMP()) {
  print("Clustering pathways...")

  if (perc_analyte_overlap <= 0 || perc_analyte_overlap >= 1 ||
    perc_pathway_overlap <= 0 || perc_pathway_overlap >= 1) {
    warning("No Clustering. perc_analyte_overlap and percent_pathway_overlap must bee in the range of (0,1), exclusive (not exactly 0 or 1).")
    return(fishers_df)
  }

  if (is.null(fishers_df$fishresults) || nrow(fishers_df$fishresults) < 1) {
    warning("The contained input pathway dataframe is empty (fishers_df$fishresults). Returning input result without clustering.")
    return(fishers_df)
  }

  analyte_type <- fishers_df$analyte_type
  fishers_df <- fishers_df$fishresults
  list_pathways <- fishers_df %>% dplyr::pull("pathwayId")
  list_pathways <- sapply(list_pathways, shQuote)
  list_pathways <- paste(list_pathways, collapse = ",")
  query <- paste0(
    "SELECT pathwayRampId, sourceId from pathway where sourceId in (",
    list_pathways,
    ")"
  )

  idkey <- runQuery(sql = query, db = db) %>%
    dplyr::rename("pathwayId" = "sourceId") ##  %>%
  ## dplyr::rename("rampId" = "pathwayRampId")

  rampToSource <- function(x) {
    out <- with(idkey, {
      idkey %>%
        dplyr::filter(`pathwayRampId` == x) %>%
        dplyr::pull("pathwayId")
    })
    return(out)
  }

  fishers_df <-
    fishers_df %>%
    dplyr::left_join(idkey, by = "pathwayId")
  if (nrow(fishers_df) == 0) {
    return(NULL)
  } else if (nrow(fishers_df) == 1) {
    fishers_df$cluster_assignment <- "Did not cluster"
    fishers_df$rampids <- fishers_df$pathwayRampId
    fishers_df$pathwayRampId <- NULL
    output <- list(fishresults = fishers_df, analyte_type = analyte_type, cluster_list = "Did not cluster")
    return(output)
  } else {
    # similarity_matrix_list<-loadOverlapMatrices()
    similarity_matrix_gene <- db@dbSummaryObjCache$genes_result
    similarity_matrix_analyte <- db@dbSummaryObjCache$analyte_result
    similarity_matrix_metab <- db@dbSummaryObjCache$metabolites_result
    if (analyte_type == "both") {
      # similarity_matrix = similarity_matrix_list[["analyte"]]
      similarity_matrix <- similarity_matrix_analyte
    } else if (analyte_type == "metabolites") {
      # similarity_matrix = similarity_matrix_list[["metab"]]
      similarity_matrix <- similarity_matrix_metab
    } else if (analyte_type == "genes") {
      # similarity_matrix = similarity_matrix_list[["gene"]]
      similarity_matrix <- similarity_matrix_gene
    } else {
      stop("analyte_type should be 'genes' or metabolites'")
    }
    pathway_list <- fishers_df[, "pathwayRampId"]

    pathway_indices <- match(pathway_list, rownames(similarity_matrix))

    if (length(which(is.na(pathway_indices))) > 0) {
      pathway_indices <- pathway_indices[-which(is.na(pathway_indices))]
    }

    pathway_matrix <- similarity_matrix[pathway_indices, pathway_indices]
    unmerged_clusters <- apply(pathway_matrix, 1, function(x) {
      # if(length(which(x>=perc_analyte_overlap))>(min_pathway_tocluster+1)){
      if (length(which(x >= perc_analyte_overlap)) > (min_pathway_tocluster - 1)) {
        return(colnames(pathway_matrix)[which(x >= perc_analyte_overlap)])
      } else {
        return(NA)
      }
    })
    # Remove the unmerged clusters
    if (length(which(is.na(unmerged_clusters))) > 0) {
      unmerged_clusters <- unmerged_clusters[-which(is.na(unmerged_clusters))]
    }

    if (length(unmerged_clusters) == 0) {
      # stop("No medoids found, make perc_analyte_overlap or min_pathway_tocluster smaller")
      cluster_list <- rep("Did not cluster", times = nrow(fishers_df))
    } else {
      # Evaluate similarity between clusters
      cluster_similarity <- matrix(0, ncol = length(unmerged_clusters), nrow = length(unmerged_clusters))
      for (i in 1:length(unmerged_clusters)) {
        for (j in 1:length(unmerged_clusters)) {
          cluster_similarity[i, j] <- length(intersect(unmerged_clusters[[i]], unmerged_clusters[[j]])) /
             ## length(unique(c(unmerged_clusters[[i]], unmerged_clusters[[j]])))
             min(c(length(unmerged_clusters[[i]]),length(unmerged_clusters[[j]])))
        }
      }
      colnames(cluster_similarity) <- rownames(cluster_similarity) <- names(unmerged_clusters)
      unmerged_cluster_similarity <- cluster_similarity

      cluster_list <- unmerged_clusters

      # Merge Clusters
      count <- 1
      while (length(which(cluster_similarity >= perc_pathway_overlap)) > nrow(cluster_similarity)) {
        cluster_similarity_mod <- cluster_similarity
        for (i in 1:nrow(cluster_similarity_mod)) {
          cluster_similarity_mod[i, i] <- 0
        }

        clusters_to_merge <- which(cluster_similarity_mod == max(cluster_similarity_mod), arr.ind = TRUE)
        clusters_to_merge <- unique(t(apply(clusters_to_merge, 1, sort)))
        for (i in 1:nrow(clusters_to_merge)) {

          # print(" ")
          # if(length(!is.na(cluster_list[[clusters_to_merge[i, 1]]])) > 1) {
          #   print(paste0("cluster list 1 is greater than 1 ", length(!is.na(cluster_list[[clusters_to_merge[i, 1]]]))))
          #   print(cluster_list[[clusters_to_merge[i, 1]]])
          # } else {
          #   print("cluster list 1 length is 1")
          # }
          #
          # if(length(!is.na(cluster_list[[clusters_to_merge[i, 1]]])) > 1) {
          #   print(paste0("cluster list 2 is greater than 1 ", length(!is.na(cluster_list[[clusters_to_merge[i, 2]]]))))
          #   print(cluster_list[[clusters_to_merge[i, 2]]])
          # } else {
          #   print("cluster list 2 length is 1")
          # }

          if (all(!is.na(cluster_list[[clusters_to_merge[i, 1]]])) && all(!is.na(cluster_list[[clusters_to_merge[i, 2]]]))) {
            cluster_list[[clusters_to_merge[i, 1]]] <- unique(unlist(cluster_list[c(clusters_to_merge[i, 1], clusters_to_merge[i, 2])]))
            cluster_list[[clusters_to_merge[i, 2]]] <- NA
          }
        }

        if (length(which(is.na(cluster_list))) > 0) {
          cluster_list <- cluster_list[-which(is.na(cluster_list))]
        }

        cluster_similarity <- matrix(0, ncol = length(cluster_list), nrow = length(cluster_list))
        for (i in 1:length(cluster_list)) {
          for (j in 1:length(cluster_list)) {
            cluster_similarity[i, j] <- length(intersect(cluster_list[[i]], cluster_list[[j]])) /
              length(unique(c(cluster_list[[i]], cluster_list[[j]])))
          }
        }

        if (nrow(cluster_similarity) == 1) {
          # stop("Clusters converged, use larger perc_pathway_overlap")
          # return(rep(1,times = nrow(fishers_df)))
          cluster_list <- rep("Did not cluster", times = nrow(fishers_df))
        }
        count <- count + 1
        if (count == length(unmerged_clusters) + 1) {
          stop("ERROR: while loop failed to terminate")
          # return(rep(1,times = nrow(fishers_df)))
          # cluster_list<-rep("Did not cluster",times = nrow(fishers_df))
        }
      }
      if (length(unique(cluster_list)) != 1) {
        colnames(cluster_similarity) <- rownames(cluster_similarity) <- paste0("cluster_", c(1:length(cluster_list)))
      }
    }
    # return(cluster_list)

    # Reformat cluster list to embed into results file
    rampids <- as.vector(fishers_df$pathwayRampId)
    # fishers_df$pathwayRampId<-NULL

    if (length(cluster_list) > 1) {
      cluster_assignment <- sapply(rampids, function(x) {
        pathway <- x
        clusters <- ""
        for (i in 1:length(cluster_list)) {
          if (pathway %in% cluster_list[[i]]) {
            clusters <- paste0(clusters, i, sep = ", ", collapse = ", ")
          }
        }
        if (clusters != "") {
          clusters <- substr(clusters, 1, nchar(clusters) - 2)
        } else {
          clusters <- "Did not cluster"
        }
        return(clusters)
      })
      fishers_df <- cbind(fishers_df, cluster_assignment)
    } else {
      fishers_df <- cbind(fishers_df, rep("Did not cluster", times = nrow(fishers_df)))
    }

    ## fishers_df$rampids <- rampids
    fishers_df <- cleanup(data = fishers_df)
    rownames(fishers_df) <- NULL

    ## Remove RaMP ids
    colnames(pathway_matrix) <- rownames(pathway_matrix) <- sapply(
      rownames(pathway_matrix),
      rampToSource
    )
    cluster_list <- lapply(cluster_list, function(x) {
      out <- sapply(x, rampToSource)
      names(out) <- NULL
      return(out)
    })
    names(cluster_list) <- paste("Cluster", 1:length(cluster_list))
    output <- list(fishresults = fishers_df, analyte_type = analyte_type, cluster_list = cluster_list, pathway_matrix = pathway_matrix)
    print("Finished clustering pathways...")
    return(output)
  }
}
