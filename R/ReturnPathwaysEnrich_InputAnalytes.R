#' Do fisher test for only one pathway from search result
#'
#' @param analytes a vector of analytes (genes or metabolites) that need to be searched
#' @param namesOrIds whether input is "names" or "ids" (default is "ids", must be the same for analytes and background)
#' @param totalGenes number of genes analyzed in the experiment (e.g. background) (default is 20000, with assumption that analyteType is "genes")
#' @param analyteType "metabolites" or "genes" (default is "metabolites")
#' @param alternative alternative hypothesis test passed on to fisher.test().  Options are two.sided, greater, or less (default is "less")
#' @param minPathwaySize the minimum number of pathway members (genes and metabolites) to include the pathway in the output (default = 5)
#' @param maxPathwaySize the maximum number of pathway memnbers (genes and metaboltes) to include the pathway in the output (default = 150)
#' @param backgroundType type of background that is input by the user.  Opions are "database" if user wants all
#' analytes from the RaMP database will be used; "file", if user wants to input a file with a list of background
#' analytes; "list", if user wants to input a vector of analyte IDs; "biospecimen", if user wants to specify a
#' biospecimen type (e.g. blood, adipose tissue, etc.) and have those biospecimen-specific analytes used.  For genes,
#' only the "database" option is used.
#' @param background background to be used for Fisher's tests.  If parameter 'backgroundType="database"', this parameter
#' is ignored (default="database"); if parameter 'backgroundType= "file"', then 'background' should be a file name (with
#' directory); if 'backgroundType="list"', then 'background' should be a vector of RaMP IDs; if 'backgroud_type="biospecimen"'
#' then users should specify one of the following: "Blood", "Adipose tissue", "Heart", "Urine", "Brain", "Liver", "Kidney",
#' "Saliva", and "Feces"
#' @param pathwayDefinitions If "RaMP" (default), use pathway definitions within RaMP-DB. Else, supply path to gmx file containing custom pathway definitions. GMX files are a tab-separated format that contain one analyte set per column, with the name of the set in the first row, and constituent analytes in subsequent rows
#' @param includeSmpdb Include pathways from smpdb/hmdb in analysis. Excluded by default since definitions are highly redundant
#' @param db a RaMP databse object
#' @return a dataframe with columns containing pathway ID, fisher's p value, user analytes in pathway, and total analytes in pathway
#' @noRd
runFisherTest <- function(analytes,
                          totalGenes = 20000,
                          namesOrIds = "ids",
                          analyteType = "metabolites",
                          alternative = "less",
                          minPathwaySize = 5, maxPathwaySize = 150,
                          backgroundType = "database", background = "database",
                          pathwayDefinitions = "RaMP", includeSmpdb = FALSE,
                          db = RaMP()) {
  if(analyteType == "genes"){
    backgroundType = "database"
    print("Using database background for genes")
  }
  now <- proc.time()
  print("Fisher Testing ......")
  if (pathwayDefinitions != "RaMP") {
    pathwaydf <- getCustomPathwayFromAnalyte(analytes=analytes,
      pathwayDefinitions=pathwayDefinitions,
      analyteType = analyteType
      )
  } else {
    pathwaydf <- getPathwayFromAnalyte(db = db, analytes = analytes,
    includeRaMPids = TRUE,
      namesOrIds = namesOrIds,
      findSynonym = FALSE,
      includeSmpdb = includeSmpdb
    )
  }

  pathwayRampId <- rampId <- c()

  pathwayRampId <- rampId <- c()

  if (analyteType == "metabolites") {
    pathwaydf <- pathwaydf[grep("RAMP_C_", pathwaydf$rampId), ]
  } else if (analyteType == "genes") {
    pathwaydf <- pathwaydf[grep("RAMP_G_", pathwaydf$rampId), ]
  }


  # moved this check until we determine if we have analytes of a given type.
  if (nrow(pathwaydf) == 0) {
    return(NULL)
  }

#  if(class(backgroundType)=="list"){
   if(is(backgroundType, "list")){
    background = unlist(background)
  }

  if (pathwayDefinitions == "RaMP"){
    if (backgroundType == "list" & analyteType == "metabolites") {
      backgrounddf <- getPathwayFromAnalyte(db = db, background,
                                            includeRaMPids = TRUE,
                                            namesOrIds = namesOrIds,
                                            includeSmpdb = includeSmpdb
                                            )
      print("Custom background specified, genes will be discarded")

    } else if (backgroundType=="file" & analyteType == "metabolites") {
      userbkg <- utils::read.table(background, header=F)[,1]
      backgrounddf <- getPathwayFromAnalyte(db = db, analytes = userbkg,
                                            includeRaMPids = TRUE,
                                            namesOrIds = namesOrIds,
                                            includeSmpdb = includeSmpdb
                                            )
      print("Custom background specified, genes will be discarded")
    } else if (backgroundType == "biospecimen" & analyteType == "metabolites") {
      biospecimen <- background
      if (biospecimen == "Adipose") {
        biospecimen <- "Adipose tissue"
      }

      backgrounddf <- db@api$getAnalytePathwaysWithOntology(biospecimen=biospecimen)

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
    } else if (backgroundType == "database") {
    # do nothing, it's handled down below in if statements
    } else {
      stop("backgroundType was not specified correctly.  Please specify one of the following options: database, file, list, biospecimen")
    }
  }else{
    if (backgroundType == "list" & analyteType == "metabolites") {
      backgrounddf <- getCustomPathwayFromAnalyte(analytes=background,
                                                  pathwayDefinitions=pathwayDefinitions,
                                                  analyteType = analyteType
                                                  )
      print("Custom background specified, genes will be discarded")
    } else if (backgroundType == "file" & analyteType == "metabolites") {
      userbkg <- utils::read.table(background, header = F)[, 1]
      backgrounddf <- getCustomPathwayFromAnalyte(analytes=userbkg,
                                                  pathwayDefinitions=pathwayDefinitions,
                                                  analyteType = analyteType
                                                  )
      print("Custom background specified, genes will be discarded")
    }else if(analyteType == "genes"){
      # do nothing, it's handled down below in if statements
    }else{
      stop("Only custom backgrounds are supported for custom pathway definitions. Please provide a 'list' or 'file' containing the analyte background")
    }
  }

  ## Check that all metabolites of interest are in the background
  if (backgroundType != "database") {
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

  # Get the total number of metabolites that are mapped to pathways in RaMP (that's the default background)
  allids <- db@api$getRampIDsAndSourcesForPathways(includeSMPDB = includeSmpdb)
  allids <- allids[!duplicated(allids), ]

  if ((analyteType == "metabolites")) {
    totanalytes <- lapply(unique(allids$pathwaySource), function(x){
      return(allids %>% dplyr::filter(.data$pathwaySource==x) %>% nrow)
    })
    names(totanalytes) <- unique(allids$pathwaySource)
  } else if (analyteType == "genes") {

    # for now we're using a fixed population size for genes
    # this can be enhanced to take a list of all measured genes
    # or use a subset of genes having pathway annotations within each source
    ## wiki_totanalytes <- react_totanalytes <- kegg_totanalytes <- totalGenes
    totanalytes <- as.list(rep(totalGenes,length(unique(allids$pathwaySource))))
    names(totanalytes) <- unique(allids$pathwaySource)
  } else {
    print("analyteType must be 'metabolites' or 'genes'")
  }

  ## Input_RampIds is a table of all analytes included in pathways represented in the user set
  ## "User" refers to significant analytes
  input_RampIds <- buildFrequencyTables(inputdf = pathwaydf, pathwayDefinitions = pathwayDefinitions, analyteType = analyteType, db = db)
  if (is.null(input_RampIds)) {
    stop("Data doesn't exist")
  } else {
    segregated_id_list <- segregateDataBySource(inputRampIds = input_RampIds)
    if(analyteType=="metabolites"){
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

    if (analyteType == "metabolites") {
      # Check to make sure that this pathway does have metabolites
      if (length(grep("RAMP_C", ids_inpath)) == 0) {
        user_in_pathway <- 0
      } else {
        user_in_pathway <- length(unique(grep("RAMP_C", ids_inpath, value = TRUE)))
        if (backgroundType != "database") {
          ids_inpath_bg <- backgrounddf[which(backgrounddf$pathwayRampId == i), "rampId"]
          bg_in_pathway <- length(unique(grep("RAMP_C", ids_inpath_bg, value = TRUE)))
        }
      }

      if(pidCount == 1) {
        tot_user_analytes <- length(grep("RAMP_C", unique(pathwaydf$rampId)))
        if (backgroundType != "database") {
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
      if (backgroundType != "database") {
        bg_in_pathway <- length(unique(backgrounddf[which(backgrounddf$pathwayRampId == i), "rampId"]))
      }
      # EM - Corrected the following line that initially counted all input analytes without regard as to whether
      # whether they were genes or metabolites.
      # user_out_pathway <- length(unique(pathwaydf$rampId)) - user_in_pathway
      user_out_pathway <- tot_user_analytes - user_in_pathway

      if (backgroundType != "database") {
        bg_in_pathway <- length(unique(backgrounddf[which(backgrounddf$pathwayRampId == i), "rampId"]))
        bg_out_pathway <- tot_bg_analytes - bg_in_pathway
      }

      if (backgroundType == "database") {
        contingencyTb[1, 1] <- tot_in_pathway - user_in_pathway
      } else {
        contingencyTb[1, 1] <- bg_in_pathway
      }

      if (backgroundType == "database") {
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
          print(analyteType)
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
    OR = oddsratio[keepers],
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
#' @param totalGenes number of genes analyzed in the experiment (e.g. background) (default is 20000, with assumption that analyteType is "genes")
#' @param minAnalyte if the number of analytes (gene or metabolite) in a pathway is
#' < minAnalyte, do not report
#' @param alternative alternative hypothesis test passed on to fisher.test().  Options are two.sided, greater, or less (default is "less")
#' @param minPathwaySize the minimum number of pathway members (genes and metabolites) to include the pathway in the output (default = 5)
#' @param maxPathwaySize the maximum number of pathway memnbers (genes and metaboltes) to include the pathway in the output (default = 150)
#' @param includeRaMPids include internal RaMP identifiers (default is "FALSE")
#' @param backgroundType type of background that is input by the user.  Opions are "database" if user wants all
#' analytes from the RaMP database to be used as background; "file", if user wnats to input a file path with a list of background
#' analytes; "list", if user wants to input a vector of analyte IDs; "biospecimen", if user wants to specify a
#' biospecimen type (e.g. blood, adipose tissue, etc.) and have those biospecimen-specific analytes used.  For genes,
#' only the "database" option is used.
#' @param background background to be used for Fisher's tests.  If parameter 'backgroundType="database"', this parameter
#' is ignored (default="database"); if parameter 'backgroundType= "file"', then 'background' should be a file name (with
#' directory); if 'backgroundType="list"', then 'background' should be a vector of RaMP IDs; if 'backgroud_type="biospecimen"'
#' then users should specify one of the following: "Blood", "Adipose tissue", "Heart", "Urine", "Brain", "Liver", "Kidney",
#' "Saliva", and "Feces"
#' @param pathwayDefinitions If "RaMP" (default), use pathway definitions within RaMP-DB. Else, supply path to gmx file containing custom pathway definitions. GMX files are a tab-separated format that contain one analyte set per column, with the name of the set in the first row, and constituent analytes in subsequent rows. Please supply a .xls or .xlsx file. If supplying pathway definitions for genes and metabolites, ensure that metabolite definitions are on tab 1, and gene definitions are on tab2.
#' @param includeSmpdb Include pathways from smpdb/hmdb in analysis. Excluded by default since definitions are highly redundant
#' @param db a RaMP databse object
#' @return a list containing two entries: [[1]] fishresults, a dataframe containing pathways with Fisher's p values
#' (raw and with FDR and Holm adjustment), number of user analytes in pathway, total number of analytes in pathway,
#' and pathway source ID/database. [[2]] analyteType, a string specifying the type of analyte input into the function ("genes", "metabolites", or "both")
#' @examples
#' \dontrun{
#' fisher.results <- runEnrichPathways(analytes = c(
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
#' fisher.results <- runEnrichPathways(analytes = analyte.list, namesOrIds = "ids")
#' }
#' @export
runEnrichPathways <- function(
    analytes,
    namesOrIds = "ids",
    totalGenes = 20000,
    minAnalyte = 2,
    alternative = "less",
    minPathwaySize = 5,
    maxPathwaySize = 150,
    includeRaMPids = FALSE,
    backgroundType = "database",
    background = "database",
    pathwayDefinitions = "RaMP",
    includeSmpdb = FALSE,
    db = RaMP()) {

  G <- M <- 0

  # Grab pathways that contain metabolites to run Fisher on metabolites
  # This will return all pathways that have at 8-120 metabolites/genes in them
  ## fishmetab <- pathwaydf[grep("RAMP_C_", pathwaydf$rampId), ]

  print("Running Fisher's tests on metabolites")

  outmetab <- runFisherTest(
    db = db,
    analytes = analytes,
    analyteType = "metabolites",
    totalGenes = totalGenes,
    minPathwaySize = minPathwaySize,
    maxPathwaySize = maxPathwaySize,
    backgroundType = backgroundType,
    background = background,
    pathwayDefinitions = pathwayDefinitions,
    includeSmpdb=includeSmpdb,
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
  if (backgroundType == "database" & pathwayDefinitions == "RaMP") {
    print("Running Fisher's tests on genes")
    outgene <- runFisherTest(
      db = db,
      analytes = analytes,
      analyteType = "genes",
      totalGenes = totalGenes,
      minPathwaySize = minPathwaySize,
      maxPathwaySize = maxPathwaySize,
      includeSmpdb=includeSmpdb,
      namesOrIds = namesOrIds
    )
    pathwaydf_gene <- outgene[[2]]
    outgene <- outgene[[1]]
  } else if (pathwayDefinitions != "RaMP") {
    outgene <- runFisherTest(
      db = db,
      analytes = analytes,
      analyteType = "genes",
      totalGenes = totalGenes,
      minPathwaySize = minPathwaySize,
      maxPathwaySize = maxPathwaySize,
      backgroundType = backgroundType,
      background = background,
      pathwayDefinitions = pathwayDefinitions,
      includeSmpdb=includeSmpdb
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
    return(list(fishresults = data.frame(), analyteType = "none_empty_pathway_mapping", result_type = "pathway_enrichment"))
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
    keepers <- which(out$Num_In_Path >= minAnalyte)
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
    keepers <- which(out$Num_In_Path >= minAnalyte)
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
    colnames(allfish)[which(colnames(allfish) == "OR.x")] <- "OR_Metab"
    colnames(allfish)[which(colnames(allfish) == "OR.y")] <- "OR_Gene"
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
    colnames(out)[ncol(out)] <- "Pval_FDR"
    holm <- stats::p.adjust(out$Pval_combined, method = "holm")
    out <- cbind(out, holm)
    colnames(out)[ncol(out)] <- "Pval_Holm"

    ## keepers <- intersect(
    ##   c(
    ##     which(out$Num_In_Path_Metab >= minAnalyte),
    ##     which(is.na(out$Num_In_Path_Metab))
    ##   ),
    ##   c(
    ##     which(out$Num_In_Path_Gene >= minAnalyte),
    ##     which(is.na(out$Num_In_Path_Gene))
    ##   )
    ## )

    keepers <- unique(
      c(
        which(out$Num_In_Path_Metab >= minAnalyte),
        which(out$Num_In_Path_Gene >= minAnalyte)
      )
    )

    # Now that p-values are calculated, only return pathways that are in the list
    # of pathways that contain user genes and metabolites
    ## pathwaydf <- getPathwayFromAnalyte(analytes,
    ##   includeRaMPids = TRUE,
    ##   namesOrIds = namesOrIds
    ##   )
    if(pathwayDefinitions!="RaMP"){
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

  analyteType <- c()
  if (G == 1 && M == 1) {
    analyteType <- "both"
  } else if (G == 1 && M == 0) {
    analyteType <- "genes"
  } else if (G == 0 && M == 1) {
    analyteType <- "metabolites"
  }

  out2$analytes <- apply(out2, 1, function(x) {
    pathwayid <- x["pathwayRampId"]
    sigpathwaydf <- pathwaydf[which(pathwaydf$pathwayRampId == pathwayid), ]
    analytes <- sigpathwaydf[, "commonName"] %>%
      paste0(collapse = "|")  # KJK - changed semicolon to pipe, since lipids sometimes have ';' in the name
    return(analytes)
  })
  if (includeRaMPids) {
    return(list(fishresults = out2, analyteType = analyteType, result_type = "pathway_enrichment"))
  } else {
    return(list(fishresults = out2 %>% cleanup(), analyteType = analyteType, result_type = "pathway_enrichment"))
  }
}


#' Function that search analytes (gene or compounds)  or a list of analytes and
#' returns associated pathways
#'
#' @param analytes a vector of analytes (genes or metabolites) that need to be searched
#' @param findSynonym find all synonyms or just return same synonym (T/F)
#' @param namesOrIds whether input is "names" or "ids" (default is "ids")
#' @param includeRaMPids include internal RaMP identifiers (default is "FALSE")
#' @param includeSmpdb Include pathways from smpdb/hmdb in analysis. Excluded by default since definitions are
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
                                  findSynonym = FALSE,
                                  namesOrIds = "ids",
                                  includeRaMPids = FALSE,
                                  includeSmpdb = FALSE,
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

  df2 <- db@api$getPathwaysForAnalytes(analytes = analytes, namesOrIds = namesOrIds, includeSMPDB = includeSmpdb)

  if (findSynonym && nrow(df2) > 0) {
    rampIds <- df2[, "rampId"]
    synonymsDf <- db@api$getSynonymsForAnalyte(rampIds = rampIds)

    if (nrow(synonymsDf) > 0) {
      df2 <- merge(df2, synonymsDf, by = "rampId")
    }
  }

  # filter by pathway size criteria

  df2 <- filterPathwaysByAnalyteCount(db = db, pathwayDataframe=df2, pathwayRampIdColName = 'pathwayRampId', minPathwaySize = minPathwaySize, maxPathwaySize = maxPathwaySize)

  if (!includeRaMPids && nrow(df2) > 0) {
    df2 <- subset(df2, select = -c(rampId, pathwayRampId))
  }

  print("finished getPathwayFromAnalyte()")
  print(paste0("Found ", nrow(df2), " associated pathways."))

  return(df2)
}


#' Utility method supporting pathway analyses when file-based pathway lists are used rather than the ramp database.
#' @param analytes a vector of analytes (genes or metabolites) that need to be searched
#' @param pathwayDefinitions If "RaMP" (default), use pathway definitions within RaMP-DB. Else, supply path to gmx file containing custom pathway definitions. GMX files are a tab-separated format that contain one analyte set per column, with the name of the set in the first row, and constituent analytes in subsequent rows
#' @param analyteType "genes" or "metabolites"
#' @return A pathwaydf compatible with runFisherTest
#' @importFrom rlang .data
#' @author Andrew Patt
#' @noRd
getCustomPathwayFromAnalyte <- function(analytes, pathwayDefinitions, analyteType) {
  print("Starting getCustomPathwayFromAnalyte()")
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

  pathwaydf <- data.frame("Analyte" = character(), "Pathway" = character())

  for (i in analytes) {
    for (j in 1:ncol(pathwayDefinitions)) {
      if (i %in% unlist(pathwayDefinitions[, j])) {
        pathwaydf <- rbind(pathwaydf, data.frame(
          "Analyte" = i,
          "Pathway" = colnames(pathwayDefinitions)[j]
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
  if (analyteType == "metabolites") {
    pathwaydf$rampId <- paste0("RAMP_C_", pathwaydf$inputId)
  } else if (analyteType == "genes") {
    pathwaydf$rampId <- paste0("RAMP_G_", pathwaydf$inputId)
  }
  pathwaydf$pathwaySource <- "custom"
  pathwaydf <- subset(pathwaydf, select = -c(.data$Pathway, .data$Analyte))

  return(pathwaydf)
}



#' Perform fuzzy multiple linkage partitioning clustering on pathways
#' identified by Fisher's test
#'
#' @param fishersDf The full result object generated by runEnrichPathways
#' @param percAnalyteOverlap Minimum overlap for pathways to be considered similar
#' (Default = 0.5)
#' @param minPathwayToCluster Minimum number of 'similar' pathways required to start
#' a cluster (medoid) (Default = 3)
#' @param percPathwayOverlap Minimum overlap for clusters to merge (Default = 0.5)
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
#' fisher.results <- runEnrichPathways(pathwaydf = pathwaydf)
#'
#' clustered.fisher.results <- findCluster(fisher.results)
#' }
#' @export
findCluster <- function(fishersDf, percAnalyteOverlap = 0.5,
                        minPathwayToCluster = 2, percPathwayOverlap = 0.5, db = RaMP()) {
  print("Clustering pathways...")

  if (percAnalyteOverlap <= 0 || percAnalyteOverlap >= 1 ||
    percPathwayOverlap <= 0 || percPathwayOverlap >= 1) {
    warning("No Clustering. percAnalyteOverlap and percent_pathway_overlap must bee in the range of (0,1), exclusive (not exactly 0 or 1).")
    return(fishersDf)
  }

  if (is.null(fishersDf$fishresults) || nrow(fishersDf$fishresults) < 1) {
    warning("The contained input pathway dataframe is empty (fishersDf$fishresults). Returning input result without clustering.")
    return(fishersDf)
  }

  analyteType <- fishersDf$analyteType
  fishersDf <- fishersDf$fishresults
  list_pathways <- fishersDf %>% dplyr::pull("pathwayId")

  idkey <- db@api$getPathwayFromSourceId(pathwaySourceIDs = list_pathways) %>%
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

  fishersDf <-
    fishersDf %>%
    dplyr::left_join(idkey, by = "pathwayId")
  if (nrow(fishersDf) == 0) {
    return(NULL)
  } else if (nrow(fishersDf) == 1) {
    fishersDf$cluster_assignment <- "Did not cluster"
    fishersDf$rampids <- fishersDf$pathwayRampId
    fishersDf$pathwayRampId <- NULL
    output <- list(fishresults = fishersDf, analyteType = analyteType, cluster_list = "Did not cluster")
    return(output)
  } else {
    # similarity_matrix_list<-loadOverlapMatrices()
    similarity_matrix_gene <- db@dbSummaryObjCache$genes_result
    similarity_matrix_analyte <- db@dbSummaryObjCache$analyte_result
    similarity_matrix_metab <- db@dbSummaryObjCache$metabolites_result
    if (analyteType == "both") {
      # similarity_matrix = similarity_matrix_list[["analyte"]]
      similarity_matrix <- similarity_matrix_analyte
    } else if (analyteType == "metabolites") {
      # similarity_matrix = similarity_matrix_list[["metab"]]
      similarity_matrix <- similarity_matrix_metab
    } else if (analyteType == "genes") {
      # similarity_matrix = similarity_matrix_list[["gene"]]
      similarity_matrix <- similarity_matrix_gene
    } else {
      stop("analyteType should be 'genes' or metabolites'")
    }
    pathway_list <- fishersDf[, "pathwayRampId"]

    pathway_indices <- match(pathway_list, rownames(similarity_matrix))

    if (length(which(is.na(pathway_indices))) > 0) {
      pathway_indices <- pathway_indices[-which(is.na(pathway_indices))]
    }

    pathway_matrix <- similarity_matrix[pathway_indices, pathway_indices]
    unmerged_clusters <- apply(pathway_matrix, 1, function(x) {
      # if(length(which(x>=percAnalyteOverlap))>(minPathwayToCluster+1)){
      if (length(which(x >= percAnalyteOverlap)) > (minPathwayToCluster - 1)) {
        return(colnames(pathway_matrix)[which(x >= percAnalyteOverlap)])
      } else {
        return(NA)
      }
    })
    # Remove the unmerged clusters
    if (length(which(is.na(unmerged_clusters))) > 0) {
      unmerged_clusters <- unmerged_clusters[-which(is.na(unmerged_clusters))]
    }

    if (length(unmerged_clusters) == 0) {
      # stop("No medoids found, make percAnalyteOverlap or minPathwayToCluster smaller")
      cluster_list <- rep("Did not cluster", times = nrow(fishersDf))
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
      while (length(which(cluster_similarity >= percPathwayOverlap)) > nrow(cluster_similarity)) {
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
          # stop("Clusters converged, use larger percPathwayOverlap")
          # return(rep(1,times = nrow(fishersDf)))
          cluster_list <- rep("Did not cluster", times = nrow(fishersDf))
        }
        count <- count + 1
        if (count == length(unmerged_clusters) + 1) {
          stop("ERROR: while loop failed to terminate")
          # return(rep(1,times = nrow(fishersDf)))
          # cluster_list<-rep("Did not cluster",times = nrow(fishersDf))
        }
      }
      if (length(unique(cluster_list)) != 1) {
        colnames(cluster_similarity) <- rownames(cluster_similarity) <- paste0("cluster_", c(1:length(cluster_list)))
      }
    }
    # return(cluster_list)

    # Reformat cluster list to embed into results file
    rampids <- as.vector(fishersDf$pathwayRampId)
    # fishersDf$pathwayRampId<-NULL

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
      fishersDf <- cbind(fishersDf, cluster_assignment)
    } else {
      fishersDf <- cbind(fishersDf, rep("Did not cluster", times = nrow(fishersDf)))
    }

    ## fishersDf$rampids <- rampids
    fishersDf <- cleanup(data = fishersDf)
    rownames(fishersDf) <- NULL

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
    output <- list(fishresults = fishersDf, analyteType = analyteType, cluster_list = cluster_list, pathway_matrix = pathway_matrix)
    print("Finished clustering pathways...")
    return(output)
  }
}
