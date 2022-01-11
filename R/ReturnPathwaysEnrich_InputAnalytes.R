#' Do fisher test for only one pathway from search result
#' clicked on highchart
#' @param analytes a vector of analytes (genes or metabolites) that need to be searched
#' @param background optional vector of all metabolites detected in study. This will be used as the background for the Fisher's contingency table for metabolites. If left "database", all metabolites in RaMP originating in the pathway database of origin are used as background for testing
#' @param NameOrIds whether input is "names" or "ids" (default is "ids", must be the same for analytes and background)
#' @param biospecimen_background If "none", test all metabolites in RaMP or custom panel as Fisher's background. Else, use background for specific biospecimens. Choices are "Blood", "Adipose", "Heart", "Urine", "Brain", "Liver", "Kidney", "Saliva", and "Feces"
#' @param total_genes number of genes analyzed in the experiment (e.g. background) (default is 20000, with assumption that analyte_type is "genes")
#' @param analyte_type "metabolites" or "genes" (default is "metabolites")
#' @param MCall T/F if true, all pathways are used for multiple comparison corrections; if false, only pathways covering user analytes will be used (default is "T")
#' @param alternative alternative hypothesis test passed on to fisher.test().  Options are two.sided, greater, or less (default is "less")
#' @return a dataframe with columns containing pathway ID, fisher's p value, user analytes in pathway, and total analytes in pathway

runFisherTest <- function(analytes, background = "database",
                          biospecimen_background = "none", total_genes = 20000,
                          NameOrIds = "ids",
                          analyte_type = "metabolites",
                          MCall = T, alternative = "less") {
  now <- proc.time()
  print("Fisher Testing ......")
  pathwaydf <- getPathwayFromAnalyte(analytes,
    includeRaMPids = TRUE,
    NameOrIds = NameOrIds,
    find_synonym=FALSE
  )

  if (analyte_type == "metabolites") {
    pathwaydf <- pathwaydf[grep("RAMP_C_", pathwaydf$rampId), ]
  } else if (analyte_type == "genes") {
    pathwaydf <- pathwaydf[grep("RAMP_G_", pathwaydf$rampId), ]
  }

  if (nrow(pathwaydf) == 0) {
    return(NULL)
  }

  if (background != "database") {
    backgrounddf <- getPathwayFromAnalyte(background,
      includeRaMPids = TRUE,
      NameOrIds = NameOrIds
    )
  }

  if (biospecimen_background != "none") {
    if (biospecimen_background == "Adipose") {
      biospecimen_background <- "Adipose tissue"
    }
    query <- paste0(
      "SELECT analytehasontology.*, ontology.*, analytehaspathway.* from analytehasontology, ontology, analytehaspathway where ontology.commonName in ('",
      biospecimen_background,
      "') and ontology.rampOntologyIdLocation = analytehasontology.rampOntologyIdLocation and analytehasontology.rampCompoundId = analytehaspathway.rampId"
    )
    con <- connectToRaMP()
    backgrounddf <- DBI::dbGetQuery(con, query)
    DBI::dbDisconnect(con)
    if (nrow(backgrounddf) == 0) {
      stop("Biospecimen background not found. Choices are 'Blood', 'Adipose', 'Heart', 'Urine', 'Brain', 'Liver', 'Kidney', 'Saliva', and 'Feces'")
    }
    pathwaydf <- with(pathwaydf, {
      pathwaydf %>%
        dplyr::filter(rampId %in% backgrounddf$rampId)
    })
    if (nrow(pathwaydf) == 0) {
      stop("There are no metabolites in your input that map to your selected biospecimen")
    }
  }

  ## Check that all metabolites of interest are in the background
  if (background != "database") {
    if (length(setdiff(pathwaydf$rampId, backgrounddf$rampId) != 0)) {
      stop("All analytes in pathwaydf must also be in backgrounddf")
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
  query <- "select * from analytehaspathway"
  con <- connectToRaMP()
  allids <- DBI::dbGetQuery(con, query)

  # Close connection, then deduplicate id list
  DBI::dbDisconnect(con)
  allids <- allids[!duplicated(allids), ]

  if ((analyte_type == "metabolites")) {

    # JCB replaced these lines. Reducing to a source, extracting compound indices, then applying to the full set of rampIds
    # caused an error in the tally of analytes
    #
    # wiki_totanalytes <- length(unique(allids$rampId[grep("RAMP_C",allids[which(allids$pathwaySource=="wiki"),"rampId"])]))
    # react_totanalytes <- length(unique(allids$rampId[grep("RAMP_C",allids[which(allids$pathwaySource=="reactome"),"rampId"])]))
    # kegg_totanalytes <- length(unique(allids$rampId[grep("RAMP_C",allids[which(allids$pathwaySource=="kegg"),"rampId"])]))

    # first extract source-specific ids, then select for compound ids from the source-specific ids
    sourceIds <- allids[which(allids$pathwaySource == "wiki"), "rampId"]
    wiki_totanalytes <- length(unique(sourceIds[grep("RAMP_C", sourceIds)]))

    sourceIds <- allids[which(allids$pathwaySource == "reactome"), "rampId"]
    react_totanalytes <- length(unique(sourceIds[grep("RAMP_C", sourceIds)]))

    ## sourceIds <- allids[which(allids$pathwaySource=="kegg"),"rampId"]
    ## kegg_totanalytes <- length(unique(sourceIds[grep("RAMP_C",sourceIds)]))

    sourceIds <- allids[which(allids$pathwaySource == "kegg"), "rampId"]
    kegg_totanalytes <- length(unique(sourceIds[grep("RAMP_C", sourceIds)]))
  } else if (analyte_type == "genes") {

    # for now we're using a fixed population size for genes
    # this can be enhanced to take a list of all measured genes
    # or use a subset of genes having pathway annotations within each source
    wiki_totanalytes <- react_totanalytes <- kegg_totanalytes <- total_genes
  } else {
    print("analyte_type must be 'metabolites' or 'genes'")
  }

  ## Input_RampIds is a table of all analytes included in pathways represented in the user set
  ## "User" refers to significant analytes
  input_RampIds <- buildFrequencyTables(pathwaydf)

  if (is.null(input_RampIds)) {
    stop("Data doesn't exist")
  } else {
    segregated_id_list <- segregateDataBySource(input_RampIds)
  }

  # Loop through each pathway, build the contingency table, and calculate Fisher's Exact
  # test p-value
  pval <- totinpath <- userinpath <- pidused <- c()
  for (i in pid) {
    ids_inpath <- pathwaydf[which(pathwaydf$pathwayRampId == i), "rampId"]


    if (analyte_type == "metabolites") {
      # Check to make sure that this pathway does have metabolites
      if (length(grep("RAMP_C", ids_inpath)) == 0) {
        user_in_pathway <- 0
      } else {
        user_in_pathway <- length(unique(grep("RAMP_C", ids_inpath, value = TRUE)))
        if (background != "database") {
          ids_inpath_bg <- backgrounddf[which(backgrounddf$pathwayRampId == i), "rampId"]
          bg_in_pathway <- length(unique(grep("RAMP_C", ids_inpath_bg, value = TRUE)))
        }
      }
      inputkegg <- segregated_id_list[[1]][1][[1]]
      inputreact <- segregated_id_list[[1]][2][[1]]
      inputwiki <- segregated_id_list[[1]][3][[1]]
      tot_user_analytes <- length(grep("RAMP_C", unique(pathwaydf$rampId)))
      if (background != "database") {
        tot_bg_analytes <- length(grep("RAMP_C", unique(backgrounddf$rampId)))
      }
      ## if(background != "database"){
      ##     inputkegg_bg <- segregated_id_list_bg[[1]][1][[1]]
      ##     inputreact_bg <- segregated_id_list_bg[[1]][2][[1]]
      ##     inputwiki_bg <- segregated_id_list_bg[[1]][3][[1]]
      ## }
    } else {
      # Check to make sure that this pathway does have genes
      if (length(grep("RAMP_G", ids_inpath)) == 0) {
        user_in_pathway <- 0
      } else {
        user_in_pathway <- length(unique(grep("RAMP_G", ids_inpath, value = TRUE)))
      }
      inputkegg <- segregated_id_list[[2]][1][[1]]
      inputreact <- segregated_id_list[[2]][2][[1]]
      inputwiki <- segregated_id_list[[2]][3][[1]]
      tot_user_analytes <- length(grep("RAMP_G", unique(pathwaydf$rampId)))
      ## tot_bg_analytes <- length(grep("RAMP_G", unique(backgrounddf$rampId)))
    }

    if ((!is.na(inputkegg$pathwayRampId[1])) && i %in% inputkegg$pathwayRampId) {
      tot_in_pathway <- inputkegg[which(inputkegg[, "pathwayRampId"] == i), "Freq"]
      total_pathway_analytes <- kegg_totanalytes
    } else if ((!is.na(inputwiki$pathwayRampId[1])) && i %in% inputwiki$pathwayRampId) {
      tot_in_pathway <- inputwiki[which(inputwiki[, "pathwayRampId"] == i), "Freq"]
      total_pathway_analytes <- wiki_totanalytes
    } else if ((!is.na(inputreact$pathwayRampId[1])) && i %in% inputreact$pathwayRampId) {
      tot_in_pathway <- inputreact[which(inputreact[, "pathwayRampId"] == i), "Freq"]
      total_pathway_analytes <- react_totanalytes
    } else {
      tot_in_pathway <- 0
    }

    if (tot_in_pathway == 0 || user_in_pathway == 0) {
      pval <- c(pval, NA)
    } else {
      tot_out_pathway <- total_pathway_analytes - tot_in_pathway
      # fill the rest of the table out

      ## user_in_pathway <- length(unique(pathwaydf[which(pathwaydf$pathwayRampId==i),"rampId"]))
      if (background != "database") {
        bg_in_pathway <- length(unique(backgrounddf[which(backgrounddf$pathwayRampId == i), "rampId"]))
      }
      # EM - Corrected the following line that initially counted all input analytes without regard as to whether
      # whether they were genes or metabolites.
      # user_out_pathway <- length(unique(pathwaydf$rampId)) - user_in_pathway
      user_out_pathway <- tot_user_analytes - user_in_pathway

      if (background != "database") {
        bg_in_pathway <- length(unique(backgrounddf[which(backgrounddf$pathwayRampId == i), "rampId"]))
        bg_out_pathway <- tot_bg_analytes - bg_in_pathway
      }

      contingencyTb[1, 1] <- ifelse(background == "database",
        tot_in_pathway - user_in_pathway,
        bg_in_pathway
      )
      contingencyTb[1, 2] <- ifelse(background == "database",
        tot_out_pathway - user_out_pathway,
        bg_out_pathway
      )
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
    } # End else tot_in_pathway is not zero

    userinpath <- c(userinpath, user_in_pathway)
    totinpath <- c(totinpath, tot_in_pathway)
    pidused <- c(pidused, i)
  } # end for loop
  # Now run fisher's tests for all other pids (all pathways not covered in dataset)
  if (MCall == T) {
    # Now run fisher's tests for all other pids
    query <- "select distinct(pathwayRampId) from analytehaspathway where pathwaySource != 'hmdb';"
    con <- connectToRaMP()
    allpids <- DBI::dbGetQuery(con, query)
    DBI::dbDisconnect(con)
    pidstorun <- setdiff(allpids[, 1], pid)
    pidstorunlist <- sapply(pidstorun, shQuote)
    pidstorunlist <- paste(pidstorunlist, collapse = ",")

    query2 <- paste0(
      "select rampId,pathwayRampId from analytehaspathway where pathwayRampId in (",
      pidstorunlist, ")"
    )
    con <- connectToRaMP()
    restcids <- DBI::dbGetQuery(con, query2) # [[1]]
    DBI::dbDisconnect(con)

    query1 <- paste0("select rampId,pathwayRampId from analytehaspathway;")

    con <- connectToRaMP()
    allcids <- DBI::dbGetQuery(con, query1) # [[1]]
    DBI::dbDisconnect(con)

    # We're collecting p-values for all pathways, now those with no analyte support at all - JCB:?

    # calculating p-values for all other pathways
    kegg_metab <- kegg_metab
    kegg_gene <- kegg_gene
    wiki_metab <- wiki_metab
    wiki_gene <- wiki_gene
    reactome_metab <- reactome_metab
    reactome_gene <- reactome_gene
    hmdb_metab <- hmdb_metab
    hmdb_gene <- hmdb_gene
    count <- 1
    pval2 <- userinpath2 <- totinpath2 <- c()

    for (i in pidstorun) {
      if ((count %% 100) == 0) {
        print(paste0("Processed ", count))
      }
      count <- count + 1
      user_in_pathway <- 0
      if (analyte_type == "metabolites") {
        if (i %in% kegg_metab$pathwayRampId) {
          tot_in_pathway <- kegg_metab[which(kegg_metab[, "pathwayRampId"] == i), "Freq"]
          total_analytes <- kegg_totanalytes
        } else if (i %in% wiki_metab$pathwayRampId) {
          tot_in_pathway <- wiki_metab[which(wiki_metab[, "pathwayRampId"] == i), "Freq"]
          total_analytes <- wiki_totanalytes
        } else if (i %in% reactome_metab$pathwayRampId) {
          tot_in_pathway <- reactome_metab[which(reactome_metab[, "pathwayRampId"] == i), "Freq"]
          total_analytes <- react_totanalytes
        } else if (i %in% hmdb_metab$pathwayRampId) {
          tot_in_pathway <- hmdb_metab[which(hmdb_metab[, "pathwayRampId"] == i), "Freq"]
          total_analytes <- NULL
        } else {
          tot_in_pathway <- 0
          total_analytes <- NULL
        }
      } else {
        if (i %in% kegg_gene$pathwayRampId) {
          tot_in_pathway <- kegg_gene[which(kegg_gene[, "pathwayRampId"] == i), "Freq"]
          total_analytes <- kegg_totanalytes
        } else if (i %in% wiki_gene$pathwayRampId) {
          tot_in_pathway <- wiki_gene[which(wiki_gene[, "pathwayRampId"] == i), "Freq"]
          total_analytes <- wiki_totanalytes
        } else if (i %in% reactome_gene$pathwayRampId) {
          tot_in_pathway <- reactome_gene[which(reactome_gene[, "pathwayRampId"] == i), "Freq"]
          total_analytes <- react_totanalytes
        } else if (i %in% hmdb_gene$pathwayRampId) {
          tot_in_pathway <- hmdb_gene[which(hmdb_gene[, "pathwayRampId"] == i), "Freq"]
          total_analytes <- NULL
        } else {
          tot_in_pathway <- 0
          total_analytes <- NULL
        }
      }

      # Check that the pathway being considered has your analyte type, if not, move on

      if (is.null(total_analytes)) {
        next
      }
      tot_out_pathway <- total_analytes - tot_in_pathway
      # fill the rest of the table out

      # JCB: Another issue 10/7/2020
      # This line was used for user_out_pathway
      # This section of code is for all pathways that have no analyte support.
      user_out_pathway <- length(unique(pathwaydf$rampId))

      # This line was commented out in production *but* now we have total_analytes set properly
      # not sure why this line was changed.
      # user_out_pathway <- total_analytes - user_in_pathway

      contingencyTb[1, 1] <- tot_in_pathway - user_in_pathway
      contingencyTb[1, 2] <- tot_out_pathway - user_out_pathway
      contingencyTb[2, 1] <- user_in_pathway
      contingencyTb[2, 2] <- user_out_pathway
      # Added try catch
      tryCatch(
        {
          result <- stats::fisher.test(contingencyTb, alternative = alternative)
        },
        error = function(e) {
          print(toString(e))
          print(i)
          print(contingencyTb)
        }
      )

      pval2 <- c(pval2, result$p.value)
      userinpath2 <- c(userinpath2, user_in_pathway)
      totinpath2 <- c(totinpath2, tot_in_pathway)
      # pidused <- c(pidused,i)
    } # end for loop
    keepers <- intersect(
      which(c(totinpath, totinpath2) >= 8),
      which(c(totinpath, totinpath2) < 100)
    )

    print(paste0("Calculated p-values for ", length(c(pval, pval2)), " pathways"))
    out <- data.frame(
      pathwayRampId = c(pidused, pidstorun)[keepers],
      Pval = c(pval, pval2)[keepers], # FDR.Adjusted.Pval=fdr,
      # Holm.Adjusted.Pval=holm,
      Num_In_Path = c(userinpath, userinpath2)[keepers],
      Total_In_Path = c(totinpath, totinpath2)[keepers]
    )
  } # end if MCall is T and all pathways are being calculated, even ones that do not represent input analytes

  else {

    # only keep pathways that have > 8 or < 100 compounds
    keepers <- intersect(
      which(c(totinpath) >= 8),
      which(c(totinpath) < 100)
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
      Num_In_Path = userinpath[keepers],
      Total_In_Path = totinpath[keepers]
    )
  } # End else if MCall (when False)

  print(dim(out))
  out <- out[!duplicated(out), ]
  print(colnames(out))
  # for user is the output needed, based on what user input
  return(out)
}

#' Do fisher test for only one pathway from search result
#' clicked on highchart
#' @param analytes a vector of analytes (genes or metabolites) that need to be searched
#' @param background optional vector of all metabolites detected in study to be used as the background for the Fisher's contingency table for
#' enrichment.  If value is "database", all metabolites in the RaMP database will be used as background.  Default: "database"
#' @param biospecimen_background If "none", test all metabolites in RaMP or custom panel as Fisher's background. Else, use background for specific biospecimens. Choices are "Blood", "Adipose", "Heart", "Urine", "Brain", "Liver", "Kidney", "Saliva", and "Feces"
#' @param NameOrIds whether input is "names" or "ids" (default is "ids", must be the same for analytes and background)
#' @param total_genes number of genes analyzed in the experiment (e.g. background) (default is 20000, with assumption that analyte_type is "genes")
#' @param min_analyte if the number of analytes (gene or metabolite) in a pathway is
#' < min_analyte, do not report
#' @param MCall T/F if true, all pathways are used for multiple comparison corrections; if false, only pathways covering user analytes will be used (default is "T")
#' @param alternative alternative hypothesis test passed on to fisher.test().  Options are two.sided, greater, or less (default is "less")
#' @return a list containing two entries: [[1]] fishresults, a dataframe containing pathways with Fisher's p values (raw and with FDR and Holm adjustment), number of user analytes in pathway, total number of analytes in pathway, and pathway source ID/database. [[2]] analyte_type, a string specifying the type of analyte input into the function ("genes", "metabolites", or "both")
#' @examples
#' \dontrun{
#' pkg.globals <- setConnectionToRaMP(
#'   dbname = "ramp2", username = "root",
#'   conpass = "", host = "localhost"
#' )
#' pathwaydf <- getPathwayFromAnalyte(c("MDM2", "TP53", "glutamate", "creatinine"),
#'   NameOrIds = "names"
#' )
#' fisher.results <- runCombinedFisherTest(pathwaydf = pathwaydf)
#' }
#' @export
runCombinedFisherTest <- function(analytes, background = "database",
                                  biospecimen_background = "none",
                                  NameOrIds = "ids",
                                  total_genes = 20000,
                                  min_analyte = 2,
                                  MCall = T, alternative = "less") {
  G <- M <- 0

  # Grab pathways that contain metabolites to run Fisher on metabolites
  # This will return all pathways that have at 8-120 metabolites/genes in them
  ## fishmetab <- pathwaydf[grep("RAMP_C_", pathwaydf$rampId), ]

  print("Running Fisher's tests on metabolites")
  outmetab <- runFisherTest(
    analytes = analytes, background = background,
    biospecimen_background = biospecimen_background,
    analyte_type = "metabolites",
    total_genes = total_genes,
    MCall = MCall
  )
  if (!is.null(outmetab)) {
    M <- 1
  }

  # Grab pathways that contain genes to run Fisher on genes
  ## fishgene <- pathwaydf[grep("RAMP_G_", pathwaydf$rampId), ]
  print("Running Fisher's tests on genes")
  outgene <- runFisherTest(
    analytes = analytes,
    analyte_type = "genes",
    total_genes = total_genes,
    MCall = MCall
  )

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
    out2 <- merge(out[keepers, ],
      pathwaydf[, c(
        "pathwayName", "pathwayRampId", "pathwaysourceId",
        "pathwaysource"
      )],
      by = "pathwayRampId"
    )
  } else if (!is.null(outgene) & is.null(outmetab)) {
    out <- outgene
    fdr <- stats::p.adjust(out$Pval, method = "fdr")
    out <- cbind(out, fdr)
    colnames(out)[ncol(out)] <- "Pval_FDR"
    holm <- stats::p.adjust(out$Pval, method = "holm")
    out <- cbind(out, holm)
    colnames(out)[ncol(out)] <- "Pval_Holm"
    keepers <- which(out$Num_In_Path >= min_analyte)
    out2 <- merge(out[keepers, ],
      pathwaydf[, c(
        "pathwayName", "pathwayRampId", "pathwaysourceId",
        "pathwaysource"
      )],
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

    keepers <- intersect(
      c(
        which(out$Num_In_Path_Metab >= min_analyte),
        which(is.na(out$Num_In_Path_Metab))
      ),
      c(
        which(out$Num_In_Path_Gene >= min_analyte),
        which(is.na(out$Num_In_Path_Gene))
      )
    )


    # Now that p-values are calculated, only return pathways that are in the list
    # of pathways that contain user genes and metabolites
    pathwaydf <- getPathwayFromAnalyte(analytes,
      includeRaMPids = TRUE,
      NameOrIds = NameOrIds
    )
    out2 <- merge(out[keepers, ],
      pathwaydf[, c(
        "pathwayName", "pathwayRampId", "pathwaysourceId",
        "pathwaysource"
      )],
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
      paste0(collapse = ";")
    return(analytes)
  })
  return(list(fishresults = out2 %>% cleanup(), analyte_type = analyte_type))
}


#' Function that search analytes (gene or compounds)  or a list of analytes and
#' returns associated pathways
#'
#' @param analytes a vector of analytes (genes or metabolites) that need to be searched
#' @param find_synonym find all synonyms or just return same synonym (T/F)
#' @param NameOrIds whether input is "names" or "ids" (default is "ids")
#' @param includeRaMPids include internal RaMP identifiers (default is "FALSE")
#' @return a list contains all metabolites as name and pathway inside.
#' @examples
#' \dontrun{
#' pkg.globals <- setConnectionToRaMP(
#'   dbname = "ramp2", username = "root",
#'   conpass = "", host = "localhost"
#' )
#' mypath <- getPathwayFromAnalyte(analytes = c("2-hydroxyglutarate", "glutamate"))
#' }
#' @export
getPathwayFromAnalyte <- function(analytes = "none",
                                  find_synonym = FALSE,
                                  NameOrIds = "ids",
                                  includeRaMPids = FALSE) {
    now <- proc.time()
  if (length(analytes) == 1) {
    if (analytes == "none") {
      return(NULL)
    }
  }

  list_metabolite <- getRaMPInfoFromAnalytes(
    analytes = analytes, NameOrIds = NameOrIds,
    PathOrChem = "path"
  )
  list_metabolite <- unique(list_metabolite$rampId)
  list_metabolite <- sapply(list_metabolite, shQuote)
  list_metabolite <- paste(list_metabolite, collapse = ",")
  if (list_metabolite == "") {
    warning("Unable to retrieve metabolites")
    return(NULL)
  }
  # Now using the RaMP compound id, retrieve associated pathway ids
  query2 <- paste0(
    "select pathwayRampId,rampId from analytehaspathway where
                      rampId in (",
    list_metabolite, ");"
  )
  con <- connectToRaMP()
  # print(query2)
  df2 <- DBI::dbGetQuery(con, query2)
  DBI::dbDisconnect(con)
  pathid_list <- df2$pathwayRampId
  pathid_list <- sapply(pathid_list, shQuote)
  pathid_list <- paste(pathid_list, collapse = ",")
  # With pathway ids, retrieve pathway information
  if (pathid_list == "") {
    warning("The input list of analytes do not map to any pathways")
    return(NULL)
  }
  query3 <- paste0(
    "select pathwayName,sourceId as pathwaysourceId,type as pathwaysource,pathwayRampId from pathway where pathwayRampId in (",
    pathid_list, ");"
  )
  con <- connectToRaMP()
  df3 <- DBI::dbGetQuery(con, query3)
  DBI::dbDisconnect(con)
  # Format output
  mdf <- merge(df3, df2, all.x = T)
  # And with rampIds (list_metabolite), get common names when Ids are input
  if (NameOrIds == "ids") {
    list_analytes <- sapply(analytes, shQuote)
    list_analytes <- paste(list_analytes, collapse = ",")
    query4 <- paste0("select sourceId,commonName,rampId from source where sourceId in (", list_analytes, ");")
    con <- connectToRaMP()
    df4 <- DBI::dbGetQuery(con, query4)
    DBI::dbDisconnect(con)
    # convert latin1 encoding to UTF-8
    df4$commonName <- sapply(as.character(df4$commonName), function(x) {
      if (stringi::stri_enc_mark(x) == "native") {
        x <- iconv(x, "latin1", "UTF-8")
      } else {
        x
      }
    })
    mdf <- merge(mdf, df4, all.x = T, by.y = "rampId")
    mdf$commonName <- tolower(mdf$commonName)
  } else { # Just take on the name
    list_analytes <- sapply(mdf$rampId, shQuote)
    list_analytes <- paste(list_analytes, collapse = ",")
    query4 <- paste0("select sourceId,commonName,rampId from source where rampId in (", list_analytes, ");")
    con <- connectToRaMP()
    df4 <- DBI::dbGetQuery(con, query4)
    DBI::dbDisconnect(con)
                                        # convert latin1 encoding to UTF-8
    df4$commonName <- sapply(as.character(df4$commonName), function(x) {
      if (stringi::stri_enc_mark(x) == "native") {
        x <- iconv(x, "latin1", "UTF-8")
      } else {
        x
      }
    })
    mdf <- merge(mdf, df4, all.x = T, by = "rampId")
  }
  ## mdf <- merge(mdf, df4, all.x = T, by.y = "rampId")
  mdf$commonName <- tolower(mdf$commonName)
  ## if (find_synonym) {
  ##   mdf <- merge(mdf, "synonym", all.x = T, by = "rampId")
  ## } else
    if(NameOrIds == "names" & !find_synonym){
      mdf <- mdf[unlist(lapply(tolower(analytes),
                               function(x) which(tolower(mdf$commonName) %in% x))),]
  }
    
  out <- mdf[!duplicated(mdf), ]
  # For now, not returning HMDB pathways because they include the 30K
  # new pathways that are mainly drug and lipid pathways (need more proper
  # structural resolution matching)
  if (includeRaMPids) {
    return(out[
      which(out$pathwaysource != "hmdb"),
      c(
        "commonName", "pathwayName", "pathwaysource",
        "pathwaysourceId", "rampId", "pathwayRampId"
      )
    ])
  } else {
    out <- out[
      which(out$pathwaysource != "hmdb"),
      c(
        "commonName", "pathwayName", "pathwaysource",
        "pathwaysourceId", "rampId", "pathwayRampId"
      )
    ] %>%
        cleanup()
    out <- out[!duplicated(out), ]
    return(out)
  }
}

#' Perform fuzzy multiple linkage partitioning clustering on pathways
#' identified by Fisher's test
#'
#' @param fishers_df The data frame generated by runFisherTest
#' @param perc_analyte_overlap Minimum overlap for pathways to be considered similar
#' (Default = 0.5)
#' @param min_pathway_tocluster Minimum number of 'similar' pathways required to start
#' a cluster (medoid) (Default = 3)
#' @param perc_pathway_overlap Minimum overlap for clusters to merge (Default = 0.5)
#'
#' @return list:[[1]] Pathway enrichment result dataframe with cluster assignment column added
#' [[2]] analyte type
#' [[3]] cluster assignment in the list form
#' @examples
#' \dontrun{
#' pkg.globals <- setConnectionToRaMP(
#'   dbname = "ramp2", username = "root",
#'   conpass = "", host = "localhost"
#' )
#' pathwaydf <- getPathwayFromAnalyte(c(
#'   "ensembl:ENSG00000135679", "hmdb:HMDB0000064",
#'   "hmdb:HMDB0000148", "ensembl:ENSG00000141510"
#' ))
#' fisher.results <- runCombinedFisherTest(pathwaydf = pathwaydf)
#' filtered.fisher.results <- FilterFishersResults(fisher.results, p_holmadj_cutoff = 0.05)
#' filteredclust.fisher.results <- findCluster(filtered.fisher.results)
#' }
#' @export
findCluster <- function(fishers_df, perc_analyte_overlap = 0.5,
                        min_pathway_tocluster = 2, perc_pathway_overlap = 0.5) {
  if (perc_analyte_overlap <= 0 || perc_analyte_overlap >= 1 ||
    perc_pathway_overlap <= 0 || perc_pathway_overlap >= 1) {
    return(NULL)
  }
  analyte_type <- fishers_df$analyte_type
  fishers_df <- fishers_df$fishresults
  list_pathways <- fishers_df %>% dplyr::pull("pathwaysourceId")
  list_pathways <- sapply(list_pathways, shQuote)
  list_pathways <- paste(list_pathways, collapse = ",")
  query <- paste0(
    "SELECT pathwayRampId, sourceId from pathway where sourceId in (",
    list_pathways,
    ")"
  )
  con <- connectToRaMP()
  idkey <- DBI::dbGetQuery(con, query) %>%
    dplyr::rename("pathwaysourceId" = "sourceId") ##  %>%
  ## dplyr::rename("rampId" = "pathwayRampId")

  rampToSource <- function(x) {
    out <- with(idkey, {
      idkey %>%
        dplyr::filter(pathwayRampId == x) %>%
        dplyr::pull("pathwaysourceId")
    })
    return(out)
  }

  DBI::dbDisconnect(con)

  fishers_df <-
    fishers_df %>%
    dplyr::left_join(idkey, by = "pathwaysourceId")
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
    similarity_matrix_gene <- genes_result
    similarity_matrix_analyte <- analyte_result
    similarity_matrix_metab <- metabolites_result
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
            length(unique(c(unmerged_clusters[[i]], unmerged_clusters[[j]])))
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
          if (!is.na(cluster_list[[clusters_to_merge[i, 1]]]) && !is.na(cluster_list[[clusters_to_merge[i, 2]]])) {
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
    fishers_df <- cleanup(fishers_df)
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
    return(output)
  }
}

#' Filter pathways by p-value cutoff for display and clustering
#' @param fishers_df The data frame generated by runFisherTest
#' @param p_holmadj_cutoff return pathways where Holm adjusted pvalues are < p_holmadj_cutoff
#' @param p_fdradj_cutoff return pathways where FDR adjusted pvalues are < p_fdradj_cutoff
#' @return list:[[1]]Dataframe with pathway enrichment results, only significant pathways
#' [[2]]analyte type
#' @examples
#' \dontrun{
#' pkg.globals <- setConnectionToRaMP(
#'   dbname = "ramp2", username = "root",
#'   conpass = "", host = "localhost"
#' )
#' pathwaydf <- getPathwayFromAnalyte(c("MDM2", "TP53", "glutamate", "creatinine"),
#'   NameOrIds = "names"
#' )
#' fisher.results <- runCombinedFisherTest(pathwaydf = pathwaydf)
#' filtered.fisher.results <- FilterFishersResults(fisher.results, p_holmadj_cutoff = 0.05)
#' }
#' @export
FilterFishersResults <- function(fishers_df, p_holmadj_cutoff = 1,
                                 p_fdradj_cutoff = 1) {

  # Check to see whether the output is from ORA performed on genes and metabolites
  # or genes or metabolites
  analyte_type <- fishers_df$analyte_type
  fishers_df <- fishers_df$fishresults

  if (length(grep("Pval_combined", colnames(fishers_df))) == 0) {
    if (p_holmadj_cutoff != 1) {
      return(list(fishresults = fishers_df[which(fishers_df[, "Pval_Holm"] <=
        p_holmadj_cutoff), ], analyte_type = analyte_type))
    } else if (p_fdradj_cutoff != 1) {
      return(list(fishresults = fishers_df[which(fishers_df[, "Pval_FDR"] <=
        p_fdradj_cutoff), ], analyte_type = analyte_type))
    } else {
      stop("Please set a cutoff for Holm Adjusted pvalues
			(p_holmadj_cutoff paramter) or FDR Adjusted pvalues
			(p_fdradj_cutoff)")
    }
  } else { # ORA was performed on both genes and metabolites:
    if (p_holmadj_cutoff != 1) {
      return(list(fishresults = fishers_df[which(fishers_df[, "Pval_combined_Holm"] <=
        p_holmadj_cutoff), ], analyte_type = analyte_type))
    } else if (p_fdradj_cutoff != 1) {
      return(list(fishresults = fishers_df[which(fishers_df[, "Pval_combined_FDR"] <=
        p_fdradj_cutoff), ], analyte_type = analyte_type))
    } else {
      stop("Please set a cutoff for Holm Adjusted pvalues
                        (p_holmadj_cutoff paramter) or FDR Adjusted pvalues
                        (p_fdradj_cutoff)")
    }
  }
}
