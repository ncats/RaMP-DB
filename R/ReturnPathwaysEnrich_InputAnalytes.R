#' Do fisher test for only one pathway from search result
#' clicked on highchart
#' @param analytes a vector of analytes (genes or metabolites) that need to be searched
#' @param NameOrIds whether input is "names" or "ids" (default is "ids", must be the same for analytes and background)
#' @param total_genes number of genes analyzed in the experiment (e.g. background) (default is 20000, with assumption that analyte_type is "genes")
#' @param analyte_type "metabolites" or "genes" (default is "metabolites")
#' @param MCall T/F if true, all pathways are used for multiple comparison corrections; if false, only pathways covering user analytes will be used (default is "F")
#' @param alternative alternative hypothesis test passed on to fisher.test().  Options are two.sided, greater, or less (default is "less")
#' @param min_path_size the minimum number of pathway members (genes and metabolites) to include the pathway in the output (default = 5)
#' @param max_path_size the maximum number of pathway memnbers (genes and metaboltes) to include the pathway in the output (default = 150)
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
#' @return a dataframe with columns containing pathway ID, fisher's p value, user analytes in pathway, and total analytes in pathway

runFisherTest <- function(analytes,
                          total_genes = 20000,
                          NameOrIds = "ids",
                          analyte_type = "metabolites",
                          MCall = F, alternative = "less", min_path_size = 5, max_path_size = 150,
                          background_type = "database", background = "database") {
  now <- proc.time()
  print("Fisher Testing ......")
  pathwaydf <- getPathwayFromAnalyte(analytes,
    includeRaMPids = TRUE,
    NameOrIds = NameOrIds,
    find_synonym = FALSE
  )

  if (analyte_type == "metabolites") {
    pathwaydf <- pathwaydf[grep("RAMP_C_", pathwaydf$rampId), ]
  } else if (analyte_type == "genes") {
    pathwaydf <- pathwaydf[grep("RAMP_G_", pathwaydf$rampId), ]
  }

  # moved this check until we determine if we have analytes of a given type.
  if (nrow(pathwaydf) == 0) {
    return(NULL)
  }

  if (class(background_type) == "list") {
    background <- unlist(background)
  }

  if (background_type == "list") {
    backgrounddf <- getPathwayFromAnalyte(background,
      includeRaMPids = TRUE,
      NameOrIds = NameOrIds
    )
    print("Custom background specified, genes will be discarded")
  } else if (background_type == "file") {
    userbkg <- read.table(background, header = F)[, 1]
    backgrounddf <- getPathwayFromAnalyte(userbkg,
      includeRaMPids = TRUE,
      NameOrIds = NameOrIds
    )
    print("Custom background specified, genes will be discarded")
  } else if (background_type == "biospecimen") {
    biospecimen <- background
    if (biospecimen == "Adipose") {
      biospecimen <- "Adipose tissue"
    }
    # Get metabolites that belong to a specific biospecimen
    query <- paste0(
      "SELECT analytehasontology.*, ontology.*, analytehaspathway.* from analytehasontology, ontology, analytehaspathway where ontology.commonName in ('",
      biospecimen,
      "') and ontology.rampOntologyId = analytehasontology.rampOntologyId and analytehasontology.rampCompoundId = analytehaspathway.rampId"
    )
    con <- connectToRaMP()
    backgrounddf <- RMariaDB::dbGetQuery(con, query)
    RMariaDB::dbDisconnect(con)
    if (nrow(backgrounddf) == 0) {
      stop("Biospecimen background not found. Choices are 'Blood', 'Adipose', 'Heart', 'Urine', 'Brain', 'Liver', 'Kidney', 'Saliva', and 'Feces'")
    }
    # only keep the input metabolites (converted into pathwaydf in line above) that are in the biospecimen type specified
    pathwaydf <- with(pathwaydf, {
      pathwaydf %>%
        dplyr::filter(rampId %in% backgrounddf$rampId)
    })
    if (nrow(pathwaydf) == 0) {
      stop("There are no metabolites in your input that map to your selected biospecimen")
    }
  } else if (background_type == "database") {
    # do nothing, it's handled down below in if statements
  } else {
    stop("background_type was not specified correctly.  Please specify one of the following options: database, file, list, biospecimen")
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
  query <- "select * from analytehaspathway"
  con <- connectToRaMP()
  allids <- RMariaDB::dbGetQuery(con, query)

  # Close connection, then deduplicate id list
  RMariaDB::dbDisconnect(con)
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
        if (background_type != "database") {
          ids_inpath_bg <- backgrounddf[which(backgrounddf$pathwayRampId == i), "rampId"]
          bg_in_pathway <- length(unique(grep("RAMP_C", ids_inpath_bg, value = TRUE)))
        }
      }
      inputkegg <- segregated_id_list[[1]][1][[1]]
      inputreact <- segregated_id_list[[1]][2][[1]]
      inputwiki <- segregated_id_list[[1]][3][[1]]
      tot_user_analytes <- length(grep("RAMP_C", unique(pathwaydf$rampId)))
      if (background_type != "database") {
        tot_bg_analytes <- length(grep("RAMP_C", unique(backgrounddf$rampId)))
      }
      ## if(background != "database"){
      ##     inputkegg_bg <- segregated_id_list_bg[[1]][1][[1]]
      ##     inputreact_bg <- segregated_id_list_bg[[1]][2][[1]]
      ##     inputwiki_bg <- segregated_id_list_bg[[1]][3][[1]]
      ## }
    } else { # if genes
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

      ## contingencyTb[1, 1] <- ifelse(background_type == "database",
      ##   tot_in_pathway - user_in_pathway,
      ##   bg_in_pathway
      ## )
      ## contingencyTb[1, 2] <- ifelse(background_type == "database",
      ##   tot_out_pathway - user_out_pathway,
      ##   bg_out_pathway
      ## )
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
    allpids <- RMariaDB::dbGetQuery(con, query)
    RMariaDB::dbDisconnect(con)
    pidstorun <- setdiff(allpids[, 1], pid)
    pidstorunlist <- sapply(pidstorun, shQuote)
    pidstorunlist <- paste(pidstorunlist, collapse = ",")

    query2 <- paste0(
      "select rampId,pathwayRampId from analytehaspathway where pathwayRampId in (",
      pidstorunlist, ")"
    )
    con <- connectToRaMP()
    restcids <- RMariaDB::dbGetQuery(con, query2) # [[1]]
    RMariaDB::dbDisconnect(con)

    query1 <- paste0("select rampId,pathwayRampId from analytehaspathway;")

    con <- connectToRaMP()
    allcids <- RMariaDB::dbGetQuery(con, query1) # [[1]]
    RMariaDB::dbDisconnect(con)

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
      which(c(totinpath, totinpath2) >= min_path_size),
      which(c(totinpath, totinpath2) < max_path_size)
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

    # only keep pathways that have >= min_path_size or < max_path_size compounds
    keepers <- intersect(
      which(c(totinpath) >= min_path_size),
      which(c(totinpath) < max_path_size)
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
  # Remove duplicate pathways between wikipathways and KEGG
  duplicate_pathways = find_duplicate_pathways()

  out <- out[-which(out$pathwayRampId %in% duplicate_pathways),]
  
  out <- out[!duplicated(out), ]
  print(dim(out))
  print(colnames(out))
  # for user is the output needed, based on what user input
  return(list(out, pathwaydf))
}

#' Do fisher test with user-supplied custom pathway definitions
#' @param analytes a vector of analytes (genes or metabolites) that need to be searched
#' @param alternative alternative hypothesis test passed on to fisher.test().  Options are two.sided, greater, or less (default is "less")
#' @param background background to be used for Fisher's tests.  If parameter 'background_type="database"', this parameter
#' is ignored (default="database"); if parameter 'background_type= "file"', then 'background' should be a file name (with
#' directory); if 'background_type="list"', then 'background' should be a vector of RaMP IDs
#' @param background_type type of background that is input by the user.  Options are "database" if user wants all
#' analytes from the RaMP database will be used; "file", if user wants to input a file with a list of background
#' analytes; "list", if user wants to input a vector of analyte IDs
#' @param pathways (String) Path to gmx file containing custom pathway definitions. GMX files are a tab-separated format that contain one analyte set per column, with the name of the set in the first row, and constituent analytes in subsequent rows
#' @return a dataframe with columns containing pathway ID, fisher's p value, user analytes in pathway, and total analytes in pathway

runCustomFisherTest <- function(analytes, background = "kegg", background_type = "database", pathways,
                                alternative = "less") {
  browser()
  tryCatch({
    pathway_definitions <- read.csv(pathways, sep = "\t")
  },
  error = function(e) {
    print("Pathway file could not be found or is improperly formatted. Please supply path to GMX file for custom pathway definitions")
  })
  now <- proc.time()
  print("Fisher Testing ......")

  pathwaydf <- data.frame(Analyte=character(),Pathway=character())

  for(i in analytes){
    for(j in 1:ncol(pathway_definitions)){
      if(i %in% unlist(pathway_definitions[,j])){
        pathwaydf<-rbind(pathwaydf, data.frame(Analyte=i,Pathway=colnames(pathway_definitions)[j]))
      }
    }
  }

  print(paste0("Found ", length(unique(pathwaydf$Analyte)),
               " out of ", length(unique(analytes)),
               " input analytes in pathway definitions file"))
  
  # moved this check until we determine if we have analytes of a given type.
  if (nrow(pathwaydf) == 0) {
    return(NULL)
  }
  
  if (class(background_type) == "list") {
    background <- unlist(background)
  }
  
  if (background_type == "list") {
    backgrounddf <- data.frame(Analyte=character(),Pathway=character())
    
    for(i in background){
      for(j in 1:ncol(pathway_definitions)){
        if(i %in% unlist(pathway_definitions[,j])){
          pathwaydf<-rbind(pathwaydf, data.frame(Analyte=i,Pathway=colnames(pathway_definitions)[j]))
        }
      }
    }
  } else if (background_type == "file") {
    userbkg <- read.table(background, header = F)[, 1]
    for(i in userbkg){
      for(j in 1:ncol(pathway_definitions)){
        if(i %in% unlist(pathway_definitions[,j])){
          pathwaydf<-rbind(pathwaydf, data.frame(Analyte=i,Pathway=colnames(pathway_definitions)[j]))
        }
      }
    }
  } else if (background_type == "database") {
    # do nothing, it's handled down below in if statements
  } else {
    stop("background_type was not specified correctly.  Please specify one of the following options: database, file, list (biospecimen is not supported for custom pathway definitions)")
  }


  ## Check that all metabolites of interest are in the background
  if (background_type != "database") {
    if (length(setdiff(pathwaydf$Analyte, backgrounddf$Analyte) != 0)) {
      stop("All analytes in set of interest must also be in background")
    }
  }

  ## Initialize empty contingency table for later
  contingencyTb <- matrix(0, nrow = 2, ncol = 2)
  colnames(contingencyTb) <- c("In Pathway", "Not In Pathway")
  rownames(contingencyTb) <- c("All Metabolites", "User's Metabolites")
  ## Get pathway ids that contain the user analytes
  pid <- unique(pathwaydf$Pathway)
  list_pid <- sapply(pid, shQuote)
  list_pid <- paste(list_pid, collapse = ",")

  # Get the total number of metabolites that are mapped to pathways in RaMP (that's the default background)
  query <- "select * from analytehaspathway"
  con <- connectToRaMP()
  allids <- RMariaDB::dbGetQuery(con, query)

  # Close connection, then deduplicate id list
  RMariaDB::dbDisconnect(con)
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
  # Left off here 5/11/22 (ACP)
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
        if (background_type != "database") {
          ids_inpath_bg <- backgrounddf[which(backgrounddf$pathwayRampId == i), "rampId"]
          bg_in_pathway <- length(unique(grep("RAMP_C", ids_inpath_bg, value = TRUE)))
        }
      }
      inputkegg <- segregated_id_list[[1]][1][[1]]
      inputreact <- segregated_id_list[[1]][2][[1]]
      inputwiki <- segregated_id_list[[1]][3][[1]]
      tot_user_analytes <- length(grep("RAMP_C", unique(pathwaydf$rampId)))
      if (background_type != "database") {
        tot_bg_analytes <- length(grep("RAMP_C", unique(backgrounddf$rampId)))
      }
      ## if(background != "database"){
      ##     inputkegg_bg <- segregated_id_list_bg[[1]][1][[1]]
      ##     inputreact_bg <- segregated_id_list_bg[[1]][2][[1]]
      ##     inputwiki_bg <- segregated_id_list_bg[[1]][3][[1]]
      ## }
    } else { # if genes
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

      ## contingencyTb[1, 1] <- ifelse(background_type == "database",
      ##   tot_in_pathway - user_in_pathway,
      ##   bg_in_pathway
      ## )
      ## contingencyTb[1, 2] <- ifelse(background_type == "database",
      ##   tot_out_pathway - user_out_pathway,
      ##   bg_out_pathway
      ## )
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
    allpids <- RMariaDB::dbGetQuery(con, query)
    RMariaDB::dbDisconnect(con)
    pidstorun <- setdiff(allpids[, 1], pid)
    pidstorunlist <- sapply(pidstorun, shQuote)
    pidstorunlist <- paste(pidstorunlist, collapse = ",")

    query2 <- paste0(
      "select rampId,pathwayRampId from analytehaspathway where pathwayRampId in (",
      pidstorunlist, ")"
    )
    con <- connectToRaMP()
    restcids <- RMariaDB::dbGetQuery(con, query2) # [[1]]
    RMariaDB::dbDisconnect(con)

    query1 <- paste0("select rampId,pathwayRampId from analytehaspathway;")

    con <- connectToRaMP()
    allcids <- RMariaDB::dbGetQuery(con, query1) # [[1]]
    RMariaDB::dbDisconnect(con)

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
      which(c(totinpath, totinpath2) >= min_path_size),
      which(c(totinpath, totinpath2) < max_path_size)
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

    # only keep pathways that have >= min_path_size or < max_path_size compounds
    keepers <- intersect(
      which(c(totinpath) >= min_path_size),
      which(c(totinpath) < max_path_size)
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
  return(list(out, pathwaydf))
}

#' Do fisher test for only one pathway from search result
#' clicked on highchart
#' @param analytes a vector of analytes (genes or metabolites) that need to be searched
#' @param NameOrIds whether input is "names" or "ids" (default is "ids", must be the same for analytes and background)
#' @param total_genes number of genes analyzed in the experiment (e.g. background) (default is 20000, with assumption that analyte_type is "genes")
#' @param min_analyte if the number of analytes (gene or metabolite) in a pathway is
#' < min_analyte, do not report
#' @param MCall T/F if true, all pathways are used for multiple comparison corrections; if false, only pathways covering user analytes will be used (default is "F")
#' @param alternative alternative hypothesis test passed on to fisher.test().  Options are two.sided, greater, or less (default is "less")
#' @param min_path_size the minimum number of pathway members (genes and metabolites) to include the pathway in the output (default = 5)
#' @param max_path_size the maximum number of pathway memnbers (genes and metaboltes) to include the pathway in the output (default = 150)
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
#' @return a list containing two entries: [[1]] fishresults, a dataframe containing pathways with Fisher's p values
#' (raw and with FDR and Holm adjustment), number of user analytes in pathway, total number of analytes in pathway,
#' and pathway source ID/database. [[2]] analyte_type, a string specifying the type of analyte input into the function ("genes", "metabolites", or "both")
#' @examples
#' \dontrun{
#' pkg.globals <- setConnectionToRaMP(
#'   dbname = "ramp2", username = "root",
#'   conpass = "", host = "localhost"
#' )
#' analyte.list <- c(
#'   "chebi:15344", "chebi:10983", "chebi:15351",
#'   "uniprot:Q86V21", "uniprot:Q02338", "uniprot:Q9BUT1"
#' )
#'
#' fisher.results <- runCombinedFisherTest(analytes = analyte.list, NameOrIds = "ids")
#' }
#' @export
runCombinedFisherTest <- function(analytes,
                                  NameOrIds = "ids",
                                  total_genes = 20000,
                                  min_analyte = 2,
                                  MCall = F,
                                  alternative = "less",
                                  min_path_size = 5,
                                  max_path_size = 150,
                                  includeRaMPids = FALSE,
                                  background_type = "database",
                                  background = "database") {
  G <- M <- 0

  # Grab pathways that contain metabolites to run Fisher on metabolites
  # This will return all pathways that have at 8-120 metabolites/genes in them
  ## fishmetab <- pathwaydf[grep("RAMP_C_", pathwaydf$rampId), ]

  print("Running Fisher's tests on metabolites")
  outmetab <- runFisherTest(
    analytes = analytes,
    analyte_type = "metabolites",
    total_genes = total_genes,
    MCall = MCall,
    min_path_size = min_path_size,
    max_path_size = max_path_size,
    background_type = background_type,
    background = background
  )
  pathwaydf_metab <- outmetab[[2]]
  outmetab <- outmetab[[1]]
  if (!is.null(outmetab)) {
    M <- 1
  }

  # Grab pathways that contain genes to run Fisher on genes
  ## fishgene <- pathwaydf[grep("RAMP_G_", pathwaydf$rampId), ]
  ## Genes are not evaluated if custom background is specified
  if (background_type == "database") {
    print("Running Fisher's tests on genes")
    outgene <- runFisherTest(
      analytes = analytes,
      analyte_type = "genes",
      total_genes = total_genes,
      MCall = MCall,
      min_path_size = min_path_size,
      max_path_size = max_path_size
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
    ##   NameOrIds = NameOrIds
    ##   )
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
      paste0(collapse = ";")
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
#' @param NameOrIds whether input is "names" or "ids" (default is "ids")
#' @param includeRaMPids include internal RaMP identifiers (default is "FALSE")
#' @return a list contains all metabolites as name and pathway inside.
#' @examples
#' \dontrun{
#' pkg.globals <- setConnectionToRaMP(
#'   dbname = "ramp2", username = "root",
#'   conpass = "", host = "localhost"
#' )
#' mypath <- getPathwayFromAnalyteV2(analytes = c("2-hydroxyglutarate", "glutamate"))
#' }
#' @export
getPathwayFromAnalyte <- function(analytes = "none",
                                  find_synonym = FALSE,
                                  NameOrIds = "ids",
                                  includeRaMPids = FALSE) {
  print("Starting getPathwayFromAnalyte()")
  if (is.null(analytes) || length(analytes) == 0) {
    warning("Input analyte list is NULL or empty. Aborting getPathwayFromAnalyte()")
    return(NULL)
  }

  if (!(NameOrIds %in% c("ids", "names"))) {
    warning(paste0(
      "NameOrIds must have a value in c('ids','names')\n",
      "Supplied NameOrIds falue = ('", NameOrIds, "')\nAborting getPathwayFromAnlyte()"
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

  if (NameOrIds == "ids") {
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
                  and
                  p.type != 'hmdb'
                  group by inputId, rampId, pathwayId, p.pathwayName, p.type, p.pathwayRampId
                  order by pathwayName asc")
  } else {
    print("Working on analyte name list...")
    sql <- paste0(
      "select p.pathwayName, p.type as pathwaySource, p.sourceId as pathwayId, lower(asyn.Synonym) as inputCommonName, group_concat(distinct s.sourceId order by s.sourceId separator '; ') as sourceIds, s.rampId, p.pathwayRampId
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
    and
    p.type != 'hmdb'
    group by inputCommonName, s.rampId, pathwayId, p.pathwayName, p.type, p.pathwayRampId
    order by pathwayName asc
  "
    )
  }

  con <- connectToRaMP()
  df2 <- RMariaDB::dbGetQuery(con, sql)


  if (find_synonym && nrow(df2) > 0) {
    rampIds <- df2[, "rampId"]

    rampIds <- sapply(rampIds, shQuote)
    rampIds <- paste(rampIds, collapse = ",")

    sql <- paste0("select rampId as rampId, group_concat(distinct Synonym order by Synonym separator '; ')
     as synonyms from analytesynonym
     where rampId in (", rampIds, ") group by rampId")

    synonymsDf <- RMariaDB::dbGetQuery(con, sql)

    if (nrow(synonymsDf) > 0) {
      df2 <- merge(df2, synonymsDf, by = "rampId")
    }
  }

  RMariaDB::dbDisconnect(con)

  if (!includeRaMPids && nrow(df2) > 0) {
    df2 <- subset(df2, select = -c(rampId, pathwayRampId))
  }

  print("finished getPathwaytFromAnalyte()")
  print(paste0("Found ", nrow(df2), " associated pathways."))

  return(df2)
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
#'
#' @return list:[[1]] Pathway enrichment result with dataframe having a cluster assignment column added
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
#' clustered.fisher.results <- findCluster(fisher.results)
#' }
#' @export

findCluster <- function(fishers_df, perc_analyte_overlap = 0.5,
                        min_pathway_tocluster = 2, perc_pathway_overlap = 0.5) {
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
  con <- connectToRaMP()
  idkey <- RMariaDB::dbGetQuery(con, query) %>%
    dplyr::rename("pathwayId" = "sourceId") ##  %>%
  ## dplyr::rename("rampId" = "pathwayRampId")

  rampToSource <- function(x) {
    out <- with(idkey, {
      idkey %>%
        dplyr::filter(pathwayRampId == x) %>%
        dplyr::pull("pathwayId")
    })
    return(out)
  }

  RMariaDB::dbDisconnect(con)

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
    print("Finished clustering pathways...")
    return(output)
  }
}
