#' Function that query database to find ontology information based on
#' the given list of metabolites
#' @param mets a vector of metabolites or a metabolites delimited by new line character
#' @param namesOrIds specify the type of given data
#' @param includeRaMPids whether or not to include RaMP ids in the output (TRUE/FALSE)
#' @param minOntologySize the minimum number of metabolites in an ontology for it to be included in results. Default is 1,000
#' @param maxOntologySize the maximum number of metabolites in an ontology for it to be included in results. Default is Inf
#' @param db a RaMP database object
#' @return dataframe that contains searched ontology from given metabolites
#'
#' @examples
#' \dontrun{
#' rampDB <- RaMP()
#' getOntoFromMeta(mets = "hmdb:HMDB0071437", db=rampDB)
#' }
#' @export
getOntoFromMeta <- function(mets, namesOrIds = "ids", includeRaMPids = FALSE, minOntologySize = 1E3, maxOntologySize = Inf, db = RaMP()) {
  if (!(namesOrIds %in% c("ids", "names"))) {
    stop("Specifiy the type of given data to 'ids' or 'names'")
  }

  now <- proc.time()
  if (is.character(mets)) {
    if (grepl("\n", mets)[1]) {
      list_metabolite <- strsplit(mets, "\n")
      list_metabolite <- unlist(list_metabolite)
    } else if (grepl(",", mets)[1]) {
      list_metabolite <- strsplit(mets, "\n")
      list_metabolite <- unlist(list_metabolite)
    } else {
      list_metabolite <- mets
    }
  } else if (is.data.frame(mets)) {
    list_metabolite <- unlist(mets)
  }
  list_metabolite <- unique(list_metabolite)

  if (namesOrIds == "ids") {
    df <- db@api$getSourceDataForAnalyteIDs(analyteIDs = list_metabolite)
  } else if (namesOrIds == "names") {
    df <- db@api$getSourceDataForAnalyteNames(analyteNames = list_metabolite)
  }

  if (nrow(df) == 0) {

    if (namesOrIds == "ids") {
      message("These ids do not exist in the source table")
    } else {
      message("These names do not exist in the metabolite synonym table")
    }
    return(NULL)
  }

  rampid <- unique(df$rampId)

  temp_ontologies_df <- getOntologies(db = db) %>%
    dplyr::filter(.data$metCount <= maxOntologySize,
                  .data$metCount > minOntologySize)

  df2 <- db@api$getOntologiesForRampIDs(rampIds = rampid) %>%
    dplyr::filter(.data$rampOntologyId %in% temp_ontologies_df$rampOntologyId)

  if (nrow(df2) == 0) {
    message("No searching result because these metabolites are not linked to ontology in search parameters")
    return(NULL)
  }

  rampontoid <- unique(df2$rampOntologyId)
  df3 <- db@api$getOntologyData(rampIds = rampontoid)

  mdf <- unique(merge(df3, df2, all.x = T))
  mdf <- unique(merge(mdf, df,
    all.x = T, by.x = "rampCompoundId",
    by.y = "rampId"
  ))
  colnames(mdf)[colnames(mdf) == "commonName.x"] <- "Ontology"
  colnames(mdf)[colnames(mdf) == "commonName.y"] <- "Metabolites"

  # KJK - deduplicate rows ignoring case for Metabolites
  filtered_column_names = c("rampCompoundId", "rampOntologyId", "sourceId", "IDtype", "Ontology", "HMDBOntologyType")
  if (!includeRaMPids) {
    filtered_column_names = c("sourceId", "IDtype", "Ontology", "HMDBOntologyType")
  }

  final_column_names = c("Metabolites", filtered_column_names)
  mdf <- mdf[ , final_column_names, drop = FALSE]

  mdf$insensitive_metabolites = tolower(mdf$Metabolites)
  dedup_column_names = c("insensitive_metabolites", filtered_column_names)
  mdf <- mdf[!duplicated(mdf[ , dedup_column_names]), ]
  mdf$insensitive_metabolites <- NULL

  return(mdf)
}


#' function that query database to find mets in given ontologies
#' @param ontology a vector of ontology or ontologies delimited by new line character
#' @param minOntologySize the minimum number of metabolites in an ontology for it to be included in results. Default is 1,000
#' @param maxOntologySize the maximum number of metabolites in an ontology for it to be included in results. Default is Inf
#' @param curate filter searchable ontologies to only those visible to filter by when searching the HMDB website. Primarily for front-end UI use.
#' @param db a RaMP database object
#' @return dataframe that contains searched mets from given ontology
#' @examples
#' \dontrun{
#' ontologies.of.interest <- c("Colon", "Liver", "Lung")
#'
#' new.metabolites <- RaMP::getMetaFromOnto(db = rampDB, ontology = ontologies.of.interest)
#' }
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @export
getMetaFromOnto <- function(ontology, minOntologySize = 1E3, maxOntologySize = Inf, curate = F, db = RaMP()) {
  print("Retreiving Metabolites for input ontology terms.")
  now <- proc.time()
  if (is.character(ontology)) {
    if (grepl("\n", ontology)[1]) {
      list_ontology <- strsplit(ontology, "\n")
      list_ontology <- unlist(list_ontology)
    } else if (grepl(",", ontology)[1]) {
      list_ontology <- strsplit(ontology, ",")
      list_ontology <- unlist(list_ontology)
    } else {
      list_ontology <- ontology
    }
  } else if (is.data.frame(ontology)) {
    list_ontology <- unlist(ontology)
  }

  list_ontology <- unique(list_ontology)

  allontos <- getOntologies(db = db) %>%
    dplyr::filter(.data$metCount <= maxOntologySize,
                  .data$metCount > minOntologySize)

  #This if statement is meant for use with the front end.
  #The idea is to limit the ontologies searchable in the RaMP UI to be the same ones HMDB lets the user filter by on their website.
  #This is hard-coded as a filter here and will need to be updated manually in the current form.
  #The default 'curate' argument should be false.
  #Adam 10/3/2024
  if (curate) {
    allontos <- allontos %>% dplyr::filter(.data$commonName %in% c('Blood',
                                                'Urine',
                                                'Saliva',
                                                'Cerebrospinal fluid',
                                                'Feces',
                                                'Sweat',
                                                'Breast milk',
                                                'Bile',
                                                'Amniotic fluid',
                                                'Exogenous',
                                                'Endogenous',
                                                'Food',
                                                'Plant',
                                                'Microbe',
                                                'Cosmetic',
                                                'Drug',
                                                'Cell membrane',
                                                'Cytoplasm',
                                                'Nucleus',
                                                'Mitochondria'))
  }

  matched_ontos <- list_ontology[list_ontology %in% allontos$commonName]

  # Only proceed if df has anything returned
  if (length(matched_ontos) > 0) {
    print(paste0("Found ", length(matched_ontos), " ontology term matches."))

    mdf_final <- db@api$getMetabolitesForOntology(ontologyList = matched_ontos)

    mdf_final <- unique(mdf_final)
    mdf_final <- mdf_final[, c(4, 5, 3, 2)]
    colnames(mdf_final) <- c("ontologyTerm", "ontologyCategory", "metNames", "metIds")

    print(paste0("Found ", nrow(mdf_final), " metabolites associated with the input ontology terms."))
    print("Finished getting metabolies from ontology terms.")

    return(mdf_final)
  } else {
    warning("The input ontology terms were not found in RaMP within input parameters.\nRun the getOntologies() function to see available ontology terms.")
    return(NA)
  }
}


#' Enrichment analysis for metabolite ontology mappings
#' @importFrom rlang .data
#' @param mets a vector of metabolites that need to be searched
#' @param namesOrIds whether input is "names" or "ids" (default is "ids", must be the same for mets and background)
#' @param alternative alternative hypothesis test passed on to fisher.test().  Options are two.sided, greater, or less (default is "less")
#' @param minMets minimum number of mets per pathway (pathways with < minMets mets are filtered out).
#' @param minOntologySize the minimum number of ontology members (genes and metabolites) to include the ontology in the output (default = 5)
#' @param maxOntologySize the maximum number of ontology members (genes and metabolites) to include the ontology in the output (default = 150)
#' @param includeRaMPids whether or not to include RaMP IDs in the output (TRUE/FALSE)
#' @param backgroundType type of background that is input by the user.  Options are "database" if user wants all
#' mets from the RaMP database will be used; "file", if user wants to input a file with a list of background
#' mets; "list", if user wants to input a vector of analyte IDs; "biospecimen", if user wants to specify a
#' biospecimen type (e.g. blood, adipose tissue, etc.) and have those biospecimen-specific mets used.  For genes,
#' only the "database" option is used.
#' @param background background to be used for Fisher's tests.  If parameter 'backgroundType="database"', this parameter
#' is ignored (default="database"); if parameter 'backgroundType= "file"', then 'background' should be a file name (with
#' directory); if 'backgroundType="list"', then 'background' should be a vector of RaMP IDs; if 'backgroundType="biospecimen"'
#' then users should specify one of the following: "Blood", "Adipose tissue", "Heart", "Urine", "Brain", "Liver", "Kidney",
#' "Saliva", and "Feces"
#' @param db a RaMP database object
#' @return a dataframe with columns containing pathway ID, fisher's p value, user mets in pathway, and total mets in pathway
#' @examples
#' \dontrun{
#' ontologies.enriched <- runEnrichOntologies(mets = c("hmdb:HMDB0000033","hmdb:HMDB0000052","hmdb:HMDB0000094", "hmdb:HMDB0000161","hmdb:HMDB0000168","hmdb:HMDB0000191","hmdb:HMDB0000201","chemspider:10026", "hmdb:HMDB0006059", "Chemspider:6405", "CAS:5657-19-2","hmdb:HMDB0002511", "chemspider:20171375", "CAS:133-32-4","CAS:5746-90-7", "CAS:477251-67-5", "hmdb:HMDB0000695", "chebi:15934", "CAS:838-07-3", "hmdb:HMDBP00789", "hmdb:HMDBP00283", "hmdb:HMDBP00284", "hmdb:HMDBP00850"))
#' }
#' @export
#' @importFrom methods is

runEnrichOntologies <- function(mets,
                            namesOrIds = "ids",
                            alternative = "less", minMets = 2,
                            minOntologySize = 5, maxOntologySize = 1500,
                            includeRaMPids = FALSE,
                            backgroundType = "database", background = "database",
                            db = RaMP()) {
  now <- proc.time()
  print("Fisher Testing ......")

  ontologydf <- getOntoFromMeta(
    db = db, mets = mets,
    includeRaMPids = TRUE,
    namesOrIds = namesOrIds
  )
  ontologyRampId <- rampId <- c()


  # moved this check until we determine if we have mets of a given type.
  if (nrow(ontologydf) == 0) {
    return(NULL)
  }

  #  if(class(backgroundType)=="list"){
  if (is(backgroundType, "list")) {
    background <- unlist(background)
  }

  if (backgroundType == "list") {
    backgrounddf <- getOntoFromMeta(
      db = db, mets = background,
      includeRaMPids = TRUE,
      namesOrIds = namesOrIds
    )
  } else if (backgroundType == "file") {
    userbkg <- utils::read.table(background, header = F)[, 1]
    backgrounddf <- getOntoFromMeta(
      db = db, mets = userbkg,
      includeRaMPids = TRUE,
      namesOrIds = namesOrIds
    )
  } else if (backgroundType == "biospecimen") {
    biospecimen <- background
    if (biospecimen == "Adipose") {
      biospecimen <- "Adipose tissue"
    }

    backgrounddf <- db@api$getAnalytesFromOntology(biospecimen = biospecimen)

    if (nrow(backgrounddf) == 0) {
      stop("Biospecimen background not found. Choices are 'Blood', 'Adipose', 'Heart', 'Urine', 'Brain', 'Liver', 'Kidney', 'Saliva', and 'Feces'")
    }

    # only keep the input metabolites (converted into ontologydf in line above) that are in the biospecimen type specified
    ontologydf <- with(ontologydf, {
      ontologydf %>%
        dplyr::filter(.data$rampCompoundId %in% backgrounddf$rampCompoundId)
    })
    if (nrow(ontologydf) == 0) {
      stop("There are no metabolites in your input that map to your selected biospecimen")
    }
  } else if (backgroundType == "database") {
    # do nothing, it's handled down below in if statements
  } else {
    stop("backgroundType was not specified correctly.  Please specify one of the following options: database, file, list, biospecimen")
  }


  ## Check that all metabolites of interest are in the background
  if (backgroundType != "database") {
    if (length(setdiff(ontologydf$rampId, backgrounddf$rampId) != 0)) {
      stop("All mets in set of interest must also be in background")
    }
  }

  ## Initialize empty contingency table for later
  contingencyTb <- matrix(0, nrow = 2, ncol = 2)
  colnames(contingencyTb) <- c("In Ontology", "Not In Ontology")
  rownames(contingencyTb) <- c("All Metabolites", "User's Metabolites")

  # Get the total number of metabolites that are mapped to ontologys in RaMP (that's the default background)
  totanalytes <- db@api$getMetaboliteWithOntologyCount()

  ## Input_RampIds is a table of all mets included in ontologys represented in the user set
  ## "User" refers to significant mets

  ## Get pathway ids that contain the user mets
  pid <- unique(ontologydf$rampOntologyId)

  ## Retrieve compound ids associated with background pathways and count
  input_RampIds <- db@api$getRampIDsForOntologies(ontologyIDs = pid)

  if (is.null(input_RampIds)) {
    stop("Data doesn't exist")
  } else {
    # data frames for metabolites with pathwayRampID, Freq based  on Source(kegg, reactome, wiki)
    unique_input_RampId_C <- unique(input_RampIds[, c(
      "rampCompoundId",
      "rampOntologyId"
    )])

    freq_unique_input_RampId_C <- as.data.frame(table(unique_input_RampId_C[, "rampOntologyId"]))

    names(freq_unique_input_RampId_C)[1] <- "rampOntologyId"

    input_metab <- freq_unique_input_RampId_C
    names(input_metab) <- c("ontology", "Freq")
  }

  # Loop through each ontology, build the contingency table, and calculate Fisher's Exact
  # test p-value
  pidCount <- 0
  pval <- oddsratio <- totinontology <- userinontology <- pidused <- c()
  for (i in pid) {
    ids_inontology <- ontologydf[which(ontologydf$rampOntologyId == i), "rampCompoundId"]

    pidCount <- pidCount + 1
    # Check to make sure that this ontology does have metabolites
    if (length(grep("RAMP_C", ids_inontology)) == 0) {
      user_in_ontology <- 0
    } else {
      user_in_ontology <- length(unique(grep("RAMP_C", ids_inontology, value = TRUE)))
      if (backgroundType != "database") {
        ids_inontology_bg <- backgrounddf[which(backgrounddf$rampOntologyId == i), "rampId"]
        bg_in_ontology <- length(unique(grep("RAMP_C", ids_inontology_bg, value = TRUE)))
      }
    }

    if (pidCount == 1) {
      tot_user_analytes <- length(grep("RAMP_C", unique(ontologydf$rampCompoundId)))
      if (backgroundType != "database") {
        tot_bg_analytes <- length(grep("RAMP_C", unique(backgrounddf$rampCompoundId)))
      }
    }
    total_ontology_analytes <- totanalytes
    tot_in_ontology <- input_metab %>%
      dplyr::filter(.data$ontology == i) %>%
      dplyr::pull(.data$Freq)
    if (tot_in_ontology == 0 || user_in_ontology == 0) {
      pval <- c(pval, NA)
      oddsratio <- c(oddsratio, NA)
    } else {
      tot_out_ontology <- total_ontology_analytes - tot_in_ontology
      # fill the rest of the table out

      ## user_in_ontology <- length(unique(ontologydf[which(ontologydf$rampOntologyId==i),"rampId"]))
      if (backgroundType != "database") {
        bg_in_ontology <- length(unique(backgrounddf[which(backgrounddf$rampOntologyId == i), "rampCompoundId"]))
      }

      user_out_ontology <- tot_user_analytes - user_in_ontology

      if (backgroundType != "database") {
        bg_in_ontology <- length(unique(backgrounddf[which(backgrounddf$rampOntologyId == i), "rampCompoundId"]))
        bg_out_ontology <- tot_bg_analytes - bg_in_ontology
      }

      if (backgroundType == "database") {
        contingencyTb[1, 1] <- tot_in_ontology - user_in_ontology
      } else {
        contingencyTb[1, 1] <- bg_in_ontology
      }
      if (backgroundType == "database") {
        contingencyTb[1, 2] <- tot_out_ontology - user_out_ontology
      } else {
        contingencyTb[1, 2] <- bg_out_ontology
      }
      contingencyTb[2, 1] <- user_in_ontology
      contingencyTb[2, 2] <- user_out_ontology
      # Put the test into a try catch in case there's an issue, we'll have some details on the contingency matrix
      tryCatch(
        {
          result <- stats::fisher.test(contingencyTb, alternative = alternative)
        },
        error = function(e) {
          print(toString(e))
          print(i)
          print(contingencyTb)
          print(tot_in_ontology)
          print(tot_out_ontology)
          print(user_in_ontology)
          print(user_out_ontology)
          print(ontologydf)
        }
      )
      pval <- c(pval, result$p.value)
      oddsratio <- c(oddsratio, result$estimate)
    } # End else tot_in_ontology is not zero

    userinontology <- c(userinontology, user_in_ontology)
    totinontology <- c(totinontology, tot_in_ontology)
    pidused <- c(pidused, i)
  } # end for loop

  print("")
  print(now - proc.time())
  print("")

  # only keep ontologys that have >= minOntologySize or < maxOntologySize compounds
  keepers <- intersect(
    which(c(totinontology) >= minOntologySize),
    which(c(totinontology) < maxOntologySize)
  )

  # hist(totinontology,breaks=1000)
  print(paste0("Keeping ", length(keepers), " ontologys"))
  # fdr <- stats::p.adjust(c(pval,pval2)[keepers],method="fdr")
  # holm <- stats::p.adjust(c(pval,pval2)[keepers],method="holm")
  print(paste0("Calculated p-values for ", length(pval), " ontologys"))

  # format output (retrieve ontology name for each unique source id first
  out <- data.frame(
    rampOntologyId = pidused[keepers],
    Pval = pval[keepers], # FDR.Adjusted.Pval=fdr,
    # Holm.Adjusted.Pval=holm,
    OR = oddsratio[keepers],
    Num_In_Ontology = userinontology[keepers],
    Total_In_Ontology = totinontology[keepers]
  )

  fdr <- stats::p.adjust(out$Pval, method = "fdr")
  out <- cbind(out, fdr)
  colnames(out)[ncol(out)] <- "Pval_FDR"
  holm <- stats::p.adjust(out$Pval, method = "holm")
  out <- cbind(out, holm)
  colnames(out)[ncol(out)] <- "Pval_Holm"
  keepers <- which(out$Num_In_Ontology >= minMets)
  out <- merge(
    ontologydf[, c(
      "Ontology", "rampOntologyId",
      "HMDBOntologyType"
    )],
    out[keepers, ],
    by = "rampOntologyId"
  )

  out <- out[!duplicated(out), ]
  out$analytes <- apply(out, 1, function(x) {
    ontologyid <- x["rampOntologyId"]
    sigontologydf <- ontologydf[which(ontologydf$rampOntologyId == ontologyid), ]
    analytes <- sigontologydf[, "Metabolites"] %>%
      paste0(collapse = "|")
    return(analytes)
  })

  # for user is the output needed, based on what user input
  if (includeRaMPids) {
    return(list(fishertresults = out, result_type = "ontology_enrichment"))
  } else {
    return(list(fishertresults = cleanup(data = out), result_type = "ontology_enrichment" ))
  }
}
