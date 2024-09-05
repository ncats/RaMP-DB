#' Function that query database to find ontology information based on
#' the given list of analytes
#' @param analytes a vector of analytes or a analytes delimited by new line character
#' @param namesOrIds specify the type of given data
#' @param includeRaMPids whether or not to include RaMP ids in the output (TRUE/FALSE)
#' @param db a RaMP database object
#' @return dataframe that contains searched ontology from given analytes
#'
#' @examples
#' \dontrun{
#' pkg.globals <- setConnectionToRaMP(
#'   dbname = "ramp2", username = "root",
#'   conpass = "", host = "localhost"
#' )
#' getOntoFromMeta("hmdb:HMDB0071437")
#' }
#' @export
getOntoFromMeta <- function(analytes, namesOrIds = "ids", includeRaMPids = FALSE, db = RaMP()) {
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

  allontos <- getOntologies(db = db)
  matched_ontos <- unlist(lapply(
    list_ontology,
    function(x) grep(paste0("^", x, "$"), allontos$commonName)
  ))

  # Only proceed if df has anything returned
  if (length(matched_ontos) > 0) {
    print(paste0("Found ", length(matched_ontos), " ontology term matches."))

    ontologyList <- paste0("'", paste(list_ontology, collapse = "','"), "'")

    sql <- paste0("select rampId,
          group_concat(distinct s.sourceId order by s.sourceId separator '; ') as source_ids,
          group_concat(distinct s.commonName order by s.commonName separator '; ') as common_names, o.commonName, o.HMDBOntologyType
          from source s, analytehasontology ao, ontology o where ao.rampOntologyId in (
          select distinct rampOntologyId from ontology where commonName in (", ontologyList, "))
          and o.rampOntologyId = ao.rampOntologyId and s.rampId = ao.rampCompoundId
          group by o.commonName, s.rampId, o.HMDBOntologyType")

    if (.is_sqlite(db)) {
      sql <- paste0("select rampId,
          group_concat(distinct s.sourceId COLLATE NOCASE) as source_ids,
          group_concat(distinct s.commonName COLLATE NOCASE) as common_names, o.commonName, o.HMDBOntologyType
          from source s, analytehasontology ao, ontology o where ao.rampOntologyId in (
          select distinct rampOntologyId from ontology where commonName in (", ontologyList, "))
          and o.rampOntologyId = ao.rampOntologyId and s.rampId = ao.rampCompoundId
          group by o.commonName, s.rampId, o.HMDBOntologyType")
    }

    mdf_final <- RaMP::runQuery(sql, db)

    mdf_final <- unique(mdf_final)
    mdf_final <- mdf_final[, c(4, 5, 3, 2)]
    colnames(mdf_final) <- c("ontologyTerm", "ontologyCategory", "metNames", "metIds")

    print(paste0("Found ", nrow(mdf_final), " metabolites associated with the input ontology terms."))
    print("Finished getting metabolies from ontology terms.")

    return(mdf_final)
  } else {
    warning("The input ontology terms were not found in RaMP.\nRun the getOntologies() function to see available ontology terms.")
    return(NA)
  }
}


#' Enrichment analysis for metabolite ontology mappings
#' @importFrom rlang .data
#' @param analytes a vector of analytes (genes or metabolites) that need to be searched
#' @param namesOrIds whether input is "names" or "ids" (default is "ids", must be the same for analytes and background)
#' @param alternative alternative hypothesis test passed on to fisher.test().  Options are two.sided, greater, or less (default is "less")
#' @param min_analyte minimum number of analytes per pathway (pathways with < min_analyte analytes are filtered out).
#' @param min_ontology_size the minimum number of ontology members (genes and metabolites) to include the ontology in the output (default = 5)
#' @param max_ontology_size the maximum number of ontology memnbers (genes and metaboltes) to include the ontology in the output (default = 150)
#' @param includeRaMPids whether or not to include RaMP IDs in the output (TRUE/FALSE)
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
#' @param db a RaMP database object
#' @return a dataframe with columns containing pathway ID, fisher's p value, user analytes in pathway, and total analytes in pathway
#' @export
#' @importFrom methods is

runOntologyTest <- function(analytes,
                            namesOrIds = "ids",
                            alternative = "less", min_analyte = 2,
                            min_ontology_size = 5, max_ontology_size = 1500,
                            includeRaMPids = FALSE,
                            background_type = "database", background = "database",
                            db = RaMP()) {
  now <- proc.time()
  print("Fisher Testing ......")

  ontologydf <- getOntoFromMeta(
    db = db, analytes = analytes,
    includeRaMPids = TRUE,
    namesOrIds = namesOrIds
  )
  ontologyRampId <- rampId <- c()


  # moved this check until we determine if we have analytes of a given type.
  if (nrow(ontologydf) == 0) {
    return(NULL)
  }

  #  if(class(background_type)=="list"){
  if (is(background_type, "list")) {
    background <- unlist(background)
  }

  if (background_type == "list") {
    backgrounddf <- getOntoFromMeta(
      db = db, background,
      includeRaMPids = TRUE,
      namesOrIds = namesOrIds
    )
  } else if (background_type == "file") {
    userbkg <- utils::read.table(background, header = F)[, 1]
    backgrounddf <- getOntoFromMeta(
      db = db, analytes = userbkg,
      includeRaMPids = TRUE,
      namesOrIds = namesOrIds
    )
  } else if (background_type == "biospecimen") {
    biospecimen <- background
    if (biospecimen == "Adipose") {
      biospecimen <- "Adipose tissue"
    }


    # Get metabolites that belong to a specific biospecimen
    # query <- paste0(
    #   "SELECT analytehasontology.*, ontology.*, analytehasontology.* from analytehasontology, ontology, analytehasontology where ontology.commonName in ('",
    #   biospecimen,
    #   "') and ontology.rampOntologyId = analytehasontology.rampOntologyId and analytehasontology.rampCompoundId = analytehasontology.rampId"
    # )

    # less data pull-back
    query <- paste0(
      "SELECT analytehasontology.* from analytehasontology, ontology, analytehasontology where ontology.commonName in ('",
      biospecimen,
      "') and analytehasontology.rampOntologyId = ontology.rampOntologyId and analytehasontology.rampCompoundId = analytehasontology.rampId"
    )

    backgrounddf <- RaMP::runQuery(query, db)

    if (nrow(backgrounddf) == 0) {
      stop("Biospecimen background not found. Choices are 'Blood', 'Adipose', 'Heart', 'Urine', 'Brain', 'Liver', 'Kidney', 'Saliva', and 'Feces'")
    }

    # only keep the input metabolites (converted into ontologydf in line above) that are in the biospecimen type specified
    ontologydf <- with(ontologydf, {
      ontologydf %>%
        dplyr::filter(.data$rampId %in% .data$backgrounddf$rampId)
    })
    if (nrow(ontologydf) == 0) {
      stop("There are no metabolites in your input that map to your selected biospecimen")
    }
  } else if (background_type == "database") {
    # do nothing, it's handled down below in if statements
  } else {
    stop("background_type was not specified correctly.  Please specify one of the following options: database, file, list, biospecimen")
  }


  ## Check that all metabolites of interest are in the background
  if (background_type != "database") {
    if (length(setdiff(ontologydf$rampId, backgrounddf$rampId) != 0)) {
      stop("All analytes in set of interest must also be in background")
    }
  }

  ## Initialize empty contingency table for later
  contingencyTb <- matrix(0, nrow = 2, ncol = 2)
  colnames(contingencyTb) <- c("In Ontology", "Not In Ontology")
  rownames(contingencyTb) <- c("All Metabolites", "User's Metabolites")
  ## Get ontology ids that contain the user analytes
  pid <- unique(ontologydf$rampOntologyId)
  list_pid <- sapply(pid, shQuote)
  list_pid <- paste(list_pid, collapse = ",")

  # Get the total number of metabolites that are mapped to ontologys in RaMP (that's the default background)
  # added conditional to not pull hmdb ids
  query <- "select distinct rampCompoundId from analytehasontology;"

  allids <- RaMP::runQuery(query, db)
  allids <- allids[!duplicated(allids), ]

  totanalytes <- length(allids)

  ## Input_RampIds is a table of all analytes included in ontologys represented in the user set
  ## "User" refers to significant analytes

  ## Get pathway ids that contain the user analytes
  pid <- unique(ontologydf$rampOntologyId)
  list_pid <- sapply(pid, shQuote)
  list_pid <- paste(list_pid, collapse = ",")

  ## Retrieve compound ids associated with background pathways and count
  query <- paste0(
    "select * from analytehasontology where rampOntologyId in (",
    list_pid, ")"
  )

  input_RampIds <- RaMP::runQuery(query, db)

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
      if (background_type != "database") {
        ids_inontology_bg <- backgrounddf[which(backgrounddf$rampOntologyId == i), "rampId"]
        bg_in_ontology <- length(unique(grep("RAMP_C", ids_inontology_bg, value = TRUE)))
      }
    }

    if (pidCount == 1) {
      tot_user_analytes <- length(grep("RAMP_C", unique(ontologydf$rampCompoundId)))
      if (background_type != "database") {
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
      if (background_type != "database") {
        bg_in_ontology <- length(unique(backgrounddf[which(backgrounddf$rampOntologyId == i), "rampCompoundId"]))
      }

      user_out_ontology <- tot_user_analytes - user_in_ontology

      if (background_type != "database") {
        bg_in_ontology <- length(unique(backgrounddf[which(backgrounddf$rampOntologyId == i), "rampCompoundId"]))
        bg_out_ontology <- tot_bg_analytes - bg_in_ontology
      }

      if (background_type == "database") {
        contingencyTb[1, 1] <- tot_in_ontology - user_in_ontology
      } else {
        contingencyTb[1, 1] <- bg_in_ontology
      }
      if (background_type == "database") {
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

  # only keep ontologys that have >= min_ontology_size or < max_ontology_size compounds
  keepers <- intersect(
    which(c(totinontology) >= min_ontology_size),
    which(c(totinontology) < max_ontology_size)
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
    Odds_Ratio = oddsratio[keepers],
    Num_In_Ontology = userinontology[keepers],
    Total_In_Ontology = totinontology[keepers]
  )

  fdr <- stats::p.adjust(out$Pval, method = "fdr")
  out <- cbind(out, fdr)
  colnames(out)[ncol(out)] <- "Pval_FDR"
  holm <- stats::p.adjust(out$Pval, method = "holm")
  out <- cbind(out, holm)
  colnames(out)[ncol(out)] <- "Pval_Holm"
  keepers <- which(out$Num_In_Ontology >= min_analyte)
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
      paste0(collapse = ";")
    return(analytes)
  })

  # for user is the output needed, based on what user input
  if (includeRaMPids) {
    return(out)
  } else {
    return(out %>% cleanup())
  }
}
