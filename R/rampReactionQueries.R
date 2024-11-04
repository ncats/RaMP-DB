# RaMP Reaction Queries

#' getReactionsForAnalytes returns all associated reactions for a collection of input compound ids
#'
#' @param analytes list of analytes. chebi and/or uniprot ids are required.
#' @param onlyHumanMets boolean to only return pathways containing only human metabolites (ChEBI ontology) (dev
#' in progress)
#' @param humanProtein boolean to only control pathways catalyzed by a human proteins (having human Uniprot)
#' (dev in progress)
#' @param includeTransportRxns if TRUE, returns metabolic and transport reactions, default is TRUE
#' @param rxnDirs character vector of length > 1, specifying reaction directions to return  c("UN", "LR", "RL",
#' "BD", "ALL"), default = c("UN").
#' @param includeRxnURLs if TRUE, urls to Rhea.org will be delivered in the result dataframe for each reaction
#' @param db a RaMP database object
#'
#' @return a list of 3 dataframes [[1]] met2rxn, reaction information for all metabolite inputs (chebi),  [[2]] prot2rxn, reaction information for all protein inputs (uniprot), and [[3]] metProteinCommonReactions, common reactions found across both metabolite and protein inputs (overlap between met2rxn and prot2rxn)
#' @examples
#' \dontrun{
#' analytes.of.interest <- c('chebi:58115', 'chebi:456215', 'chebi:58245', 'chebi:58450', 'chebi:17596',
#' 'chebi:16335', chebi:16750', 'chebi:172878', 'chebi:62286', 'chebi:77897', 'uniprot:P30566','uniprot:P30520',
#' 'uniprot:P00568', 'uniprot:P23109', 'uniprot:P22102', 'uniprot:P15531')
#'
#' reactionsLists <- getReactionsForAnalytes(analytes = analytes.of.interest,
#'     includeTransportRxns = T, humanProtein = T, db = rampDB )
#' }
#' @export
#'
getReactionsForAnalytes <- function( analytes, onlyHumanMets=F, humanProtein=T, includeTransportRxns=T,
	rxnDirs=c("UN"), includeRxnURLs=F, db = RaMP() ) {

  message("Running getReactionsForAnalytes()")

  genes <- data.frame()
  metabolites <- data.frame()
  mdf <- data.frame()
  gdf <- data.frame()

  idLists = checkReactionInputIds("getReactionsForAnalytes", analytes)

  if((length(idLists$proteinIds) + length(idLists$metIds) == 0)) {
    message("Wrong ID types. chebi and/or uniprot ids are required. Returning empty dataframe.")
    return(data.frame())
  }

  mets <- idLists$metIds
  proteins <- idLists$proteinIds

  if(!is.null(mets) && length(mets)>0) {
    mdf <- db@api$getRheaRxnPartnersFromMetIDs(metaboliteIDs=mets, onlyHumanMets=onlyHumanMets, humanProtein=humanProtein, includeTransportRxns=includeTransportRxns, rxnDirs=rxnDirs)
  }

  if(!is.null(proteins) && length(proteins)>0) {
    gdf <- db@api$getRheaRxnPartnersFromGeneIDs(geneIDs=proteins, onlyHumanMets=onlyHumanMets, humanProtein=humanProtein, includeTransportRxns=includeTransportRxns, rxnDirs=rxnDirs)
  }

  resultList <- list()
  resultList[['met2rxn']] <- mdf
  resultList[['prot2rxn']] <- gdf
  resultList[['metProteinCommonReactions']] <- data.frame()

  # find common reactions
  if(nrow(mdf) > 0 && nrow(gdf) > 0) {
    mRampRxnIds <- unlist(mdf$rxn_source_id)
    gRampRxnIds <- unlist(gdf$rxn_source_id)
    commonRxnIds <- intersect(mRampRxnIds, gRampRxnIds)
    if(length(commonRxnIds) > 0) {
      mRxn <- mdf[mdf$rxn_source_id %in% commonRxnIds,]
      gRxn <- gdf[gdf$rxn_source_id %in% commonRxnIds,]

      # one row per reaction, include concat mets, concat proteins
      # probably a more streamlined way to do this...
      commonRampRxnList <- c()
      commonRxnList <- c()
      mSourceIds <- c()
      mRxnSourceNames <- c()
      gRxnSourceIds <- c()
      gRxnSourceName <- c()
      gRxnUniprot <- c()
      commonRxnLabel <- c()
      commonRxnHTMLEq <- c()
      hasHumanProtein <- c()
      onlyHumanMets <- c()
      isTransport <- c()
      rxnDir <- c()
      for(rxn in commonRxnIds) {
        commonRampRxnList <- c(commonRampRxnList, rxn)
        commonRxnList <- c(commonRxnList, paste0(unique(mRxn$rxn_source_id[mRxn$rxn_source_id == rxn]), collapse = ','))
        mSourceIds <- c(mSourceIds, paste0(unique(mRxn$met_source_id[mRxn$rxn_source_id == rxn]), collapse = ','))
        mRxnSourceNames <- c(mRxnSourceNames, paste0(unique(mRxn$common_name[mRxn$rxn_source_id == rxn]), collapse = ','))

        #gRxnSourceIds <- c(gRxnSourceIds, paste0(unique(gRxn$sourceId[gRxn$rxn_source_id == rxn]), collapse = ','))
        gRxnUniprot <- c(gRxnUniprot, paste0(unique(gRxn$uniprot[gRxn$rxn_source_id == rxn]), collapse = ','))
        gRxnSourceName <- c(gRxnSourceName, paste0(unique(gRxn$common_name[gRxn$rxn_source_id == rxn]), collapse = ','))

        rxnDir <- c(rxnDir, paste0(unique(mRxn$direction[mRxn$rxn_source_id == rxn]), collapse = ','))

        commonRxnLabel <- c(commonRxnLabel, paste0(unique(mRxn$label[mRxn$rxn_source_id == rxn]), collapse = ','))
        commonRxnHTMLEq <- c(commonRxnHTMLEq, paste0(unique(mRxn$html_equation[mRxn$rxn_source_id == rxn]), collapse = ','))

        isTransport <- as.integer(c(isTransport, paste0(unique(mRxn$is_transport[mRxn$rxn_source_id == rxn]), collapse = ',')))
        hasHumanProtein <- as.integer(c(hasHumanProtein, paste0(unique(mRxn$has_human_prot[mRxn$rxn_source_id == rxn]), collapse = ',')))
        justHumanMets <- as.integer(c(onlyHumanMets, paste0(unique(mRxn$only_human_mets[mRxn$rxn_source_id == rxn]), collapse = ',')))
      }
      commonReactions <- data.frame(list('metabolites' = mSourceIds, 'met_names' = mRxnSourceNames,
                                         'uniprot' = gRxnUniprot, 'proteinNames' = gRxnSourceName,
                                         'reactionId' = commonRxnList, 'rxn_dir' = rxnDir,'rxn_label' = commonRxnLabel,
                                         'rxn_html_label' = commonRxnHTMLEq,
                                         'is_transport' = isTransport,
                                         'has_human_protein' = hasHumanProtein,
                                         'only_human_mets' = justHumanMets))



      resultList[['metProteinCommonReactions']] <- commonReactions
    }
  }

  # let's add analyte hits per pathway to help sort by rxns with greater support
  if(nrow(resultList$met2rxn) > 0) {
    metCountDf <- data.frame(table(resultList$met2rxn$rxn_source_id))
    colnames(metCountDf) <- c("rxn_source_id", "metHitCount")
    resultList$met2rxn <- merge(x=resultList$met2rxn, y=metCountDf,by.x="rxn_source_id", by.y="rxn_source_id")
    resultList$met2rxn <- resultList$met2rxn[order(-resultList$met2rxn$metHitCount, resultList$met2rxn$rxn_source_id, resultList$met2rxn$substrate_product),]
    colnames(resultList$met2rxn) <- c("reactionId", "metSourceId", "substrateProductFlag", "isCofactor", "metName", 'isTransport', "label", "direction", "equation", "htmlEquation", "ecNumber",
                                      "hasHumanProtein", "onlyHumanMets", "metHitCount")
    resultList$met2rxn <- resultList$met2rxn[,c(1,14,2, 5, 3, 4, 7, 9, 10, 8, 11, 6, 12:13)]
    if(includeRxnURLs) {
      resultList$met2rxn$rheaRxnURL <- getReactionRheaURLs(unlist(resultList$met2rxn$reactionId))
    }
  }

  if(nrow(resultList$prot2rxn) > 0) {
    colnames(resultList$prot2rxn) <- c("uniprot", "proteinName", "reactionId", "isTransport", "label", "direction", "equation", "htmlEquation", "ecNumber", "hasHumanProtein", "onlyHumanMets")
    resultList$prot2rxn <- resultList$prot2rxn[,c(3, 1:2, 5, 7, 6, 9, 4, 10, 11)]
    if(includeRxnURLs) {
      resultList$prot2rxn$rheaRxnURL <- getReactionRheaURLs(unlist(resultList$prot2rxn$reactionId))
    }
  }

  if(nrow(resultList$metProteinCommonReactions) > 0) {
    colnames(resultList$metProteinCommonReactions) <- c("metabolites", "metNames", "uniprot", "proteinName", "reactionId", "direction", "label", "htmlLabel", "isTransport", "hasHumanProtein", "onlyHumanMets")
    resultList$metProteinCommonReactions <- resultList$metProteinCommonReactions[,c(5, 1:4, 7:8, 6, 9:11)]
    if(includeRxnURLs) {
      resultList$metProteinCommonReactions$rheaRxnURL <- getReactionRheaURLs(unlist(resultList$metProteinCommonReactions$reactionId))
    }
  }

  met2rxn <- resultList$met2rxn
  prot2rxn <- resultList$prot2rxn
  if(any(met2rxn$metHitCount>1)==TRUE)
  {
    met2rxn_edit <- met2rxn[which(met2rxn$metHitCount>1),]
    met2rxn <- met2rxn[which(met2rxn$metHitCount==1),]

    met2rxn_edit <- split(met2rxn_edit, met2rxn_edit$reactionId)

    for (i in 1:length(met2rxn_edit))
    {
      metSourceId <- paste0(c(met2rxn_edit[[i]]$metSourceId), collapse = "|")
      metName <- paste0(c(met2rxn_edit[[i]]$metName), collapse = "|")
      line2add <- met2rxn_edit[[i]][1,]
      line2add$metSourceId <- metSourceId
      line2add$metName <- metName

      met2rxn <- rbind(met2rxn, line2add)

    }
  }

  if(any(duplicated(prot2rxn$reactionId))==TRUE)
  {
    duplicate_reactionId <- prot2rxn$reactionId[which(duplicated(prot2rxn$reactionId))]

    prot2rxn_edit <- prot2rxn[which(is.element(prot2rxn$reactionId, duplicate_reactionId)),]
    prot2rxn <- prot2rxn[-which(is.element(prot2rxn$reactionId, duplicate_reactionId)),]
    prot2rxn_edit <- split(prot2rxn_edit, prot2rxn_edit$reactionId)

    for (i in 1:length(prot2rxn_edit))
    {
      uniprot <- paste0(c(prot2rxn_edit[[i]]$uniprot), collapse = "|")
      proteinName <- paste0(c(prot2rxn_edit[[i]]$proteinName), collapse = "|")
      line2add <- prot2rxn_edit[[i]][1,]
      line2add$uniprot <- uniprot
      line2add$proteinName <- proteinName

      prot2rxn <- rbind(prot2rxn, line2add)

    }
  }

  resultList$met2rxn <- met2rxn
  resultList$prot2rxn <- prot2rxn

  message("Finished getReactionsForAnalytes()")

  return(resultList)
}

#' getReactionClassesForAnalytes returns reactions class and EC numbers for a collection of input compound ids
#'
#' @param analytes list of analyte ids
#' @param multiRxnParticipantCount minimum number of analytes to report a reaction class, default = 1
#' @param humanProtein require reactions to have a human protein (enzyme or transporter), default True
#' @param concatResults returns all reaction class levels in one dataframe rather than a list of 3 dataframes
#' @param includeReactionIDs adds the list of reaction ids for each reaction class
#' @param useIdMapping allows one to fuzzy match input ids to chebi and uniprot to support reaction queries.
#' @param db a RaMP database object
#'
#' @return returns a three dataframes of reaction EC class information, one for each EC level
#' @examples
#' \dontrun{
#'  analytes.of.interest = c('chebi:58115', 'chebi:456215', 'chebi:58245', 'chebi:58450',
#' 'chebi:17596', 'chebi:16335', 'chebi:16750', 'chebi:172878',
#' 'chebi:62286', 'chebi:77897', 'uniprot:P30566','uniprot:P30520',
#' 'uniprot:P00568', 'uniprot:P23109', 'uniprot:P22102', 'uniprot:P15531')
#'
#' reaction.classes <- getReactionClassesForAnalytes(analytes = analytes.of.interest, db = RaMP())
#' }
#' @export
 getReactionClassesForAnalytes <- function(analytes, multiRxnParticipantCount=1, humanProtein=TRUE,
                                           concatResults=F, includeReactionIDs=F, useIdMapping = F, db = RaMP()) {

  print("Starting reaction class query...")

  if(!useIdMapping) {
    ids = checkReactionInputIds("getReactionClassesForAnalytes", analytes)
    analytes = c(ids$metIds, ids$proteinIds)

    if(length(analytes) == 0) {
      message("None of the ids are chebi or uniprot ids.")
      message("Returning empty list. Use the parameter useIDMapping = T to enable id mapping to chebi and uniprot ids.")
      return(list())
    }
  }

  # base reaction stats, reactions per reaction class
  rxnClassStats <- getReactionClassStats(db=db, humanProtein=humanProtein)
  print("Passed the getReactionClassStats")
  # base reaction stats for analytes, distinct analyte counts (proteins or mets) per reaction class
  analyteClassStats <- getReactionClassStatsOnAnalytes(db=db, humanProtein=humanProtein)

  # if we want to enforce reactions that have more than
  # one participant in the input list
  #
  # a bit more work to query direct analyte to reaction ids
  # who would do it like this? :) (eye roll)
  if(useIdMapping) {
    if(multiRxnParticipantCount > 1) {

      analyte2Rxn = getReactionsForAnalytes(db=db, analytes=analytes, humanProtein = humanProtein)
      rxnParticipantData <- getReactionParticipantCounts(analyte2Rxn, multiRxnParticipantCount)

      if(rxnParticipantData[['total_rxns_retained']] == 0) {

        message("The input analytes map to ", rxnParticipantData$total_mapped_rxns," reactions.")
        message("The input analytes to reaction mapping do not support the multiRxnParticipantCount cutoff of ",multiRxnParticipantCount, ".")
        message("Returing empty result. Check input id prefix format and possibley reduce multiRxnParticipantCount cutoff.")

        if(concatResults) {
          return(data.frame())
        }
        else {
          return(list())
        }
      }

      # reaction list is already filtered
      keeperRxns <- unlist(rxnParticipantData$merged_rxn2analyte_count$rxn_source_id)

      metResult <- db@api$getReactionsForAnalytes(analytes = analytes, analyteType = 'metabolite', useIdMapping = useIdMapping, keeperRxns = keeperRxns, humanProtein = humanProtein)
      proteinResult <-db@api$getReactionsForAnalytes(analytes = analytes, analyteType = 'gene', useIdMapping = useIdMapping, keeperRxns = keeperRxns, humanProtein = humanProtein)
    } else {
      metResult <- db@api$getReactionsForAnalytes(analytes = analytes, analyteType = 'metabolite', useIdMapping = useIdMapping, keeperRxns = NULL, humanProtein = humanProtein)
      proteinResult <-db@api$getReactionsForAnalytes(analytes = analytes, analyteType = 'gene', useIdMapping = useIdMapping, keeperRxns = NULL, humanProtein = humanProtein)
    }

    metRxns <- metResult[,c('rxn_class', 'rxn_class_hierarchy', 'met_reactions')]
    proteinRxns <- proteinResult[,c('rxn_class_hierarchy', 'protein_reactions')]

    mergedRxnData <- merge(x=metRxns, y=proteinRxns, by.x='rxn_class_hierarchy', by.y='rxn_class_hierarchy', all.x=T, all.y=T)
    mergedRxnData$met_reactions[is.na(mergedRxnData$met_reactions)] <- ""
    mergedRxnData$protein_reactions[is.na(mergedRxnData$protein_reactions)] <- ""

    if(nrow(mergedRxnData) == 0) {

      message("The input analytes map to 0 reactions.")
      message("Returing empty result. Check input id prefix format.")

      if(concatResults) {
        return(data.frame())
      }
      else {
        return(list())
      }
    }

  } else {

    # not using the source table and ID mapping
    if(multiRxnParticipantCount > 1) {

      analyte2Rxn = RaMP::getReactionsForAnalytes(db=db, analytes=analytes, humanProtein = humanProtein)
      rxnParticipantData <- getReactionParticipantCounts(analyte2Rxn, multiRxnParticipantCount)

      if(rxnParticipantData[['total_rxns_retained']] == 0) {

        message("The input analytes map to ", rxnParticipantData$total_mapped_rxns," reactions.")
        message("The input analytes to reaction mapping do not support the multiRxnParticipantCount cutoff of ",multiRxnParticipantCount, ".")
        message("Returing empty result. Check input id prefix format and possibley reduce multiRxnParticipantCount cutoff.")

        if(concatResults) {
          return(data.frame())
        }
        else {
          return(list())
        }
      }

      # reaction list is already filtered
      keeperRxns <- unlist(rxnParticipantData$merged_rxn2analyte_count$rxn_source_id)

      metResult <- db@api$getReactionsForAnalytes(analytes = analytes, analyteType = 'metabolite', useIdMapping = useIdMapping, keeperRxns = keeperRxns, humanProtein = humanProtein)
      proteinResult <-db@api$getReactionsForAnalytes(analytes = analytes, analyteType = 'gene', useIdMapping = useIdMapping, keeperRxns = keeperRxns, humanProtein = humanProtein)
    } else {
      metResult <- db@api$getReactionsForAnalytes(analytes = analytes, analyteType = 'metabolite', useIdMapping = useIdMapping, keeperRxns = NULL, humanProtein = humanProtein)
      proteinResult <-db@api$getReactionsForAnalytes(analytes = analytes, analyteType = 'gene', useIdMapping = useIdMapping, keeperRxns = NULL, humanProtein = humanProtein)
    }

    metRxns <- metResult[,c('rxn_class', 'rxn_class_hierarchy', 'met_reactions')]
    proteinRxns <- proteinResult[,c('rxn_class_hierarchy', 'protein_reactions')]

    mergedRxnData <- merge(x=metRxns, y=proteinRxns, by.x='rxn_class_hierarchy', by.y='rxn_class_hierarchy', all.x=T, all.y=T)
    mergedRxnData$met_reactions[is.na(mergedRxnData$met_reactions)] <- ""
    mergedRxnData$protein_reactions[is.na(mergedRxnData$protein_reactions)] <- ""

    if(nrow(mergedRxnData) == 0) {

      message("The input analytes map to 0 reactions.")
      message("Returing empty result. Check input id prefix format.")

      if(concatResults) {
        return(data.frame())
      }
      else {
        return(list())
      }
    }

  }


  mergedRxnData$rxn_count <- 0

  for(i in 1:nrow(mergedRxnData)) {
    mergedRxnData$rxn_count[i]<- length(union(unlist(strsplit(mergedRxnData$met_reactions[i],',')), unlist(strsplit(mergedRxnData$protein_reactions[i],','))))
  }

  combinedResult <- merge(x=metResult, y=proteinResult, by.x=c('rxn_class_hierarchy', 'rxn_class', 'rxn_class_ec', 'ec_level'), by.y=c('rxn_class_hierarchy', 'rxn_class', 'rxn_class_ec', 'ec_level'), all.x=T, all.y=T)
  combinedResult$rxn_count.x[is.na(combinedResult$rxn_count.x)] <- 0
  combinedResult$rxn_count.y[is.na(combinedResult$rxn_count.y)] <- 0

  combinedResult$met_count[is.na(combinedResult$met_count)] <- 0
  combinedResult$protein_count[is.na(combinedResult$protein_count)] <- 0
  combinedResult$protein_reactions[is.na(combinedResult$protein_reactions)] <- ""
  combinedResult$met_reactions[is.na(combinedResult$met_reactions)] <- ""

  rxnStats <- combineStringLists(combinedResult$met_reactions, combinedResult$protein_reactions)

  combinedResult <- cbind(combinedResult, rxnStats)

  combo <- combinedResult[, (!colnames(combinedResult) %in% c("rxn_count.x", "rxn_count.y", "met_count.y", "met_reactions", "protein_reactions"))]

  combo <- combo[order(combo$ec_level, combo$rxn_count, combo$met_count, combo$protein_count, decreasing = T),]

  # now we need to merge in overall reaction and analyte stats
  combo <- merge(x=combo, y=analyteClassStats[['metStats']], by.x='rxn_class_hierarchy', by.y='rxn_class_hierarchy', all.x=F, all.y=F)
  combo <- merge(x=combo, y=analyteClassStats[['protStats']], by.x='rxn_class_hierarchy', by.y='rxn_class_hierarchy', all.x=T, all.y=F)
  combo <- merge(x=combo, y=rxnClassStats, by.x=c('rxn_class', 'rxn_class_hierarchy'), by.y=c('rxn_class', 'rxn_class_hierarchy'), all.x=T, all.y=F)




  # reorder columns...
  if(includeReactionIDs) {
    combo <- combo[,c(1,2,4,5,6,13,7,16,8,19,9)]
    colnames(combo) = c("rxnClass", "rxnClassHierarcy", "ecNumber", "classLevel", "metCount", "totalMetsInRxnClass",
                        "proteinCount", "totalProteinsInRxnClass", "reactionCount", "totalRxnsInClass", "reactionIds")
  } else {
    combo <- combo[,c(1,2,4,5,6,13,7,16,8,19)]
    colnames(combo) = c("rxnClass", "rxnClassHierarcy", "ecNumber", "classLevel", "metCount", "totalMetsInRxnClass",
                        "proteinCount", "totalProteinsInRxnClass", "reactionCount", "totalRxnsInClass")
  }


  # need to include level 1 stats for all classes
  levelOneClasses = unlist(combo[combo$classLevel == 1,]$rxnClass)
  mStats <- analyteClassStats$metStats
  pStats <- analyteClassStats$protStats

  mStats <- mStats[!(mStats$rxn_class %in% levelOneClasses) & mStats$stat_ec_level == 1,]
  pStats <- pStats[!(pStats$rxn_class %in% levelOneClasses) & pStats$stat_ec_level == 1,]
  rStats <- rxnClassStats[!(rxnClassStats$rxn_class %in% levelOneClasses) & rxnClassStats$stat_ec_level == 1,]

  # only add extra stats if there are level 1 classes that are not hit by analytes
  if(nrow(mStats > 0)) {
    extraStats <- merge(x=mStats, y=pStats, by.x="rxn_class_hierarchy", by.y="rxn_class_hierarchy", all.x=T, all.y=T, no.dups=T)
    extraStats <- merge(x=extraStats, y=rStats, by.x="rxn_class_hierarchy", by.y="rxn_class_hierarchy", all.x=T, all.y=T, no.dups=T)

    extraStats <- extraStats[,c(2,1,4,3,5,9,13)]

    extraStats$metCount <- 0
    extraStats$proteinCount <- 0
    extraStats$reactionCount <- 0

    if(includeReactionIDs) {
      extraStats$reactionIds <- ""

      extraStats <- extraStats[,c(1:4,8,5,9,6,10,7,11)]
      colnames(extraStats) <- c(c("rxnClass", "rxnClassHierarcy", "ecNumber", "classLevel", "metCount", "totalMetsInRxnClass",
                                  "proteinCount", "totalProteinsInRxnClass", "reactionCount", "totalRxnsInClass", "reactionIds"))
    } else {
      extraStats <- extraStats[,c(1:4,8,5,9,6,10,7)]
      colnames(extraStats) <- c(c("rxnClass", "rxnClassHierarcy", "ecNumber", "classLevel", "metCount", "totalMetsInRxnClass",
                                  "proteinCount", "totalProteinsInRxnClass", "reactionCount", "totalRxnsInClass"))
    }

    combo <- rbind(combo, extraStats)
  }


#  analyteNullStats <- analyteClassStats$metStats[analyteClassStats$metStats$rxn_class]
#  nullClassStats <-

  c2 <- combo
  c3 <- c2[,1:3]
  c4 <- c2[,4:10]
  c4 <- data.frame(sapply(c4, as.integer))

  combo[,4:10] <- sapply(combo[,4:10], as.integer)


  if(!concatResults) {
    result <- list()
    result[['class_ec_level_1']] <- combo[combo$classLevel == 1,]
    result[['class_ec_level_2']] <- combo[combo$classLevel == 2,]
    result[['class_ec_level_3']] <- combo[combo$classLevel == 3,]
    result[['class_ec_level_4']] <- combo[combo$classLevel == 4,]
  } else {
    result = combo
    result = result[with(result, order(classLevel, -reactionCount)),]
  }

  print("Completed reaction class query...")

  return(result)
}


#' getReactionParticipants returns protein information for a list of reaction ids.
#' This utility method can help extend information from previous queries.
#' For instance, if a user queries for reactions related to a list of metabolites,
#' this method can be used to return proteins on some subset of reaction ids to find related proteins.
#'
#' @param reactionList list of reaction ids
#' @param db a RaMP databse object
#'
#' @return returns a dataframe object with reaction to reaction participant mappings.
#' @export
getReactionParticipants <- function( reactionList = c(), db = RaMP()) {

  metResult <- db@api$getRxnMetParticipants(reactionList = reactionList)
  metResult$participant_role <- 'substrate'
  metResult$participant_role_id <- 3
  metResult$participant_role[metResult$is_product == 1] <- 'product'
  metResult$participant_role_id[metResult$is_product == 1] <- 4
  metResult$participant_role[metResult$is_cofactor == 1] <- 'cofactor'
  metResult$participant_role_id[metResult$is_cofactor == 1] <- 2

  metResult$reaction_type <- 'biochemical'

  proteinResult <- db@api$getRxnGeneParticipants(reactionList = reactionList)
  proteinResult$is_cofactor <- NA
  proteinResult$is_product <- NA
  proteinResult$iso_smiles <- NA
  proteinResult$participant_role <- 'enzyme'
  proteinResult$participant_role_id <- '1'

  proteinResult$reaction_type <- 'biochemical'

  rxnResult <- db@api$getRxnIsTransport(reactionList = reactionList)

  if(nrow(rxnResult) > 0) {
    transportRxns <- unique(unlist(rxnResult$rxn_source_id))
    proteinResult$participant_role[proteinResult$reaction_id %in% transportRxns] <- 'transporter'

    # set reaction type
    proteinResult$reaction_type[proteinResult$reaction_id %in% transportRxns] <- 'transport'
    metResult$reaction_type[metResult$reaction_id %in% transportRxns] <- 'transport'
  }

  result <- rbind(metResult, proteinResult)

  result <- result[with(result, order(reaction_id, participant_role_id)),]

  # reorder columns
  colOrder <- c(1,9, 7, 2, 3, 6)
  result <- result[,colOrder]

  colnames(result) <-  c("reactionId", "reactionType", "participantRole", "participantId", "participantName", "isoSmiles" )

  return(result)
}

#' getReactionDetails returns general reaction information for a list of reaction ids.
#' This utility methed can help extend information from previous queries.
#' For instance, if a user queries for reactions related to a list of analytes, or filtered on reactions,
#' this method can be used to return general reaction info on some subset of reaction ids of interest.
#'
#' @param reactionList list of reaction ids
#' @param db a RaMP database object
#'
#' @return returns a dataframe object with reaction information.
#' @export
getReactionDetails <- function( reactionList = c(), db = RaMP()) {
  result <- db@api$getReactionDetails(reactionIDs = reactionList)

  colnames(result) <- c("reactionId", "direction", "isTransport", "hasHumanProtein", "ecNumber", "label", "equation", "htmlEquation")

  return(result)
}


#' getAnalyteReactionAssociations returns mets associated with input proteins, and proteins associated with input metabolites
#' Associations are made through shared Rhea reactions
#'
#' @param analytes list of analyte ids
#' @param includeRheaRxnDetails returns additional columns that includes info in the reaction connecting mets and proteins
#' @param humanProtein require reactions to have a human protein (enzyme or transporter), default True
#' @param db a RaMP database object
#'
#' @return returns a three dataframes of reaction EC classe information, one for each EC level
#' @export
getRheaAnalyteReactionAssociations <- function( analytes, includeRheaRxnDetails=F, humanProtein=TRUE, db = RaMP()) {

    m2p <- getRheaEnzymesAndTransportersForMetabolites(db = db, analytes=analytes, includeRheaRxnDetails=includeRheaRxnDetails, humanProtein=humanProtein)
    p2m <- getRheaMetabolitesForProteins(db = db, analytes=analytes, includeRheaRxnDetails=includeRheaRxnDetails, humanProtein=humanProtein)

    if(nrow(m2p) > 0 && nrow(p2m) > 0) {
      result = rbind(m2p, p2m)
    } else if(nrow(m2p) > 0) {
      result = m2p
    } else {
      result = p2m
    }

    if(nrow(result) > 0) {
      if(!includeRheaRxnDetails) {
        result <- result[order(result$input_common_name, result$rxn_partner_ids),]
      } else {
        result <- result[order(result$reactionId, result$substrateProductFlag),]
      }
    }
    return(result)
}


#'Enrichment analysis for analyte-reaction class mappings
#'
#' @param analytes a vector of analyte ids (genes or metabolites) that need to be searched. ID types accepted: chebi and uniprot
#' @param humanProtein require reactions to have a human protein (enzyme or transporter), default True
#' @param alternative alternative hypothesis test passed on to fisher.test(). Options are two.sided, greater, or less (default is "less")
#' @param includeRaMPids include internal RaMP identifiers (default is "FALSE")
#' @param db a RaMP database object
#'
#' @return returns a list of [[1]] EC_Level1Stats, dataframe with columns containing reaction class ID, fisher's p value, user analytes in reaction class, and total analytes in reaction class at the enzyme class level 1 , [[2]] EC_Level2Stats, dataframe with columns containing reaction class ID, fisher's p value, user analytes in reaction class, and total analytes in reaction class at the enzyme class level 2, [[3]] analyteType,  a string specifying the type of analyte input into the function ("genes", "metabolites", or "both"), and [[4]] result_type, a string specifying enrichment of reaction classes was performed
#' @export
runEnrichReactionClass <- function( analytes,
                             humanProtein=TRUE,
                             alternative = "less",
                             includeRaMPids = FALSE,
                             db = RaMP())
{

  analytes_split <- matrix(unlist(strsplit(analytes,':')), ncol = 2L, byrow = TRUE)
  analytes_split <- split(data.frame(analytes_split), analytes_split[,1])

  if (length(analytes_split)>2)
  {
    stop("Please make sure that the input is contains only chebi and uniprot ids")
  }

  approved_ids <- c("chebi", "uniprot")

  if (length(analytes_split) == 2 & length(intersect(c("chebi", "uniprot"),
                                                     names(analytes_split)))!=2)
  {
    stop("Please make sure that the input is contains only chebi and uniprot ids")
  }

  if (all(names(analytes_split) %in% approved_ids) == FALSE)
  {
    stop("Please make sure that the input is contains only chebi and uniprot ids")
  }

  if (length(analytes_split)==2)
  {
    metabAnalytes <- paste(analytes_split$chebi$X1, analytes_split$chebi$X2, sep = ":")
    protAnalytes <- paste(analytes_split$uniprot$X1, analytes_split$uniprot$X2, sep = ":")

    reactionClassdf <- getReactionClassesForAnalytes(analytes,
                                                     humanProtein=humanProtein,
                                                     db = db)

    analyte_type = "both"

    if(length(reactionClassdf) == 0)
    {
      opt <- options(show.error.messages = FALSE)
      on.exit(options(opt))
      stop()
    }

    print("Running Fisher's tests")

    ec_level_1_stats <- runFisherReaction(reactionClassdf$class_ec_level_1, metabAnalytes = metabAnalytes, protAnalytes = protAnalytes, alternative = alternative, humanProtein=humanProtein, db = db)
    ec_level_2_stats <- runFisherReaction(reactionClassdf$class_ec_level_2, metabAnalytes = metabAnalytes, protAnalytes = protAnalytes, alternative = alternative, humanProtein=humanProtein, db = db)

    ec_level_1_stats <- cbind("rxnClass" = reactionClassdf$class_ec_level_1$rxnClass,
                              "ecNumber" = reactionClassdf$class_ec_level_1$ecNumber,
                              "Pval_Metab" = ec_level_1_stats$pvals_mets,
                              "Metab_OR" = as.numeric(ec_level_1_stats$oddsratio_mets),
                              "Pval_Prot" = ec_level_1_stats$pvals_prot,
                              "Prot_OR" = as.numeric(ec_level_1_stats$oddsratio_prot))

    ec_level_2_stats <- cbind("rxnClass" = reactionClassdf$class_ec_level_2$rxnClass,
                              "ecNumber" = reactionClassdf$class_ec_level_2$ecNumber,
                              "Pval_Metab" = ec_level_2_stats$pvals_mets,
                              "Metab_OR" = as.numeric(ec_level_2_stats$oddsratio_mets),
                              "Pval_Prot" = ec_level_2_stats$pvals_prot,
                              "Prot_OR" = as.numeric(ec_level_2_stats$oddsratio_prot))

    #Calculate adjusted pvals independently for each EC level
    ec_level_1_adjusted_stats <- adjusted_stats(ec_level_1_stats, analyte_type = "both")
    ec_level_2_adjusted_stats <- adjusted_stats(ec_level_2_stats, analyte_type = "both")

  }
  else if (length(analytes_split)==1)
  {
    if (names(analytes_split[1]) == "uniprot")
    {
      analyte_type = "uniprot"
      protAnalytes <- paste(analytes_split$uniprot$X1, analytes_split$uniprot$X2, sep = ":")

      reactionClassdf <- getReactionClassesForAnalytes(analytes,
                                                       humanProtein=humanProtein,
                                                       db = db)

      if(length(reactionClassdf) == 0)
      {
        opt <- options(show.error.messages = FALSE)
        on.exit(options(opt))
        stop()
      }

      print("Running Fisher's tests")

      ec_level_1_stats <- runFisherReaction(reactionClassdf$class_ec_level_1, metabAnalytes = NULL, protAnalytes = protAnalytes, alternative = alternative, humanProtein=humanProtein, db = db)
      ec_level_2_stats <- runFisherReaction(reactionClassdf$class_ec_level_2, metabAnalytes = NULL, protAnalytes = protAnalytes, alternative = alternative, humanProtein=humanProtein, db = db)

      ec_level_1_stats <- cbind("rxnClass" = reactionClassdf$class_ec_level_1$rxnClass,
                                "ecNumber" = reactionClassdf$class_ec_level_1$ecNumber,
                                "Pval_Prot" = ec_level_1_stats$pvals_prot,
                                "Prot_OR" = as.numeric(ec_level_1_stats$oddsratio_prot))

      ec_level_2_stats <- cbind("rxnClass" = reactionClassdf$class_ec_level_2$rxnClass,
                                "ecNumber" = reactionClassdf$class_ec_level_2$ecNumber,
                                "Pval_Prot" = ec_level_2_stats$pvals_prot,
                                "Prot_OR" = as.numeric(ec_level_2_stats$oddsratio_prot))

      #Calculate adjusted pvals independently for each EC level
      ec_level_1_adjusted_stats <- adjusted_stats(ec_level_1_stats, analyte_type = "uniprot")
      ec_level_2_adjusted_stats <- adjusted_stats(ec_level_2_stats, analyte_type ="uniprot")
    }
    else if (names(analytes_split[1]) == "chebi")
    {
      analyte_type = "chebi"
      metabAnalytes <- paste(analytes_split$chebi$X1, analytes_split$chebi$X2, sep = ":")

      reactionClassdf <- getReactionClassesForAnalytes(analytes,
                                                       humanProtein=humanProtein,
                                                       db = db)
      if(length(reactionClassdf) == 0)
      {
        opt <- options(show.error.messages = FALSE)
        on.exit(options(opt))
        stop()
      }

      print("Running Fisher's tests")

      ec_level_1_stats <- runFisherReaction(reactionClassdf$class_ec_level_1, metabAnalytes = metabAnalytes, protAnalytes = NULL, alternative = alternative, humanProtein=humanProtein, db = db)
      ec_level_2_stats <- runFisherReaction(reactionClassdf$class_ec_level_2, metabAnalytes = metabAnalytes, protAnalytes = NULL, alternative = alternative, humanProtein=humanProtein, db = db)


      ec_level_1_stats <- cbind("rxnClass" = reactionClassdf$class_ec_level_1$rxnClass,
                                "ecNumber" = reactionClassdf$class_ec_level_1$ecNumber,
                                "Pval_Metab" = ec_level_1_stats$pvals_mets,
                                "Metab_OR" = as.numeric(ec_level_1_stats$oddsratio_mets))

      ec_level_2_stats <- cbind("rxnClass" = reactionClassdf$class_ec_level_2$rxnClass,
                                "ecNumber" = reactionClassdf$class_ec_level_2$ecNumber,
                                "Pval_Metab" = ec_level_2_stats$pvals_mets,
                                "Metab_OR" = as.numeric(ec_level_2_stats$oddsratio_mets))

      #Calculate adjusted pvals independently for each EC level

      ec_level_1_adjusted_stats <- adjusted_stats(ec_level_1_stats, analyte_type ="chebi")
      ec_level_2_adjusted_stats <- adjusted_stats(ec_level_2_stats, analyte_type ="chebi")
    }
  }


  return(list(EC_Level1Stats = ec_level_1_adjusted_stats, EC_Level2Stats = ec_level_2_adjusted_stats , analyteType = analyte_type, result_type = "reactionClass_enrichment"))
}


#####################################################
###################################
#
# Supporting utility methods
#
##

#' getRheaEnzymesAndTransportersForMetabolites returns proteins that are either enzymes or transporters
#' that are reaction participants with the input metabolite ids.
#'
#' @param analytes analyte id vector
#' @param includeRheaRxnDetails returns additional columns that inlucdes info in the reaction connecting analytes and proteins
#' @param humanProtein require reactions to have a human protein (enzyme or transporter), default True
#' @param db a RaMP database object
#'
#' @return returns a dataframe object with input metabolites, and reaction and associated protein information.
#' @noRd
getRheaEnzymesAndTransportersForMetabolites <- function( analytes, includeRheaRxnDetails=FALSE,
	humanProtein=TRUE, db = RaMP()) {

  chebiIds <- analytes[grepl("chebi", analytes)]

  if(length(chebiIds) == 0) {
    print("There are no ChEBI metabolite IDs in the input. Skipping metabolite to protein query step.")
    return(data.frame())
  }

  rxns <- getReactionsForAnalytes(db=db, analytes=chebiIds, humanProtein = humanProtein, includeTransportRxns = T)


  if(nrow(rxns$met2rxn) > 0) {
    reactions <- unique(unlist(rxns$met2rxn$reactionId))
    rxnParticipants <- getReactionParticipants(db=db, reactionList = reactions)
    rxnParticipants <- rxnParticipants[rxnParticipants$participantRole %in% c("enzyme", "transporter"),]

    rxns <- rxns$met2rxn[,c(2, 3, 4, 5,6, 1, 8)]

    result <- merge(x=rxns, y=rxnParticipants, by.x="reactionId", by.y="reactionId", all.x=T, all.y=T)

    if(!includeRheaRxnDetails) {
      keeperCols = c('metSourceId', 'metName', 'participantName', 'participantId')
      result <- result[, keeperCols]
      result$relation <- "met2protein"
      result <- result[,c(5,1,2,3,4)]
      colnames(result) <- c("query_relation", "input_analyte", "input_common_name", "rxn_partner_common_name", "rxn_partner_ids")
      result <- unique(result)
    } else {
      result$relation <- "met2protein"
      result <- result[,c(13, 3, 4, 10, 11, 9, 5, 1, 8, 7)]
      colnames(result) <- c("relation","inputId", "inputCommonName", "reactionPartnerId", "reactionPartnerCommonName", "partnerRole", "substrateProductFlag", "reactionId", "reactionType", "rxnEquation")
      result <- result[order(result$inputId, result$reactionPartnerCommonName, result$substrateProductFlag),]
      result <- unique(result)
    }

  } else {
    print("None of the input ids mapped to reactions in RaMP. Empty result.")
    result = data.frame()
  }

  return(result)
}


#' getRheaEnzymesAndTransportersForMetabolites returns proteins that are either enzymes or transporters
#' that are reaction participants with the input metabolite ids.
#'
#' @param analytes analyte id vector
#' @param includeRheaRxnDetails returns additional columns that inlucdes info in the reaction connecting mets and proteins
#' @param humanProtein require reactions to have a human protein (enzyme or transporter), default True
#' @param db a RaMP database object
#'
#' @return returns a dataframe object with input metabolites, and reaction and associated protein informatin.
#' @noRd
getRheaMetabolitesForProteins <- function( analytes, includeRheaRxnDetails=F, humanProtein=T, db = RaMP()) {

  uniprotIds <- analytes[grepl("uniprot", analytes)]

  if(length(uniprotIds) == 0) {
    print("There are no uniprot accessions in the input. Skipping protein to metabolite query step.")
    return(data.frame())
  }

  rxns <- getReactionsForAnalytes(db=db, analytes=uniprotIds, humanProtein = humanProtein, includeTransportRxns = T)

  if(nrow(rxns$prot2rxn) > 0) {
    reactions <- unique(unlist(rxns$prot2rxn$reactionId))
    rxnParticipants <- getReactionParticipants(db=db, reactionList = reactions)
    rxnParticipants <- rxnParticipants[rxnParticipants$participantRole %in% c("substrate", "product", "cofactor"),]

    rxns <- rxns$prot2rxn[,c(2, 3, 4, 5,6, 2, 1, 8)]

    result <- merge(x=rxns, y=rxnParticipants, by.x="reactionId", by.y="reactionId", all.x=T, all.y=T)

    result <- result[,-ncol(result)]

    if(!includeRheaRxnDetails) {
      keeperCols = c('uniprot', 'proteinName', 'participantName', 'participantId')
      result <- result[, keeperCols]
      result$relation <- "protein2met"
      result <- result[,c(5,1,2,3,4)]
      colnames(result) <- c("query_relation", "input_analyte", "input_common_name", "rxn_partner_common_name", "rxn_partner_ids")
      result <- unique(result)
    } else {
      result$relation <- "protein2met"
      result$partnerRole <- "metabolite"
      result <- result[,c(13,2,3,11,12,14,10,1,9,5)]

      colnames(result) <- c("relation","inputId", "inputCommonName", "reactionPartnerId", "reactionPartnerCommonName", "partnerRole", "substrateProductCofactor", "reactionId", "reactionType", "rxnEquation")

      # keep substrateProductCofator descriptive, not integer flag...
      # result$substrateProductFlag[result$substrateProductFlag == 'substrate'] <- 0
      # result$substrateProductFlag[result$substrateProductFlag == 'product'] <- 1
      # result$substrateProductFlag[result$substrateProductFlag == 'cofactor'] <- -1
      result <- result[order(result$inputId, result$reactionPartnerCommonName, result$substrateProductCofactor),]
      result <- unique(result)
    }

  } else {
    print("None of the input ids mapped to reactions in RaMP. Empty result.")
    result = data.frame()
  }
  return(result)
}



#' Utility method to combine two list to tally unique counts and combine lists into a string.
#' This utility script specifically supports accounting in reaction class queries.
#' @param x list1
#' @param y list2
#' @param sep what to split on
#' @noRd
 combineStringLists <- function(x, y, sep=",") {
   reactions <- c()
   rxn_count <- c()
   for(i in 1:length(x)) {
     if(x[i] != "") {
       data <- paste0(x[i], ",", y[i])
     } else {
       data <- y[i]
     }
     dataSplit <- strsplit(data, split=sep)
     dataSplit <- unlist(unique(dataSplit))

     size = length(dataSplit)
     data <- paste0(dataSplit, collapse=", ")
     reactions <- c(reactions, data)
     rxn_count <- c(rxn_count, size)
   }
   return(data.frame(rxn_count, reactions))
 }

#' getRampSourceInfoFromAnalyteIDs Utility method to extract source table information from analyte ids
#'
#' @param db RaMP Database object
#' @param analytes list of analyte ids
#' @return returns a dataframe of ramp analyte source information
#' @noRd
getRampSourceInfoFromAnalyteIDs <- function(db = RaMP(), analytes) {
  df <- db@api$getSourceInfoForAnalyteIDs(analytes)
  return(df)
}

#' Utility method that evalutates the mapping counts for analytes to reaction ids
#'
#' @param analyte2Rxn result object from getReactionsForAnalytes
#' @param minRxnParticipantCountFilter if > 1, the set of reactions is reduced to those with having this number of mapped analytes
#' @noRd
getReactionParticipantCounts <- function(analyte2Rxn, minRxnParticipantCountFilter=1) {

  metCounts <- data.frame(table(unlist(analyte2Rxn$met2rx$reactionId)))
  proteinCounts <- data.frame(table(unlist(analyte2Rxn$protein2rxn$reactionId)))

  # need to assess status, do we have mets and proteins to reactions, one or the other
  if(nrow(metCounts) > 0 && nrow(proteinCounts) > 0) {

    mergedCounts <- merge(metCounts, proteinCounts, by.x = 'Var1', by.y='Var1', all.x=T, all.y=T)
    mergedCounts$Freq.x[is.na(mergedCounts$Freq.x)] <- 0
    mergedCounts$Freq.y[is.na(mergedCounts$Freq.y)] <- 0
    mergedCounts$partCount <- mergedCounts$Freq.x + mergedCounts$Freq.y
    mergedCounts <- mergedCounts[,c(1,4)]

  } else if(nrow(metCounts) > 0) {
    mergedCounts <- metCounts
  } else if(nrow(proteinCounts) > 0) {
    mergedCounts <- proteinCounts
  } else {
    # return an empty dataframe if there are not reaction mappings
    result <- list()
    result[['merged_rxn2analyte_count']] <- data.frame()
    result[['total_mapped_rxns']] <- 0
    result[['total_rxns_retained']] <- 0
    return(result)
  }

  colnames(mergedCounts) <- c('rxn_source_id', 'partCount')
  filteredMergedCounts <- mergedCounts[mergedCounts$partCount >= minRxnParticipantCountFilter,]

  totalRxns <- nrow(mergedCounts)
  totalRxnReatained <- nrow(filteredMergedCounts)

  result <- list()
  result[['merged_rxn2analyte_count']] <- filteredMergedCounts
  result[['total_mapped_rxns']] <- totalRxns
  result[['total_rxns_retained']] <- totalRxnReatained

  return(result)
}

#' Utility method that returns the reaction source ids for a given result from getReactionsForAnalytes
#'
#' @param reactions input result from the getReactionsForAnalytes
#' @noRd
getReactionSourceIdsFromReactionQuery <- function(reactions) {
  rxns <- c()
  if(nrow(reactions$met2rxn) > 0) {
    rxns <- c(rxns, unlist(reactions$met2rxn$rxn_source_id))
  }
  if(nrow(reactions$prot2rxn) > 0) {
    rxns <- c(rxns, unlist(reactions$prot2rxn$rxn_source_id))
  }
  return(unique(rxns))
}


#' Utility method that returns the number of rhea reaction count in each reaction class.
#'
#' @param humanProtein boolean value indicating if reactions should be constrained to those having human proteins
#' @param db a RaMP database object
#' @noRd
getReactionClassStats <- function(humanProtein=TRUE, db=RaMP()) {
  result <- db@api$getReactionClassStats(humanProtein = humanProtein)
  return(result)
}

#' Utility method that returns the number of rhea reaction count in each reaction class.
#'
#' @param humanProtein boolean value indicating if reactions should be constrained to those having human proteins
#' @param db a RaMP database object
#' @noRd
getReactionClassStatsOnAnalytes <- function( humanProtein=T, db=RaMP()) {
  metRes <- db@api$getReactionClassStats(analyteType = 'metabolite', humanProtein = humanProtein)
  protRes <- db@api$getReactionClassStats(analyteType = 'gene', humanProtein = humanProtein)
  return(list(metStats=metRes, protStats=protRes))
}

#' Utility method that returns rhea reaction urls
#'
#' @param reactionList list of prefixed (rhea:) reaction ids
#' @noRd
getReactionRheaURLs <- function(reactionList) {
  baseURL = "https://www.rhea-db.org/"
  reactionList <- gsub(":","/",reactionList)
  urls <- paste0(baseURL, reactionList)
  return(urls)
}


checkReactionInputIds <- function(callingFunction, idList) {
  idCount = length(idList)

  metIds = idList[grepl("chebi:", idList)]
  proteinIds = idList[grepl("uniprot:", idList)]

  chebiCount <- length(metIds)
  uniprotCount <- length(proteinIds)

  message("Reporting Function: ", callingFunction)
  message("The input list has ",idCount," IDs.")
  message("The input list has ",chebiCount," chebi IDs.")
  message("The input list has ",uniprotCount," uniprot IDs.")
  invalidCount = idCount - (chebiCount + uniprotCount)
  if(invalidCount > 0) {
    message(invalidCount, " IDs will not be analyzed.")
    message("Reaction queries support ChEBI IDs, prefixed with 'chebi:'")
    message("and uniprot IDs, prefixed with 'uniprot:'")
  }
  ids <- list(metIds = metIds, proteinIds = proteinIds)
  return(ids)
}

#' Utility method that runs Fisher's test of individual ec level classes
#'
#' @param ecLevelDf EC level output from getReactionClassesForAnalytes
#' @param metabAnalytes metabolites input
#' @param protAnalytes protein input
#' @noRd

runFisherReaction <- function(ecLevelDf, metabAnalytes, protAnalytes, alternative = alternative, humanProtein = humanProtein, db = db)
{
  ## Initialize empty contingency table for later
  contingencyTb <- matrix(0, nrow = 2, ncol = 2)
  colnames(contingencyTb) <- c("In Reaction Class", "Not In Reaction Class")
  rownames(contingencyTb) <- c("All Analyte", "User's Analyte")

  output <- list()

  ## First Metabolites
  if (is.null(metabAnalytes) == FALSE)
  {
    pidCount <- 0
    pval_mets <- oddsratio_mets <- totinpath <- userinpath <- pidused <- c()
    for (i in 1:nrow(ecLevelDf))
    {
      pidCount <- pidCount + 1

      tot_mets = db@api$getCountOfChebiIdsInECReactions(humanProtein)
      tot_in_reactionClass <- ecLevelDf$totalMetsInRxnClass[i]
      tot_out_reactionClass <- tot_mets - tot_in_reactionClass
      user_in_reactionClass <- ecLevelDf$metCount[i]
      user_out_reactionClass <- length(metabAnalytes) - ecLevelDf$metCount[i]

      contingencyTb[1, 1] <- tot_in_reactionClass - user_in_reactionClass
      contingencyTb[1, 2] <- tot_out_reactionClass - user_out_reactionClass
      contingencyTb[2, 1] <- user_in_reactionClass
      contingencyTb[2, 2] <- user_out_reactionClass

      # Put the test into a try catch in case there's an issue, we'll have some details on the contingency matrix
      tryCatch(
        {
          result <- stats::fisher.test(contingencyTb, alternative = alternative)
        },
        error = function(e) {
          print(toString(e))
          print(i)
          print(contingencyTb)
          print(tot_in_reactionClass)
          print(tot_out_reactionClass)
          print(user_in_reactionClass)
          print(user_out_reactionClass)
        }
      )
      pval_mets <- c(pval_mets, result$p.value)
      oddsratio_mets <- c(oddsratio_mets, result$estimate)
    }

    output$pvals_mets = pval_mets
    output$oddsratio_mets = oddsratio_mets

  }

  ## Second Protiens
  if (is.null(protAnalytes) == FALSE)
  {
    pidCount <- 0
    pval_prot <- oddsratio_prot <- totinpath <- userinpath <- pidused <- c()
    for (i in 1:nrow(ecLevelDf))
    {
      pidCount <- pidCount + 1
      tot_prot = db@api$getCountOfUniprotIdsInECReactions(humanProtein)
      tot_in_reactionClass <- ecLevelDf$totalProteinsInRxnClass[i]
      tot_out_reactionClass <- tot_prot - tot_in_reactionClass
      user_in_reactionClass <- ecLevelDf$proteinCount[i]
      user_out_reactionClass <- length(protAnalytes) - ecLevelDf$proteinCount[i]

      contingencyTb[1, 1] <- tot_in_reactionClass - user_in_reactionClass
      contingencyTb[1, 2] <- tot_out_reactionClass - user_out_reactionClass
      contingencyTb[2, 1] <- user_in_reactionClass
      contingencyTb[2, 2] <- user_out_reactionClass

      # Put the test into a try catch in case there's an issue, we'll have some details on the contingency matrix
      tryCatch(
        {
          result <- stats::fisher.test(contingencyTb, alternative = alternative)
        },
        error = function(e) {
          print(toString(e))
          print(i)
          print(contingencyTb)
          print(tot_in_reactionClass)
          print(tot_out_reactionClass)
          print(user_in_reactionClass)
          print(user_out_reactionClass)
        }
      )
      pval_prot <- c(pval_prot, result$p.value)
      oddsratio_prot <- c(oddsratio_prot, result$estimate)
    }

    output$pvals_prot = pval_prot
    output$oddsratio_prot = oddsratio_prot
  }

  return(output)
}

#' Utility method that fixes statistics values from zero inputs
#'
#' @param stats dataframe of statistic results
#' @param pValFieldName column name of p-value results
#' @param orFieldName column name of odds ratio results
#' @noRd

fix_invalid_stats <- function(stats, pValFieldName, orFieldName) {
  for ( i in 1:nrow(stats))
  {
    if (stats[[pValFieldName]][i] == 1 && stats[[orFieldName]][i] == 'Inf')
    {
      stats[[pValFieldName]][i] <- NA
      stats[[orFieldName]][i] <- NA
    }
  }
  return(stats)
}


#' Utility method that adds adjusted p-values
#'
#' @param stats dataframe of statistic results
#' @param pValFieldName column name of p-value results
#' @noRd

add_adjusted_stats <- function(stats, pValFieldName) {
  fdr <- stats::p.adjust(stats[[pValFieldName]], method = "fdr")
  stats <- cbind(stats, fdr)
  colnames(stats)[ncol(stats)] <- "Pval_FDR"
  holm <- stats::p.adjust(stats[[pValFieldName]], method = "holm")
  stats <- cbind(stats, holm)
  colnames(stats)[ncol(stats)] <- "Pval_Holm"
  return(stats)
}


#' Utility method that returns adjusted p-values when input is only metabolites or only proteins
#'
#' @param reactionClassStats modified statistics output from runFisherReaction
#' @param analyte_type whether the input is metabolites or proteins
#' @noRd
adjusted_stats <- function(reactionClassStats, analyte_type)
{
  reactionClassStats <- as.data.frame(reactionClassStats)
  reactionClassStats[-1:-2] <- sapply(reactionClassStats[-1:-2],as.numeric)

  if (analyte_type == "chebi")
  {
    reactionClassStats <- fix_invalid_stats(reactionClassStats, 'Pval_Metab', 'Metab_OR')
    reactionClassStats <- add_adjusted_stats(reactionClassStats, 'Pval_Metab')

  } else if (analyte_type == "uniprot")
  {
    reactionClassStats <- fix_invalid_stats(reactionClassStats, 'Pval_Prot', 'Prot_OR')
    reactionClassStats <- add_adjusted_stats(reactionClassStats, 'Pval_Prot')

  } else if (analyte_type == "both")
  {
    reactionClassStats <- fix_invalid_stats(reactionClassStats, 'Pval_Metab', 'Metab_OR')
    reactionClassStats <- fix_invalid_stats(reactionClassStats, 'Pval_Prot', 'Prot_OR')

    # Calculate combined p-values for pathways that have both genes and metabolites
    gm <- intersect(which(!is.na(reactionClassStats$Pval_Metab)), which(!is.na(reactionClassStats$Pval_Prot)))
    combpval <- stats::pchisq(-2 * (log((reactionClassStats$Pval_Metab[gm])) + log(reactionClassStats$Pval_Prot[gm])),
                              df = 2, lower.tail = FALSE
    )

    g <- which(is.na(reactionClassStats$Pval_Metab))
    gpval <- reactionClassStats$Pval_Prot[g]
    m <- which(is.na(reactionClassStats$Pval_Prot))
    mpval <- reactionClassStats$Pval_Metab[m]

    out <- rbind(reactionClassStats[gm, ], reactionClassStats[g, ], reactionClassStats[m, ])
    out <- cbind(out, c(combpval, gpval, mpval))
    colnames(out)[ncol(out)] <- "Pval_combined"

    reactionClassStats <- add_adjusted_stats(out, 'Pval_combined')
  }

  return(reactionClassStats)
}


