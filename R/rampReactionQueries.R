# RaMP Reaction Queries

#' getReactionsForAnalytes
#'
#' @param db a RaMP databse object
#' @param analytes list of analytes
#' @param namesOrIds indicates if input analyte list contains identifiers or analyte names
#' @param onlyHumanMets boolean to only return pathways containing only human metabolites (ChEBI ontology) (dev in progress)
#' @param humanProtein boolean to only control pathways catalyzed by a human proteins (having human Uniprot) (dev in progress)
#' @param includeTransportRxns if TRUE, returns metabolic and transport reactions
#' @param rxnDirs character vector of length > 1, specifying reaction directions to return  c("UN", "LR", "RL", "BD", "ALL"), default = c("UN").
#'
#' @return a list of reaction information on each input analyte, separate data.frame for metabolites, genes, and common reactions
#' @export
#'
getReactionsForAnalytes <- function(db = RaMP(), analytes, namesOrIds='ids', onlyHumanMets=F, humanProtein=F, includeTransportRxns=F, rxnDirs=c("UN")) {

  genes <- data.frame()
  metabolites <- data.frame()
  mdf <- data.frame()
  gdf <- data.frame()
  if(namesOrIds == 'ids') {
    analyteSourceInfo <- getRampSourceInfoFromAnalyteIDs(db = db, analytes = analytes)
    resultGenes <- analyteSourceInfo[analyteSourceInfo$geneOrCompound == 'gene',]
    resultMets <- analyteSourceInfo[analyteSourceInfo$geneOrCompound == 'compound',]

    if(!is.null(resultMets) && nrow(resultMets)>0) {
      mets <- resultMets[,c('sourceId', 'rampId')]
      rampIds <- unique(unlist(mets$rampId))
      mdf <- getReactionsForRaMPCompoundIds(db = db, rampCompoundIds=rampIds, onlyHumanMets=onlyHumanMets, humanProtein=humanProtein, includeTransportRxns=includeTransportRxns, rxnDirs=rxnDirs)
      if(nrow(mdf) > 0) {
        mdf <- merge(mets, mdf, by.x='rampId', by.y='ramp_cmpd_id')
        mdf <- subset(mdf, select=-c(rampId, ramp_rxn_id))
      }
    } else {
      # name-based queries
    }

    if(!is.null(resultGenes) && nrow(resultGenes)>0) {
      genes <- resultGenes[,c('sourceId', 'rampId')]
      rampIds <- unique(unlist(genes$rampId))
      gdf <- getReactionsForRaMPGeneIds(db = db, rampGeneIds=rampIds, onlyHumanMets=onlyHumanMets, humanProtein=humanProtein, includeTransportRxns)
      if(nrow(gdf) > 0) {
        gdf <- merge(genes, gdf, by.x='rampId', by.y='ramp_gene_id')
        gdf <- subset(gdf, select=-c(rampId, ramp_rxn_id))
      }
    }

    # if(!is.null(resultMets) && ncol(resultMets)>0) {
    #   metatabolites <- resultMets[,c('sourceId', 'rampId')]
    # }

  } else { # query on names
    # do we have a helper function to get ramp ids from names?

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
      mRxnSourceIds <- c()
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
        mRxnSourceIds <- c(mRxnSourceIds, paste0(unique(mRxn$sourceId[mRxn$rxn_source_id == rxn]), collapse = ','))
        mRxnSourceNames <- c(mRxnSourceNames, paste0(unique(mRxn$met_name[mRxn$rxn_source_id == rxn]), collapse = ','))

        gRxnSourceIds <- c(gRxnSourceIds, paste0(unique(gRxn$sourceId[gRxn$rxn_source_id == rxn]), collapse = ','))
        gRxnUniprot <- c(gRxnUniprot, paste0(unique(gRxn$uniprot[gRxn$rxn_source_id == rxn]), collapse = ','))
        gRxnSourceName <- c(gRxnSourceName, paste0(unique(gRxn$protein_name[gRxn$rxn_source_id == rxn]), collapse = ','))

        rxnDir <- c(rxnDir, paste0(unique(mRxn$direction[mRxn$rxn_source_id == rxn]), collapse = ','))

        commonRxnLabel <- c(commonRxnLabel, paste0(unique(mRxn$label[mRxn$rxn_source_id == rxn]), collapse = ','))
        commonRxnHTMLEq <- c(commonRxnHTMLEq, paste0(unique(mRxn$html_equation[mRxn$rxn_source_id == rxn]), collapse = ','))

        isTransport <- as.integer(c(isTransport, paste0(unique(mRxn$is_transport[mRxn$rxn_source_id == rxn]), collapse = ',')))
        hasHumanProtein <- as.integer(c(hasHumanProtein, paste0(unique(mRxn$has_human_prot[mRxn$rxn_source_id == rxn]), collapse = ',')))
        onlyHumanMets <- as.integer(c(onlyHumanMets, paste0(unique(mRxn$only_human_mets[mRxn$rxn_source_id == rxn]), collapse = ',')))
      }
      commonReactions <- data.frame(list('metabolites' = mRxnSourceIds, 'met_names' = mRxnSourceNames ,'genes' = gRxnSourceIds,
                                         'uniprot' = gRxnUniprot, 'proteinNames' = gRxnSourceName,
                                         'reactionId' = commonRxnList, 'rxn_dir' = rxnDir,'rxn_label' = commonRxnLabel,
                                         'rxn_html_label' = commonRxnHTMLEq,
                                         'is_transport' = isTransport,
                                         'has_human_protein' = hasHumanProtein,
                                         'only_human_mets' = onlyHumanMets))

      resultList[['metProteinCommonReactions']] <- commonReactions
    }
  }

  return(resultList)
}

#' getReactionsForRaMPCompoundIds returns reactions for a collection of ramp compound ids
#'
#' @param rampCompoundIds list of ramp compound ids
#' @param onlyHumanMets boolean to only return pathways containing only human metabolites (ChEBI ontology) (dev in progress)
#' @param humanProtein boolean to only control pathways catalyzed by a human proteins (having human Uniprot) (dev in progress)
#' @param includeTransportRxns if TRUE, returns metabolic and transport reactions
#'
#' @return returns a dataframe of reaction information for each ramp compound id
#'
getReactionsForRaMPCompoundIds <- function(db = RaMP(), rampCompoundIds, onlyHumanMets=F, humanProtein=F, includeTransportRxns=F, rxnDirs=c("UN")) {

  idStr <- listToQueryString(rampCompoundIds)
  query <- paste0("select mr.ramp_rxn_id, mr.ramp_cmpd_id, mr.met_source_id, mr.substrate_product, mr.is_cofactor, mr.met_name,
  rxn.rxn_source_id, rxn.is_transport, rxn.label, rxn.direction, rxn.equation, rxn.html_equation, rxn.ec_num, rxn.has_human_prot, rxn.only_human_mets
  from reaction2met mr, reaction rxn
  where mr.ramp_cmpd_id in (",idStr,") and rxn.ramp_rxn_id = mr.ramp_rxn_id")

  if(length(rxnDirs) == 1) {
    query <- paste0(query, " and rxn.direction = '",rxnDirs[1],"'")
  } else if(length(rxnDirs)>1) {
    query <- paste0(query, " and rxn.direction in (",listToQueryString(rxnDirs),")")
  } else {
    print("rxnDirs must be of length > 0")
  }

  if(humanProtein) {
    query <- paste0(query," and rxn.has_human_prot = 1")
  }

  if(onlyHumanMets) {
    query <- paste0(query, " and rxn.only_human_mets = 1")
  }

  if(!includeTransportRxns) {
    query <- paste0(query, " and rxn.is_transport = 0")
  }

  df <- RaMP::runQuery(query, db)

  return(df)
}



#' getReactionsForRaMPGeneIds returns reactions for a collection of ramp compound ids
#'
#' @param rampGeneIds list of ramp compound ids
#' @param onlyHumanMets boolean to only return pathways containing only human metabolites (ChEBI ontology) (dev in progress)
#' @param humanProtein boolean to only control pathways catalyzed by a human proteins (having human Uniprot) (dev in progress)
#' @param includeTransportRxns if TRUE, returns metabolic and transport reactions
#'
#' @return returns a dataframe of reaction information for each ramp compound id
#'
getReactionsForRaMPGeneIds <- function(db = RaMP(), rampGeneIds, onlyHumanMets=F, humanProtein=F, includeTransportRxns=F, rxnDirs=c("UN")) {

  idStr <- listToQueryString(rampGeneIds)
  query <- paste0("select gr.ramp_rxn_id, gr.ramp_gene_id, gr.uniprot, gr.protein_name, rxn.rxn_source_id,
                  rxn.is_transport, rxn.label, rxn.direction, rxn.equation, rxn.html_equation,
                  rxn.ec_num, rxn.has_human_prot, rxn.only_human_mets from reaction2protein gr,
                  reaction rxn where gr.ramp_gene_id in (",idStr,") and rxn.ramp_rxn_id = gr.ramp_rxn_id")

  if(length(rxnDirs) == 1) {
    query <- paste0(query, " and rxn.direction = '",rxnDirs[1],"'")
  } else if(length(rxnDirs)>1) {
    query <- paste0(query, " and rxn.direction in (",listToQueryString(rxnDirs),")")
  } else {
    print("rxnDirs must be of length > 0")
  }

  if(humanProtein) {
    query <- paste0(query, " and rxn.has_human_prot = 1")
  }

  if(onlyHumanMets) {
    query <- paste0(query, " and rxn.only_human_mets = 1")
  }

  if(!includeTransportRxns) {
    query <- paste0(query, " and rxn.is_transport = 0")
  }

  df <- RaMP::runQuery(query, db)

  return(df)
}


#' getReactionClassesForAnalytes returns reactions class and EC numbers for a collection of input compound ids
#'
#' @param db a RaMP database object
#' @param analytes list of analyte ids
#' @param multiRxnParticipantCount minimum number of analytes to report a reaction class, default = 1
#' @param humanProtein require reactions to have a human protein (enzyme or transporter), default True
#'
#' @return returns a three dataframes of reaction EC classe information, one for each EC level
#' @export
 getReactionClassesForAnalytes <- function(db = RaMP(), analytes, multiRxnParticipantCount=1, humanProtein=TRUE, concatResults=F) {

  print("Starting reaction class query...")

  analytesStr <- RaMP:::listToQueryString(analytes)

  # if we want to enforce reactions that have more than
  # one participant in the input list
  #
  # a bit more work to query direct analyte to reaction ids
  # who would do it like this? :) (eye roll)
  if(multiRxnParticipantCount > 1) {

    analyte2Rxn = RaMP::getReactionsForAnalytes(db=db, analytes=analytes, humanProtein = humanProtein)
    rxnParticipantData <- RaMP:::getReactionParticpantCounts(analyte2Rxn, multiRxnParticipantCount)

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
    keeperRxnStr <- RaMP:::listToQueryString(keeperRxns)

    metQuery <- paste0("select distinct c.rxn_class_hierarchy, c.rxn_class, c.rxn_class_ec, c.ec_level, count(distinct(r.rxn_source_id)) as rxn_count, count(distinct(r.ramp_cmpd_id)) as met_count, group_concat(distinct(r.rxn_source_id)) as met_reactions
        from reaction_ec_class c, reaction2met r, source s
        where s.sourceId in (",analytesStr,") and r.rxn_source_id in (",keeperRxnStr,")
        and r.ramp_cmpd_id = s.rampId
        and c.rxn_source_id = r.rxn_source_id
        group by c.rxn_class_hierarchy, c.rxn_class, c.rxn_class_ec, c.ec_level
        order by ec_level asc, rxn_count desc")

    metResult <- RaMP::runQuery(sql=metQuery, db = db)

    proteinQuery <- paste0("select distinct c.rxn_class_hierarchy, c.rxn_class, c.rxn_class_ec, c.ec_level, count(distinct(r.rxn_source_id)) as rxn_count, count(distinct(r.ramp_gene_id)) as protein_count, group_concat(distinct(r.rxn_source_id)) as protein_reactions
          from reaction_ec_class c, reaction2protein r, source s
          where s.sourceId in (",analytesStr,") and r.rxn_source_id in (",keeperRxnStr,")
          and r.ramp_gene_id = s.rampId
          and c.rxn_source_id = r.rxn_source_id
          group by c.rxn_class_hierarchy, c.rxn_class, c.rxn_class_ec, c.ec_level
          order by ec_level asc, rxn_count desc")

    proteinResult <- RaMP::runQuery(sql=proteinQuery, db=db)

  } else {

    metQuery <- paste0("select distinct c.rxn_class_hierarchy, c.rxn_class, c.rxn_class_ec, c.ec_level, count(distinct(r.rxn_source_id)) as rxn_count,
        count(distinct(r.ramp_cmpd_id)) as met_count, group_concat(distinct(r.rxn_source_id))
        as met_reactions from reaction_ec_class c, reaction2met r, source s
        where s.sourceId in (",analytesStr,")
        and r.ramp_cmpd_id = s.rampId
        and c.rxn_source_id = r.rxn_source_id
        group by c.rxn_class_hierarchy, c.rxn_class, c.rxn_class_ec, c.ec_level
        order by ec_level asc, rxn_count desc")

    metResult <- RaMP::runQuery(sql=metQuery, db = db)

    proteinQuery <- paste0("select distinct c.rxn_class_hierarchy, c.rxn_class, c.rxn_class_ec, c.ec_level, count(distinct(r.rxn_source_id)) as rxn_count, count(distinct(r.ramp_gene_id)) as protein_count, group_concat(distinct(r.rxn_source_id)) as protein_reactions
          from reaction_ec_class c, reaction2protein r, source s
          where s.sourceId in (",analytesStr,")
          and r.ramp_gene_id = s.rampId
          and c.rxn_source_id = r.rxn_source_id
          group by c.rxn_class_hierarchy, c.rxn_class, c.rxn_class_ec, c.ec_level
          order by ec_level asc, rxn_count desc")

    proteinResult <- RaMP::runQuery(sql=proteinQuery, db=db)
  }

  metRxns <- metResult[,c('rxn_class_hierarchy', 'rxn_class', 'met_reactions')]
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

  rxnStats <- RaMP:::combineStringLists(combinedResult$met_reactions, combinedResult$protein_reactions)

  combinedResult <- cbind(combinedResult, rxnStats)

  combo <- combinedResult[, (!colnames(combinedResult) %in% c("rxn_count.x", "rxn_count.y", "met_count.y", "met_reactions", "protein_reactions"))]

  combo <- combo[order(combo$ec_level, combo$rxn_count, combo$met_count, combo$protein_count, decreasing = T),]

  #colnames(combo)[1] <- 'rxn_class'

  if(!concatResults) {
    result <- list()
    result[['class_ec_level_1']] <- combo[combo$ec_level == 1,]
    result[['class_ec_level_2']] <- combo[combo$ec_level == 2,]
    result[['class_ec_level_3']] <- combo[combo$ec_level == 3,]
  } else {
    result = combo
    result = result[order(result$ec_level, -result$rxn_count),]
  }

  print("Completed reaction class query...")

  return(result)
}


#' getReactionParticipants returns protein information for a list of reaction ids.
#' This utility method can help extend information from previous queries.
#' For instance, if a user queries for reactions related to a list of metabolites,
#' this method can be used to return proteins on some subset of reaction ids to find related proteins.
#'
#' @param db a RaMP databse object
#' @param reactionList list of reaction ids
#'
#' @return returns a dataframe object with reaction to reaction participant mappings.
#' @export
getReactionParticipants <- function(db = RaMP(), reactionList = c()) {
  reactionListStr <- listToQueryString(reactionList)


  sql = paste0('select rm.rxn_source_id as reaction_id, rm.met_source_id as participant_id, rm.met_name as participant_name, rm.is_cofactor as is_cofactor, rm.substrate_product as is_product, cp.iso_smiles as iso_smiles ',
               'from reaction2met rm, chem_props cp ',
               'where rxn_source_id in (', reactionListStr,") and cp.chem_source_id = rm.met_source_id;")

  metResult <- runQuery(sql = sql, db=db)
  metResult$participant_role <- 'substrate'
  metResult$participant_role_id <- 3
  metResult$participant_role[metResult$is_product == 1] <- 'product'
  metResult$participant_role_id[metResult$is_product == 1] <- 4
  metResult$participant_role[metResult$is_cofactor == 1] <- 'cofactor'
  metResult$participant_role_id[metResult$is_cofactor == 1] <- 2

  metResult$reaction_type <- 'biochemical'

  sql = paste0('select rxn_source_id as reaction_id, uniprot as participant_id, protein_name as participant_name from reaction2protein where rxn_source_id in (', reactionListStr,");")

  proteinResult <- runQuery(sql = sql, db=db)
  proteinResult$is_cofactor <- NA
  proteinResult$is_product <- NA
  proteinResult$iso_smiles <- NA
  proteinResult$participant_role <- 'enzyme'
  proteinResult$participant_role_id <- '1'

  proteinResult$reaction_type <- 'biochemical'

  sql = paste0('select rxn_source_id, is_transport from reaction where rxn_source_id in (',reactionListStr,') and is_transport = 1')

  rxnResult <- runQuery(sql = sql, db=db)

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

  return(result)
}

#' getReactionDetails returns general reaction information for a list of reaction ids.
#' This utility methed can help extend information from previous queries.
#' For instance, if a user queries for reactions related to a list of analytes, or filtered on reactions,
#' this method can be used to return general reaction info on some subset of reaction ids of interest.
#'
#' @param db a RaMP databse object
#' @param reactionList list of reaction ids
#'
#' @return returns a dataframe object with reaction information.
#' @export
getReactionDetails <- function(db = RaMP(), reactionList = c()) {
  reactionListStr <- listToQueryString(reactionList)
  sql = paste0('select rxn_source_id, direction, is_transport, has_human_prot, ec_num, label, equation, html_equation from reaction where rxn_source_id in (', reactionListStr,");")

  result <- runQuery(sql = sql, db=db)
  return(result)
}


#####################################################
###################################
#
# Supporting utility methods
#
##

 #' Utility method to combine two list to tally unique counts and combine lists into a string.
 #' This utility script specifically supports accounting in reaction class queries.
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
#' @param analytes list of analyte ids
#'
#' @return returns a dataframe of ramp analyte source information
#'
getRampSourceInfoFromAnalyteIDs <- function(db = RaMP(), analytes) {

  analyteStr <- listToQueryString(analytes)

  query = paste("select distinct sourceId, rampId, geneOrCompound from source where sourceId in (",analyteStr,")")

  df <- RaMP::runQuery(query, db)

  return(df)
}


#' listToQueryString utility method to convert an id list to a comma separate string, with single quoted values.
#'
#' @param analytes list of analytes (can be names or ids)
#'
#' @return comma separated list of single quoted analyte ids or names
#'
listToQueryString <- function(analytes) {
  analyteStr <- paste0("'", paste0(analytes, collapse = "','"), "'", sep="")
  return (analyteStr)
}


#' Utility method that evalutates the mapping counts for analytes to reaction ids
#'
#' @param analyte2Rxn result object from getReactionsForAnalytes
#' @param minRxnParticipantFilter if > 1, the set of reactions is reduced to those with having this number of mapped analytes
#'
getReactionParticpantCounts <- function(analyte2Rxn, minRxnParticipantCountFilter=1) {

  metCounts <- data.frame(table(unlist(analyte2Rxn$met2rx$rxn_source_id)))
  proteinCounts <- data.frame(table(unlist(analyte2Rxn$protein2rxn$rxn_source_id)))

  # need to assess status, do we have mets and proteins to reactions, one or the other
  if(nrow(metCounts) > 0 && nrow(proteinCounts) > 0) {

    mergedCounts <- merge(metCounts, proteinCounts, by.x = 'Var1', by.y='Var1', all.x=T, all.y=T)
    mergedCounts$Freq.x[is.na(mergedCounts$Freq.x)] <- 0
    mergedCounts$Freq.y[is.na(mergedCounts$Freq.y)] <- 0
    mergedCounts$partCount <- mergedCounts$Freq.x + mergeCounts$Freq.y
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
#'
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



