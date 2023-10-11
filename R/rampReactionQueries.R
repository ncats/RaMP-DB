# RaMP Reaction Queries

#' getReactionsForAnalytes
#'
#' @param analytes list of analytes
#' @param analyteType analyte type, 'metabolites' (default), 'genes' or 'both'
#' @param namesOrIds indicates if input analyte list contains identifiers or analyte names
#' @param onlyHumanMets boolean to only return pathways containing only human metabolites (ChEBI ontology) (dev in progress)
#' @param humanProtein boolean to only control pathways catalyzed by a human proteins (having human Uniprot) (dev in progress)
#' @param includeTransportRxns if TRUE, returns metabolic and transport reactions
#' @param rxnDirs character vector of length > 1, specifying reaction directions to return  c("UN", "LR", "RL", "BD", "ALL"), default = c("UN").
#'
#' @return a list of reaction information on each input analyte, separate data.frame for metabolites, genes, and common reactions
#' @export
#'
#' @examples
getReactionsForAnalytes <- function(db = RaMP(), analytes, analyteType='metabolites', namesOrIds='ids', onlyHumanMets=F, humanProtein=F, includeTransportRxns=F, rxnDirs=c("UN")) {

  genes <- data.frame()
  metabolites <- data.frame()
  mdf <- data.frame()
  gdf <- data.frame()
  if(namesOrIds == 'ids') {
    analyteSourceInfo <- getRampSourceInfoFromAnalyteIDs(db = db, analytes = analytes)
    resultGenes <- analyteSourceInfo[analyteSourceInfo$geneOrCompound == 'gene',]
    resultMets <- analyteSourceInfo[analyteSourceInfo$geneOrCompound == 'compound',]

    if(!is.null(resultMets) && nrow(resultMets)>0) {

      print("Retrieving reactions for compounds")

      mets <- resultMets[,c('sourceId', 'rampId')]
      rampIds <- unique(unlist(mets$rampId))
      mdf <- getReactionsForRaMPCompoundIds(db = db, rampCompoundIds=rampIds, onlyHumanMets=onlyHumanMets, humanProtein=humanProtein, includeTransportRxns=includeTransportRxns, rxnDirs=rxnDirs)
      if(nrow(mdf) > 0) {
        mdf <- merge(mets, mdf, by.x='rampId', by.y='ramp_cmpd_id')
      }
    }

    if(!is.null(resultGenes) && nrow(resultGenes)>0) {

      print("Retrieving reactions for genes/proteins")

      genes <- resultGenes[,c('sourceId', 'rampId')]
      rampIds <- unique(unlist(genes$rampId))
      gdf <- getReactionsForRaMPGeneIds(rampGeneIds=rampIds, onlyHumanMets=onlyHumanMets, humanProtein=humanProtein, includeTransportRxns)
      if(nrow(gdf) > 0) {
        gdf <- merge(genes, gdf, by.x='rampId', by.y='ramp_gene_id')
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
    mRampRxnIds <- unlist(mdf$ramp_rxn_id)
    gRampRxnIds <- unlist(gdf$ramp_rxn_id)
    commonRxnIds <- intersect(mRampRxnIds, gRampRxnIds)
    if(length(commonRxnIds) > 0) {
      mRxn <- mdf[mdf$ramp_rxn_id %in% commonRxnIds,]
      gRxn <- gdf[gdf$ramp_rxn_id %in% commonRxnIds,]

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
        commonRxnList <- c(commonRxnList, paste0(unique(mRxn$rxn_source_id[mRxn$ramp_rxn_id == rxn]), collapse = ','))
        mRxnSourceIds <- c(mRxnSourceIds, paste0(unique(mRxn$sourceId[mRxn$ramp_rxn_id == rxn]), collapse = ','))
        mRxnSourceNames <- c(mRxnSourceNames, paste0(unique(mRxn$met_name[mRxn$ramp_rxn_id == rxn]), collapse = ','))

        gRxnSourceIds <- c(gRxnSourceIds, paste0(unique(gRxn$sourceId[gRxn$ramp_rxn_id == rxn]), collapse = ','))
        gRxnUniprot <- c(gRxnUniprot, paste0(unique(gRxn$uniprot[gRxn$ramp_rxn_id == rxn]), collapse = ','))
        gRxnSourceName <- c(gRxnSourceName, paste0(unique(gRxn$protein_name[gRxn$ramp_rxn_id == rxn]), collapse = ','))

        rxnDir <- c(rxnDir, paste0(unique(mRxn$direction[mRxn$ramp_rxn_id == rxn]), collapse = ','))

        commonRxnLabel <- c(commonRxnLabel, paste0(unique(mRxn$label[mRxn$ramp_rxn_id == rxn]), collapse = ','))
        commonRxnHTMLEq <- c(commonRxnHTMLEq, paste0(unique(mRxn$html_equation[mRxn$ramp_rxn_id == rxn]), collapse = ','))

        isTransport <- c(isTransport, paste0(unique(mRxn$is_transport[mRxn$ramp_rxn_id == rxn]), collapse = ','))
        hasHumanProtein <- c(hasHumanProtein, paste0(unique(mRxn$has_human_prot[mRxn$ramp_rxn_id == rxn]), collapse = ','))
        onlyHumanMets <- c(onlyHumanMets, paste0(unique(mRxn$only_human_mets[mRxn$ramp_rxn_id == rxn]), collapse = ','))
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
#' @examples
getReactionsForRaMPCompoundIds <- function(db = RaMP(), rampCompoundIds, onlyHumanMets=F, humanProtein=F, includeTransportRxns=F, rxnDirs=c("UN")) {

  idStr <- listToQueryString(rampCompoundIds)
  query <- paste0("select mr.ramp_rxn_id, mr.ramp_cmpd_id, mr.met_source_id, mr.substrate_product, mr.is_cofactor, mr.met_name,
  mr.ramp_rxn_id, rxn.rxn_source_id, rxn.is_transport, rxn.label, rxn.direction, rxn.equation, rxn.html_equation, rxn.ec_num, rxn.has_human_prot, rxn.only_human_mets
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
#' @examples
getReactionsForRaMPGeneIds <- function(db = RaMP(), rampGeneIds, onlyHumanMets=F, humanProtein=F, includeTransportRxns=F, rxnDirs=c("UN")) {

  idStr <- listToQueryString(rampGeneIds)
  query <- paste0("select gr.ramp_rxn_id, gr.ramp_gene_id, gr.uniprot, gr.protein_name,
  gr.ramp_rxn_id, rxn.rxn_source_id, rxn.is_transport, rxn.label, rxn.direction, rxn.equation, rxn.html_equation, rxn.ec_num, rxn.has_human_prot, rxn.only_human_mets
  from reaction2protein gr, reaction rxn
  where gr.ramp_gene_id in (",idStr,") and rxn.ramp_rxn_id = gr.ramp_rxn_id")

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


#####################################################
###################################
#
# general dev note... perhaps the functions below should be moved to rampQueryHelpers
#
##



#' getRampSourceInfoFromAnalyteIDs Utility method to extract source table information from analyte ids
#'
#' @param analytes list of analyte ids
#'
#' @return returns a dataframe of ramp analyte source information
#'
#' @examples
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
#' @examples
listToQueryString <- function(analytes) {
  analyteStr <- paste0("'", paste0(analytes, collapse = "','"), "'", sep="")
  return (analyteStr)
}
