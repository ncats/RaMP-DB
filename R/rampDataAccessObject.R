dbHasAnalyteCommonName <- function(db) {
  query <- "PRAGMA table_info(analyte);"
  table_info <- RaMP::runQuery(query, db)
  return("common_name" %in% table_info$name)
}

setupVersionSupport <- function(db) {
  db@versionSupport[["analyte.common_name"]] <- dbHasAnalyteCommonName(db)
}

supportsCommonName <- function(db) {
  return (db@versionSupport[["analyte.common_name"]])
}

DataAccessObject <- R6::R6Class(
  "DataAccessObject",
  public = list(
    db = NULL,
    initialize = function(db = NULL) {
      self$db <- db
    },
    getRxnPartnersFromMetIDs = function(metaboliteIDs) {
      queryFunction <- if (supportsCommonName(self$db)) rxnPartnersFromMetIDsQuery else rxnPartnersFromMetIDsQueryOld
      return (RaMP::runQuery(queryFunction(metaboliteIDs), self$db))
    },
    getRxnPartnersFromGeneIDs = function(geneIDs) {
      queryFunction <- if (supportsCommonName(self$db)) rxnPartnersFromGeneIDsQuery else rxnPartnersFromGeneIDsQueryOld
      return (RaMP::runQuery(queryFunction(geneIDs), self$db))
    },
    getRxnPartnersFromMetNames = function(metaboliteNames) {
      queryFunction <- if (supportsCommonName(self$db)) rxnPartnersFromMetNamesQuery else rxnPartnersFromMetNamesQueryOld
      return (RaMP::runQuery(queryFunction(metaboliteNames), self$db))
    },
    getRxnPartnersFromGeneNames = function(geneNames) {
      queryFunction <- if (supportsCommonName(self$db)) rxnPartnersFromGeneNamesQuery else rxnPartnersFromGeneNamesQueryOld
      return (RaMP::runQuery(queryFunction(geneNames), self$db))
    },
    getRheaRxnPartnersFromMetIDs = function(metaboliteIDs, onlyHumanMets=F, humanProtein=T, includeTransportRxns=F, rxnDirs=c("UN")) {
      idStr <- listToQueryString(metaboliteIDs)
      query <-  if (supportsCommonName(self$db)) rheaRxnPartnersFromMetIDsQuery(idStr) else rheaRxnPartnersFromMetIDsQueryOld(idStr)
      query <- private$addConstraintsToRxnPartnersQuery(query, onlyHumanMets, humanProtein, includeTransportRxns, rxnDirs)
      df <- RaMP::runQuery(query, self$db)
      return(df)
    },
    getRheaRxnPartnersFromGeneIDs = function(geneIDs, onlyHumanMets=F, humanProtein=T, includeTransportRxns=F, rxnDirs=c("UN")) {
      idStr <- listToQueryString(geneIDs)
      query <- if (supportsCommonName(self$db)) rheaRxnPartnersFromGeneIDsQuery(idStr) else rheaRxnPartnersFromGeneIDsQueryOld(idStr)
      query <- private$addConstraintsToRxnPartnersQuery(query, onlyHumanMets, humanProtein, includeTransportRxns, rxnDirs)
      df <- RaMP::runQuery(query, self$db)
      return(df)
    },
    getRxnMetParticipants = function(rxnString) {
      queryFunction <- if (supportsCommonName(self$db)) rxnMetParticipantsQuery else rxnMetParticipantsQueryOld
      return (RaMP::runQuery(queryFunction(rxnString), self$db))
    },
    getRxnGeneParticipants = function(rxnString) {
      queryFunction <- if (supportsCommonName(self$db)) rxnGeneParticipantsQuery else rxnGeneParticipantsQueryOld
      return (RaMP::runQuery(queryFunction(rxnString), self$db))
    },
    getRxnIsTransport = function(rxnString) {
      return (RaMP::runQuery(rxnTransportQuery(rxnString), self$db))
    },
    getPathwayNames = function() {
      return (RaMP::runQuery(getPathwayNamesQuery(), self$db))
    }
  ),
  private = list(
    addConstraintsToRxnPartnersQuery = function(query, onlyHumanMets=F, humanProtein=T, includeTransportRxns=F, rxnDirs=c("UN")) {
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
      return(query)
    }
  )
)
getPathwayNamesQuery <- function() {
  return ("select pathwayName from pathway
            where pathwayCategory != 'smpdb3'
            order by pathwayName;")
}
rxnPartnersFromMetIDsQuery <- function(metaboliteIDs) {
  return (paste0("select cmp_source.sourceId as input_analyte,
                               cmp_analyte.common_name as input_common_name,
                               gene_analyte.common_name as rxn_partner_common_name,
                               gene_source.rampId
                      from catalyzed
                          join source gene_source on catalyzed.rampGeneId = gene_source.rampId
                          join analyte gene_analyte on gene_source.rampId = gene_analyte.rampId
                          join source cmp_source on catalyzed.rampCompoundId = cmp_source.rampId
                          join analyte cmp_analyte on cmp_source.rampId = cmp_analyte.rampId
                      where cmp_source.sourceId in (",metaboliteIDs,")
                      group by gene_source.rampId, cmp_source.sourceId"))
}
rxnPartnersFromMetIDsQueryOld <- function(metaboliteIDs) {
  return (paste0("select c.sourceId as input_analyte, group_concat(distinct c.commonName COLLATE NOCASE) as input_common_names,
  group_concat(distinct g.commonName COLLATE NOCASE) as rxn_partner_common_name,
  group_concat(distinct g.sourceId COLLATE NOCASE) as rxn_partner_ids,
  g.rampId from catalyzed r
  join source g on r.rampGeneId = g.rampId
  join source c on r.rampCompoundId = c.rampId
  where c.sourceId in (",metaboliteIDs,") group by g.rampId, c.sourceId"))
}

rxnPartnersFromGeneIDsQuery <- function(geneIDs) {
  return (paste0("select gene_source.sourceId as input_analyte,
                                gene_analyte.common_name as input_common_name,
                                cmp_analyte.common_name as rxn_partner_common_name,
                                cmp_source.rampId
                          from catalyzed
                             join source gene_source on catalyzed.rampGeneId = gene_source.rampId
                             join analyte gene_analyte on gene_source.rampId = gene_analyte.rampId
                             join source cmp_source on catalyzed.rampCompoundId = cmp_source.rampId
                             join analyte cmp_analyte on cmp_source.rampId = cmp_analyte.rampId
                          where gene_source.sourceId in (", geneIDs,")
                          group by cmp_source.rampId, gene_source.sourceId"))
}
rxnPartnersFromGeneIDsQueryOld <- function(geneIDs) {
  return (paste0("select g.sourceId as input_analyte, group_concat(distinct g.commonName COLLATE NOCASE) as input_common_names,
  group_concat(distinct c.commonName COLLATE NOCASE) as rxn_partner_common_name,
  group_concat(distinct c.sourceId COLLATE NOCASE) as rxn_partner_ids,
  c.rampId from catalyzed r
  join source g on r.rampGeneId = g.rampId
  join source c on r.rampCompoundId = c.rampId
  where g.sourceId in (", geneIDs,") group by c.rampId, g.sourceId"))
}

rxnPartnersFromMetNamesQuery <- function(metaboliteNames) {
  return (paste0("select synonym.Synonym as input_analyte,
                                group_concat(distinct cmp_analyte.common_name COLLATE NOCASE) as input_common_name,
                                group_concat(distinct gene_analyte.common_name COLLATE NOCASE) as rxn_partner_common_name,
                                 gene_source.rampId
                          from catalyzed
                                   join source gene_source on catalyzed.rampGeneId = gene_source.rampId
                                   join analyte gene_analyte on gene_source.rampId = gene_analyte.rampId
                                   join source cmp_source on catalyzed.rampCompoundId = cmp_source.rampId
                                   join analyte cmp_analyte on cmp_source.rampId = cmp_analyte.rampId
                                   join analytesynonym synonym on synonym.rampId = catalyzed.rampCompoundId
                          where synonym.Synonym in (",metaboliteNames,")
                          group by gene_source.rampId, synonym.Synonym;"))
}
rxnPartnersFromMetNamesQueryOld <- function(metaboliteNames) {
  return (paste0("select s.Synonym as input_analyte, group_concat(distinct c.commonName COLLATE NOCASE) as input_common_names,
  group_concat(distinct g.commonName COLLATE NOCASE) as rxn_partner_common_name,
  group_concat(distinct g.sourceId COLLATE NOCASE) as rxn_partner_ids,
  g.rampId from catalyzed r
  join source g on r.rampGeneId = g.rampId
  join source c on r.rampCompoundId = c.rampId
  join analytesynonym s on s.rampId = r.rampCompoundId
  where s.Synonym in (",metaboliteNames,") group by g.rampId, s.Synonym"))
}

rxnPartnersFromGeneNamesQuery <- function(geneNames) {
  return (paste0("select synonym.Synonym as input_analyte,
                               group_concat(distinct gene_analyte.common_name COLLATE NOCASE) as input_common_name,
                               group_concat(distinct cmp_analyte.common_name COLLATE NOCASE) as rxn_partner_common_name,
                               cmp_source.rampId
                        from catalyzed
                                 join source gene_source on catalyzed.rampGeneId = gene_source.rampId
                                 join analyte gene_analyte on gene_source.rampId = gene_analyte.rampId
                                 join source cmp_source on catalyzed.rampCompoundId = cmp_source.rampId
                                 join analyte cmp_analyte on cmp_source.rampId = cmp_analyte.rampId
                                 join analytesynonym synonym on synonym.rampId = catalyzed.rampGeneId
                        where synonym.Synonym in (", geneNames,")
                        group by cmp_source.rampId, synonym.Synonym;"))
}
rxnPartnersFromGeneNamesQueryOld <- function(geneNames) {
  return (paste0("select s.Synonym as input_analyte, group_concat(distinct g.commonName COLLATE NOCASE) as input_common_names,
  group_concat(distinct c.commonName COLLATE NOCASE) as rxn_partner_common_name,
  group_concat(distinct c.sourceId COLLATE NOCASE) as rxn_partner_ids,
  c.rampId from catalyzed r
  join source g on r.rampGeneId = g.rampId
  join source c on r.rampCompoundId = c.rampId
  join analytesynonym s on s.rampId = r.rampGeneId
  where s.Synonym in (", geneNames,") group by c.rampId, s.Synonym"))
}

rheaRxnPartnersFromMetIDsQuery <- function(metaboliteIDs) {
  return (paste0("select mr.met_source_id,
                         mr.substrate_product,
                         mr.is_cofactor,
                         analyte.common_name,
                         rxn.rxn_source_id,
                         rxn.is_transport,
                         rxn.label,
                         rxn.direction,
                         rxn.equation,
                         rxn.html_equation,
                         rxn.ec_num,
                         rxn.has_human_prot,
                         rxn.only_human_mets
                  from reaction2met mr,
                       reaction rxn,
                       analyte
                  where mr.met_source_id in (",metaboliteIDs,")
                    and rxn.rxn_source_id = mr.rxn_source_id
                    and mr.ramp_cmpd_id = analyte.rampId"))
}
rheaRxnPartnersFromMetIDsQueryOld <- function(metaboliteIDs) {
  return (paste0("select mr.met_source_id, mr.substrate_product, mr.is_cofactor, mr.met_name,
  rxn.rxn_source_id, rxn.is_transport, rxn.label, rxn.direction, rxn.equation, rxn.html_equation, rxn.ec_num, rxn.has_human_prot, rxn.only_human_mets
  from reaction2met mr, reaction rxn
  where mr.met_source_id in (",metaboliteIDs,") and rxn.rxn_source_id = mr.rxn_source_id"))
}

rheaRxnPartnersFromGeneIDsQuery <- function(geneIDs) {
  return (paste0("select gr.uniprot,
                         analyte.common_name,
                         rxn.rxn_source_id,
                         rxn.is_transport,
                         rxn.label,
                         rxn.direction,
                         rxn.equation,
                         rxn.html_equation,
                         rxn.ec_num,
                         rxn.has_human_prot,
                         rxn.only_human_mets
                  from reaction2protein gr,
                       reaction rxn, analyte
                  where gr.uniprot in (",geneIDs,")
                    and rxn.rxn_source_id = gr.rxn_source_id
                    and gr.ramp_gene_id = analyte.rampId"))
}
rheaRxnPartnersFromGeneIDsQueryOld <- function(geneIDs) {
  return (paste0("select gr.uniprot, gr.protein_name, rxn.rxn_source_id,
                  rxn.is_transport, rxn.label, rxn.direction, rxn.equation, rxn.html_equation,
                  rxn.ec_num, rxn.has_human_prot, rxn.only_human_mets from reaction2protein gr,
                  reaction rxn where gr.uniprot in (",geneIDs,") and rxn.rxn_source_id = gr.rxn_source_id"))
}

rxnMetParticipantsQuery <- function(rxnString) {
  return (paste0("select rm.rxn_source_id     as reaction_id,
                       rm.met_source_id     as participant_id,
                       analyte.common_name  as participant_name,
                       rm.is_cofactor       as is_cofactor,
                       rm.substrate_product as is_product,
                       cp.iso_smiles        as iso_smiles
                from reaction2met rm,
                     chem_props cp,
                     analyte
                where rxn_source_id in (", rxnString,")
                  and cp.chem_source_id = rm.met_source_id
                  and rm.ramp_cmpd_id = analyte.rampId"))
}
rxnMetParticipantsQueryOld <- function(rxnString) {
  return (paste0('select rm.rxn_source_id as reaction_id, rm.met_source_id as participant_id, rm.met_name as participant_name, rm.is_cofactor as is_cofactor, rm.substrate_product as is_product, cp.iso_smiles as iso_smiles ',
                 'from reaction2met rm, chem_props cp ',
                 'where rxn_source_id in (', rxnString,") and cp.chem_source_id = rm.met_source_id;"))
}

rxnGeneParticipantsQuery <- function(rxnString) {
  return (paste0("select rxn_source_id       as reaction_id,
                       uniprot             as participant_id,
                       analyte.common_name as participant_name
                from   reaction2protein,
                       analyte
                where  rxn_source_id in (", rxnString, ")
                       and reaction2protein.ramp_gene_id = analyte.rampId;"))
}
rxnGeneParticipantsQueryOld <- function(rxnString) {
  return (paste0('select rxn_source_id as reaction_id, uniprot as participant_id, protein_name as participant_name from reaction2protein where rxn_source_id in (', rxnString,");"))
}

rxnTransportQuery <- function(rxnString) {
  return(paste0('select rxn_source_id, is_transport from reaction where rxn_source_id in (',rxnString,') and is_transport = 1'))
}

