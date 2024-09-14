#' @importFrom R6 R6Class

dbHasAnalyteCommonName <- function(db) {
  query <- "PRAGMA table_info(analyte);"
  table_info <- RaMP::runQuery(sql = query, db = db)
  return("common_name" %in% table_info$name)
}

setupVersionSupport <- function(db) {
  db@versionSupport[["analyte.common_name"]] <- dbHasAnalyteCommonName(db = db)
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
    getValidChemProps = function() {
      sql = 'pragma table_info(chem_props)'
      ramptypes <- runQuery(sql = sql, db = self$db)
      return (unlist(ramptypes$name))
    },
    getRxnPartnersFromMetIDs = function(metaboliteIDs) {
      queryFunction <- if (supportsCommonName(db = self$db)) rxnPartnersFromMetIDsQuery else rxnPartnersFromMetIDsQueryOld
      return (RaMP::runQuery(sql = queryFunction(metaboliteIDs), db = self$db))
    },
    getRxnPartnersFromGeneIDs = function(geneIDs) {
      queryFunction <- if (supportsCommonName(db = self$db)) rxnPartnersFromGeneIDsQuery else rxnPartnersFromGeneIDsQueryOld
      return (RaMP::runQuery(sql = queryFunction(geneIDs), db = self$db))
    },
    getRxnPartnersFromMetNames = function(metaboliteNames) {
      queryFunction <- if (supportsCommonName(db = self$db)) rxnPartnersFromMetNamesQuery else rxnPartnersFromMetNamesQueryOld
      return (RaMP::runQuery(sql = queryFunction(metaboliteNames), db = self$db))
    },
    getRxnPartnersFromGeneNames = function(geneNames) {
      queryFunction <- if (supportsCommonName(db = self$db)) rxnPartnersFromGeneNamesQuery else rxnPartnersFromGeneNamesQueryOld
      return (RaMP::runQuery(sql = queryFunction(geneNames), db = self$db))
    },
    getRheaRxnPartnersFromMetIDs = function(metaboliteIDs, onlyHumanMets=F, humanProtein=T, includeTransportRxns=F, rxnDirs=c("UN")) {
      idStr <- listToQueryString(analytes = metaboliteIDs)
      query <-  if (supportsCommonName(db = self$db)) rheaRxnPartnersFromMetIDsQuery(metaboliteIDs = idStr) else rheaRxnPartnersFromMetIDsQueryOld(metaboliteIDs = idStr)
      query <- private$addConstraintsToRxnPartnersQuery(query, onlyHumanMets, humanProtein, includeTransportRxns, rxnDirs)
      df <- RaMP::runQuery(query, self$db)
      return(df)
    },
    getRheaRxnPartnersFromGeneIDs = function(geneIDs, onlyHumanMets=F, humanProtein=T, includeTransportRxns=F, rxnDirs=c("UN")) {
      idStr <- listToQueryString(analytes = geneIDs)
      query <- if (supportsCommonName(db = self$db)) rheaRxnPartnersFromGeneIDsQuery(idStr) else rheaRxnPartnersFromGeneIDsQueryOld(idStr)
      query <- private$addConstraintsToRxnPartnersQuery(query = query, onlyHumanMets = onlyHumanMets, humanProtein = humanProtein, includeTransportRxns = includeTransportRxns, rxnDirs = rxnDirs)
      df <- RaMP::runQuery(sql = query, db = self$db)
      return(df)
    },
    getRxnMetParticipants = function(rxnString) {
      queryFunction <- if (supportsCommonName(db = self$db)) rxnMetParticipantsQuery else rxnMetParticipantsQueryOld
      return (RaMP::runQuery(sql = queryFunction(rxnString), db = self$db))
    },
    getRxnGeneParticipants = function(rxnString) {
      queryFunction <- if (supportsCommonName(db = self$db)) rxnGeneParticipantsQuery else rxnGeneParticipantsQueryOld
      return (RaMP::runQuery(sql = queryFunction(rxnString), db = self$db))
    },
    getRxnIsTransport = function(rxnString) {
      return (RaMP::runQuery(sql = rxnTransportQuery(rxnString), db = self$db))
    },
    getPathwayNames = function() {
      return (RaMP::runQuery(sql = getPathwayNamesQuery(), db = self$db))
    },
    getMetaboliteIDTypes = function() {
      return (RaMP::runQuery(sql = getMetaboliteIDTypesQuery(), db = self$db))
    },
    getGeneIDTypes = function() {
      return (RaMP::runQuery(sql = getGeneIDTypesQuery(), db = self$db))
    },
    getMetaboliteClassSources = function() {
      return (RaMP::runQuery(sql = getMetaboliteClassSourcesQuery(), db = self$db))
    },
    getMetaboliteClassTypes = function() {
      return (RaMP::runQuery(sql = getMetaboliteClassTypesQuery(), db = self$db))
    },
    getAllMetaboliteClasses = function() {
      return (RaMP::runQuery(sql = getAllMetaboliteClassesQuery(), db = self$db))
    },
    getMetaboliteClassesForType = function(classType) {
      return (RaMP::runQuery(sql = getMetaboliteClassesForTypeQuery(classType = classType), db = self$db))
    },
    getOntologies = function() {
      return (RaMP::runQuery(sql = getOntologiesQuery(), db = self$db))
    },
    getSourceDataForAnalyteIDs = function(analyteIDs) {
      return (RaMP::runQuery(sql = getSourceDataForAnalyteIDsQuery(analyteIDs = analyteIDs), db = self$db))
    },
    getSourceDataForAnalyteNames = function(analyteNames) {
      return (RaMP::runQuery(sql = getSourceDataForAnalyteNamesQuery(analyteNames = analyteNames), db = self$db))
    },
    getOntologiesForRampIDs = function(rampIDs) {
      return (RaMP::runQuery(sql = getOntologiesForRampIDsQuery(rampIDs = rampIDs), db = self$db))
    },
    getOntologyData = function(rampIDs) {
      return (RaMP::runQuery(sql = getOntologyDataQuery(rampIDs = rampIDs), db = self$db))
    },
    getMetabolitesForOntology = function(ontologyList) {
      queryFunction <- if (supportsCommonName(db = self$db)) getMetabolitesForOntologyQuery else getMetabolitesForOntologyQueryOld
      return (RaMP::runQuery(sql = queryFunction(ontologyList = ontologyList), db = self$db))
    },
    getAnalytesFromOntology = function(biospecimen) {
      return (RaMP::runQuery(sql = getAnalytesFromOntologyQuery(biospecimen = biospecimen), db = self$db))
    },
    getMetaboliteWithOntologyCount = function() {
      return (RaMP::runQuery(sql = getMetaboliteWithOntologyCountQuery(), db = self$db)$count)
    },
    getRampIDsForOntologies = function(ontologyIDs) {
      return (RaMP::runQuery(sql = getRampIDsForOntologiesQuery(ontologyIDs = ontologyIDs), db = self$db))
    },
    getMetaboliteSourceIdsForOntology = function(biospecimen) {
      return (RaMP::runQuery(sql = getMetaboliteSourceIdsForOntologyQuery(biospecimen = biospecimen), db = self$db))
    },
    getChemPropsForMetabolites = function(properties = properties, metaboliteIDs = metaboliteIDs) {
      return (RaMP::runQuery(sql = getChemPropsForMetabolitesQuery(properties = properties, metaboliteIDs = metaboliteIDs), db = self$db))
    },
    getReactionsForAnalytes = function(analytes, analyteType, useIdMapping, keeperRxns, humanProtein) {
      return (RaMP::runQuery(sql = getReactionsForAnalytesQuery(analytes = analytes, analyteType = analyteType, useIdMapping = useIdMapping, keeperRxns = keeperRxns, humanProtein = humanProtein), db = self$db))
    },
    getReactionDetails = function(reactionIDs) {
      return (RaMP::runQuery(sql = getReactionDetailsQuery(reactionIDs = reactionIDs), db = self$db))
    },
    getSourceInfoForAnalyteIDs = function(analyteIDs) {
      return (RaMP::runQuery(sql = getSourceInfoForAnalyteIDsQuery(analyteIDs = analyteIDs), db = self$db))
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
    })
)

getSourceInfoForAnalyteIDsQuery <- function(analyteIDs) {
  return (paste("select distinct sourceId, rampId, geneOrCompound from source where sourceId in (",listToQueryString(analyteIDs),")"))
}

getReactionDetailsQuery <- function(reactionIDs) {
  return (
    paste0('select rxn_source_id, direction, is_transport, has_human_prot, ec_num, label, equation, html_equation
           from reaction where rxn_source_id in (', listToQueryString(reactionIDs),");"))
}

getChemPropsForMetabolitesQuery <- function(properties, metaboliteIDs) {
  query <- paste("select ",properties," from chem_props",
               "where chem_source_id in (",metaboliteIDs,")")
  return (query)
}

getMetaboliteSourceIdsForOntologyQuery <- function(biospecimen) {
  return (paste0("select distinct s.rampId, s.sourceId from source s, analytehasontology ao, ontology o
      where o.commonName in ('", biospecimen, "') and o.rampOntologyId=ao.rampOntologyId and s.rampId = ao.rampCompoundId"))
}

getMetaboliteWithOntologyCountQuery <- function() {
  return ("select count(distinct rampCompoundId) as count from analytehasontology")
}

getAnalytesFromOntologyQuery <- function(biospecimen) {
  return (paste0("SELECT analytehasontology.*
                  from analytehasontology,
                       ontology
                  where ontology.commonName in ('",biospecimen,"')
                    and analytehasontology.rampOntologyId = ontology.rampOntologyId
                  "))
}

getMetabolitesForOntologyQuery <- function(ontologyList) {
  return (paste0("select source.rampId,
                       group_concat(distinct source.sourceId COLLATE NOCASE)   as source_ids,
                       analyte.common_name as common_names,
                       ontology.commonName,
                       ontology.HMDBOntologyType
                from source,
                     analyte,
                     analytehasontology,
                     ontology
                where analytehasontology.rampOntologyId in (select distinct rampOntologyId
                                            from ontology
                                            where commonName in (", ontologyList, "))
                  and ontology.rampOntologyId = analytehasontology.rampOntologyId
                  and source.rampId = analytehasontology.rampCompoundId
                  and analyte.rampId = source.rampId
                group by ontology.commonName, source.rampId, ontology.HMDBOntologyType"))
}
getMetabolitesForOntologyQueryOld <- function(ontologyList) {
  return (paste0("select rampId,
        group_concat(distinct s.sourceId COLLATE NOCASE) as source_ids,
        group_concat(distinct s.commonName COLLATE NOCASE) as common_names, o.commonName, o.HMDBOntologyType
        from source s, analytehasontology ao, ontology o where ao.rampOntologyId in (
        select distinct rampOntologyId from ontology where commonName in (", ontologyList, "))
        and o.rampOntologyId = ao.rampOntologyId and s.rampId = ao.rampCompoundId
        group by o.commonName, s.rampId, o.HMDBOntologyType"))
}

getOntologyDataQuery <- function(rampIDs) {
  return (paste0("select * from ontology where rampOntologyId in (", rampIDs, ");"))
}

getOntologiesForRampIDsQuery <- function(rampIDs) {
  return (paste0("select * from analytehasontology where rampCompoundId in (", rampIDs, ");"))
}

getRampIDsForOntologiesQuery <- function(ontologyIDs) {
  return (paste0("select * from analytehasontology where rampOntologyId in (", ontologyIDs, ")"))
}

getSourceDataForAnalyteNamesQuery <- function(analyteNames) {
  return (paste0("select * from source where rampId in (select * from (select rampId from analytesynonym where Synonym in (", analyteNames, ")) as subquery);"))
}

getSourceDataForAnalyteIDsQuery <- function(analyteIDs) {
  return (paste0("select * from source where sourceId in (", analyteIDs, ");"))
}

getOntologiesQuery <- function() {
  return ("select * from ontology")
}

getMetaboliteClassesForTypeQuery <- function(classType) {
  return (paste0("select class_level_name, class_name from metabolite_class where class_level_name = '",classType,"' group by class_level_name, class_name"))
}

getAllMetaboliteClassesQuery <- function() {
  return ("select class_level_name, class_name from metabolite_class group by class_level_name, class_name")
}

getMetaboliteClassTypesQuery <- function() {
  return ("select distinct(class_level_name) from metabolite_class order by class_level_name asc")
}

getMetaboliteClassSourcesQuery <- function() {
  return ("select distinct(source) from metabolite_class order by source asc")
}

getMetaboliteIDTypesQuery <- function() {
  return ("select distinct(IDtype) from source where geneOrCompound ='compound';")
}

getGeneIDTypesQuery <- function() {
  return ("select distinct(IDtype) from source where geneOrCompound ='gene';")
}

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

buildSimpleQuery <- function(selectClauses, distinct = FALSE, tables, whereClauses, groupByClause, orderByClause) {
  query <- paste("SELECT")
  if (distinct) {
    query <- paste(query, "DISTINCT")
  }
  query <- paste(query, paste(selectClauses, collapse = ","), "FROM", paste(tables, collapse = ","))
  if (!is.null(whereClauses) && length(whereClauses) > 0) {
    query <- paste(query, "WHERE", paste(whereClauses, collapse = " AND "))
  }
  if (!is.null(groupByClause) && nchar(groupByClause) > 0) {
    query <- paste(query, "GROUP BY", groupByClause)
  }
  if (!is.null(orderByClause) && nchar(orderByClause) > 0) {
    query <- paste(query, "ORDER  BY", orderByClause)
  }
  return (query)
}

getReactionsForAnalytesQuery <- function(analytes, analyteType, useIdMapping, keeperRxns, humanProtein) {
  analytesStr <- listToQueryString(analytes)

  selectClauses <- c(
    'c.rxn_class',
    'c.rxn_class_hierarchy',
    'c.rxn_class_ec',
    'c.ec_level',
    'count(distinct(r.rxn_source_id)) as rxn_count')

  tables <- c(
    'reaction_ec_class c',
    'reaction rxn'
  )

  whereClauses <- c(
    'rxn.rxn_source_id = r.rxn_source_id',
    'c.rxn_source_id = r.rxn_source_id'
  )

  groupByClause <- "c.rxn_class_hierarchy, c.rxn_class, c.rxn_class_ec, c.ec_level"
  orderByClause <- 'ec_level asc, rxn_count desc'

  if (analyteType == "metabolite") {
    if (useIdMapping) {
      selectClauses <- c(selectClauses, 'count(distinct(r.ramp_cmpd_id)) as met_count')
      whereClauses <- c(whereClauses, 'r.ramp_cmpd_id = s.rampId')
    } else {
      selectClauses <- c(selectClauses, 'count(distinct(r.met_source_id)) as met_count')
      whereClauses <- c(whereClauses, paste0("r.met_source_id in (",analytesStr,")"))
    }
    selectClauses <- c(selectClauses, 'group_concat(distinct(r.rxn_source_id)) as met_reactions')
    tables <- c(tables, 'reaction2met r')
  } else {
    if (useIdMapping) {
      selectClauses <- c(selectClauses, 'count(distinct(r.ramp_gene_id)) as protein_count')
      whereClauses <- c(whereClauses, 'r.ramp_gene_id = s.rampId')
    } else {
      selectClauses <- c(selectClauses, 'count(distinct(r.uniprot)) as protein_count')
      whereClauses <- c(whereClauses, paste0("r.uniprot in (",analytesStr,")"))
    }
    selectClauses <- c(selectClauses, 'group_concat(distinct(r.rxn_source_id)) as protein_reactions')
    tables <- c(tables, 'reaction2protein r')
  }
  if (!is.null(keeperRxns) && nchar(keeperRxns) > 0) {
    whereClauses <- c(whereClauses, paste0("r.rxn_source_id in (",listToQueryString(keeperRxns),")"))
  }
  if (useIdMapping) {
    tables <- c(tables, 'source s')
    whereClauses <- c(whereClauses, paste0("s.sourceId in (",analytesStr,")"))
  }

  if (humanProtein) {
    whereClauses <- c(whereClauses, 'rxn.has_human_prot = 1')
  }

  return(buildSimpleQuery(
    selectClauses = selectClauses,
    distinct = TRUE,
    tables = tables,
    whereClauses = whereClauses,
    groupByClause = groupByClause,
    orderByClause = orderByClause))
}
