#' @importFrom R6 R6Class
#' @noRd
dbHasAnalyteCommonName <- function(db) {
  query <- "PRAGMA table_info(analyte);"
  table_info <- runQuery(sql = query, db = db)
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
      sql <- 'pragma table_info(chem_props)'
      ramptypes <- runQuery(sql = sql, db = self$db)
      return (unlist(ramptypes$name))
    },
    getRxnPartnersFromMetIDs = function(metaboliteIDs) {
      idStr <- listToQueryString(ids = metaboliteIDs)
      queryFunction <- if (supportsCommonName(db = self$db)) rxnPartnersFromMetIDsQuery else rxnPartnersFromMetIDsQueryOld
      return (runQuery(sql = queryFunction(idStr), db = self$db))
    },
    getRxnPartnersFromGeneIDs = function(geneIDs) {
      idStr <- listToQueryString(ids = geneIDs)
      queryFunction <- if (supportsCommonName(db = self$db)) rxnPartnersFromGeneIDsQuery else rxnPartnersFromGeneIDsQueryOld
      return (runQuery(sql = queryFunction(idStr), db = self$db))
    },
    getRxnPartnersFromMetNames = function(metaboliteNames) {
      idStr <- listToQueryString(ids = metaboliteNames)
      queryFunction <- if (supportsCommonName(db = self$db)) rxnPartnersFromMetNamesQuery else rxnPartnersFromMetNamesQueryOld
      return (runQuery(sql = queryFunction(idStr), db = self$db))
    },
    getRxnPartnersFromGeneNames = function(geneNames) {
      idStr <- listToQueryString(ids = geneNames)
      queryFunction <- if (supportsCommonName(db = self$db)) rxnPartnersFromGeneNamesQuery else rxnPartnersFromGeneNamesQueryOld
      return (runQuery(sql = queryFunction(idStr), db = self$db))
    },
    getRheaRxnPartnersFromMetIDs = function(metaboliteIDs, onlyHumanMets=F, humanProtein=T, includeTransportRxns=F, rxnDirs=c("UN")) {
      idStr <- listToQueryString(ids = metaboliteIDs)
      query <-  if (supportsCommonName(db = self$db)) rheaRxnPartnersFromMetIDsQuery(metaboliteIDs = idStr) else rheaRxnPartnersFromMetIDsQueryOld(metaboliteIDs = idStr)
      query <- private$addConstraintsToRxnPartnersQuery(query, onlyHumanMets, humanProtein, includeTransportRxns, rxnDirs)
      df <- runQuery(query, self$db)
      return(df)
    },
    getRheaRxnPartnersFromGeneIDs = function(geneIDs, onlyHumanMets=F, humanProtein=T, includeTransportRxns=F, rxnDirs=c("UN")) {
      idStr <- listToQueryString(ids = geneIDs)
      query <- if (supportsCommonName(db = self$db)) rheaRxnPartnersFromGeneIDsQuery(idStr) else rheaRxnPartnersFromGeneIDsQueryOld(idStr)
      query <- private$addConstraintsToRxnPartnersQuery(query = query, onlyHumanMets = onlyHumanMets, humanProtein = humanProtein, includeTransportRxns = includeTransportRxns, rxnDirs = rxnDirs)
      df <- runQuery(sql = query, db = self$db)
      return(df)
    },
    getRxnMetParticipants = function(reactionList) {
      rxnString <- listToQueryString(ids = reactionList)
      queryFunction <- if (supportsCommonName(db = self$db)) rxnMetParticipantsQuery else rxnMetParticipantsQueryOld
      return (runQuery(sql = queryFunction(rxnString), db = self$db))
    },
    getRxnGeneParticipants = function(reactionList) {
      rxnString <- listToQueryString(ids = reactionList)
      queryFunction <- if (supportsCommonName(db = self$db)) rxnGeneParticipantsQuery else rxnGeneParticipantsQueryOld
      return (runQuery(sql = queryFunction(rxnString), db = self$db))
    },
    getRxnIsTransport = function(reactionList) {
      rxnString <- listToQueryString(ids = reactionList)
      return (runQuery(sql = rxnTransportQuery(rxnString), db = self$db))
    },
    getPathwayNames = function() {
      return (runQuery(sql = getPathwayNamesQuery(), db = self$db))
    },
    getMetaboliteIDTypes = function() {
      return (runQuery(sql = getMetaboliteIDTypesQuery(), db = self$db))
    },
    getGeneIDTypes = function() {
      return (runQuery(sql = getGeneIDTypesQuery(), db = self$db))
    },
    getMetaboliteClassSources = function() {
      return (runQuery(sql = getMetaboliteClassSourcesQuery(), db = self$db))
    },
    getMetaboliteClassTypes = function() {
      return (runQuery(sql = getMetaboliteClassTypesQuery(), db = self$db))
    },
    getAllMetaboliteClasses = function() {
      return (runQuery(sql = getAllMetaboliteClassesQuery(), db = self$db))
    },
    getMetaboliteClassesForType = function(classType) {
      return (runQuery(sql = getMetaboliteClassesForTypeQuery(classType = classType), db = self$db))
    },
    getOntologies = function() {
      return (runQuery(sql = getInfoFromTableQuery(table = 'ontology'), db = self$db))
    },
    getSourceDataForAnalyteIDs = function(analyteIDs) {
      return (runQuery(sql = getInfoForIDsQuery(table = 'source', matchColumn = 'sourceId', idList = analyteIDs), db = self$db))
    },
    getSourceDataForAnalyteNames = function(analyteNames) {
      return (runQuery(sql = getSourceDataForAnalyteNamesQuery(analyteNames = analyteNames), db = self$db))
    },
    getOntologiesForRampIDs = function(rampIds) {
      return (runQuery(sql = getInfoForIDsQuery(table = 'analytehasontology', matchColumn = 'rampCompoundId', idList = rampIds), db = self$db))
    },
    getOntologyData = function(rampIds) {
      return (runQuery(sql = getInfoForIDsQuery(table = 'ontology', matchColumn = 'rampOntologyId', idList = rampIds), db = self$db))
    },
    getMetabolitesForOntology = function(ontologyList) {
      queryFunction <- if (supportsCommonName(db = self$db)) getMetabolitesForOntologyQuery else getMetabolitesForOntologyQueryOld
      return (runQuery(sql = queryFunction(ontologyList = ontologyList), db = self$db))
    },
    getAnalytesFromOntology = function(biospecimen) {
      return (runQuery(sql = getAnalytesFromOntologyQuery(biospecimen = biospecimen), db = self$db))
    },
    getMetaboliteWithOntologyCount = function() {
      return (runQuery(sql = getMetaboliteWithOntologyCountQuery(), db = self$db)$count)
    },
    getRampIDsForOntologies = function(ontologyIDs) {
      return (runQuery(sql = getInfoForIDsQuery(table = 'analytehasontology', matchColumn = 'rampOntologyId', idList = ontologyIDs), db = self$db))
      },
    getMetaboliteSourceIdsForOntology = function(biospecimen) {
      return (runQuery(sql = getMetaboliteSourceIdsForOntologyQuery(biospecimen = biospecimen), db = self$db))
    },
    getChemPropsForMetabolites = function(properties = properties, metaboliteIDs = metaboliteIDs) {
      return (runQuery(sql = getChemPropsForMetabolitesQuery(properties = properties, metaboliteIDs = metaboliteIDs), db = self$db))
    },
    getReactionsForAnalytes = function(analytes, analyteType, useIdMapping, keeperRxns, humanProtein) {
      return (runQuery(sql = getReactionsForAnalytesQuery(analytes = analytes, analyteType = analyteType, useIdMapping = useIdMapping, keeperRxns = keeperRxns, humanProtein = humanProtein), db = self$db))
    },
    getReactionDetails = function(reactionIDs) {
      return (runQuery(sql = getReactionDetailsQuery(reactionIDs = reactionIDs), db = self$db))
    },
    getSourceInfoForAnalyteIDs = function(analyteIDs) {
      return (runQuery(sql = getSourceInfoForAnalyteIDsQuery(analyteIDs = analyteIDs), db = self$db))
    },
    getReactionClassStats = function(analyteType = 'all', humanProtein) {
      return (runQuery(sql = getReactionClassStatsQuery(analyteType = analyteType, humanProtein = humanProtein), db = self$db))
    },
    getAnalytesFromPathways = function(pathways, namesOrIds = 'names', match = "exact") {
      if (supportsCommonName(db = self$db)) {
        useCommonName = TRUE
      } else {
        useCommonName = FALSE
      }
      if (namesOrIds == 'names' && match == 'fuzzy') {
        df = data.frame(matrix(nrow=0, ncol=7))
        colnames(df) <- c('analyteName', 'sourceAnalyteIDs', 'geneOrCompound',
                          'pathwayName', 'pathwayId', 'pathwayCategory', 'pathwayType')

        for(p in pathways) {
          if(nchar(p)>2) {
            currSQL = getAnalytesFromPathwaysQuery(pathways = p, namesOrIds = namesOrIds, match = match, useCommonName = useCommonName)
            subdf <- runQuery(sql = currSQL, db = self$db)
            df <- rbind(df, subdf)
          }
        }
        return (df)
      }
      pathway_list = parseListArgument(idList = pathways)
      return (runQuery(sql = getAnalytesFromPathwaysQuery(pathways = pathway_list, namesOrIds = namesOrIds, match = match, useCommonName = useCommonName), db = self$db))
    },
    getAnalytePathwaysWithOntology = function(biospecimen) {
      return (runQuery(sql = getAnalytePathwaysWithOntologyQuery(biospecimen = biospecimen), db = self$db))
    },
    getRampIDsAndSourcesForPathways = function(includeSMPDB = FALSE) {
      return (runQuery(sql = getRampIDsAndSourcesForPathwaysQuery(includeSMPDB = includeSMPDB), db = self$db))
    },
    getPathwaysForAnalytes = function(analytes, namesOrIds, includeSMPDB) {
      if (supportsCommonName(db = self$db)) {
        useCommonName = TRUE
      } else {
        useCommonName = FALSE
      }
      return (runQuery(sql = getPathwaysForAnalytesQuery(analytes = analytes, namesOrIds = namesOrIds, includeSMPDB = includeSMPDB, useCommonName = useCommonName), db = self$db))
    },
    getSynonymsForAnalyte = function(rampIds) {
      return (runQuery(sql = getSynonymsForAnalyteQuery(rampIds = rampIds), db = self$db))
    },
    getPathwayFromSourceId = function(pathwaySourceIDs) {
      return (runQuery(sql = getPathwayFromSourceIdQuery(pathwaySourceIDs = pathwaySourceIDs), db = self$db))
    },
    getRaMPVersion = function(justVersion) {
      return (runQuery(sql = getRaMPVersionQuery(justVersion = justVersion), db = self$db))
    },
    getCurrentSourceVersion = function() {
      return (runQuery(sql = getCurrentSourceVersionQuery(), db = self$db))
    },
    getEntityCountsFromSources = function() {
      return (runQuery(sql = getInfoFromTableQuery('entity_status_info'), db = self$db))
    },
    getAnalyteIntersects = function(analyteType='metabolites', scope='mapped-to-pathway') {
      return (runQuery(sql = getAnalyteIntersectsQuery(analyteType=analyteType, scope=scope), db = self$db))
    },
    getSummaryData = function() {
      return (runQuery(sql = getSummaryDataQuery(), db = self$db))
    },
    getSynonymsForSynonym = function(synonymList) {
      return (runQuery(sql = getSynonymsForSynonymQuery(synonymList = synonymList), db = self$db))
    },
    getSynonymInfoForRampIDs = function(rampIds) {
      return (runQuery(sql = getInfoForIDsQuery(table = 'analytesynonym', matchColumn = 'rampId', idList = rampIds), db = self$db))
    },
    getSourceInfoForRampIDs = function(rampIds) {
      return (runQuery(sql = getInfoForIDsQuery(table = 'source', matchColumn = 'rampId', idList = rampIds), db = self$db))
    },
    getAllSourceInfoForSourceIDs = function(sourceIds) {
      return (runQuery(sql = getInfoForIDsQuery(table = 'source', matchColumn = 'sourceid', idList = sourceIds), db = self$db))
    },
    getAllPathwaysForRampIDs = function(rampIds) {
      return (runQuery(sql = getInfoForIDsQuery(table = 'analytehaspathway', matchColumn = 'rampId', idList = rampIds), db = self$db))
    },
    getAllRampIDsForAllPathwayRampIDs = function(pathwayRampIds) {
      return (runQuery(sql = getInfoForIDsQuery(table = 'analytehaspathway', matchColumn = 'pathwayRampId', idList = pathwayRampIds), db = self$db))
    },
    getPathwayInfoForRampIDs = function(pathwayRampIds) {
      return (runQuery(sql = getInfoForIDsQuery(table = 'pathway', matchColumn = 'pathwayRampId', idList = pathwayRampIds), db = self$db))
    },
    getSourceInfoFromSourceIDs = function(sourceIds) {
      return (runQuery(sql = getSourceInfoFromSourceIDsQuery(sourceIds = sourceIds), db = self$db))
    },
    getChemicalClassFromSourceIDs = function(sourceIds) {
      return (runQuery(sql = getChemicalClassFromSourceIDsQuery(sourceIds = sourceIds), db = self$db))
    },
    getRampIdsForPathways = function(pathwayType) {
      return (runQuery(sql = getRampIdsForPathwaysQuery(pathwayType = pathwayType), db = self$db))
    },
    getAnalyteCountsForPathways = function(pathwayRampIds) {
      return (runQuery(sql = getAnalyteCountsForPathwaysQuery(pathwayRampIds = pathwayRampIds), db = self$db))
    },
    getMetaboliteCountsForClasses = function() {
      return (runQuery(sql = getMetaboliteCountsForClassesQuery(), db = self$db))
    },
    getClassesForAnalytes = function(analytes, inferIdMapping, includeAnalyteName) {
      if (supportsCommonName(db = self$db)) {
        useCommonName = TRUE
      } else {
        useCommonName = FALSE
      }
      return (runQuery(sql = getClassesForAnalytesQuery(analytes = analytes, inferIdMapping = inferIdMapping, includeAnalyteName = includeAnalyteName, useCommonName = useCommonName), db = self$db))
    }
  ),
  private = list(
    addConstraintsToRxnPartnersQuery = function(query, onlyHumanMets=F, humanProtein=T, includeTransportRxns=F, rxnDirs=c("UN")) {
      if(length(rxnDirs) == 1) {
        query <- paste0(query, " and rxn.direction = '",rxnDirs[1],"'")
      } else if(length(rxnDirs)>1) {
        query <- paste0(query, " and rxn.direction in (",listToQueryString(ids = rxnDirs),")")
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

parseListArgument = function(idList) {
  if(is.character(idList)){
    if(grepl("\n",idList)[1]){
      output_list <- strsplit(idList,"\n")
      output_list <- unlist(output_list)
    } else if(grepl(",",idList)[1]){
      output_list <- strsplit(idList,"\n")
      output_list <- unlist(output_list)
    } else {
      output_list <- idList
    }
  } else if(is.data.frame(idList)){
    output_list <- unlist(idList)
  } else {
    return("Wrong Data Format")
  }
  return (formatListAsString(idList = output_list))
}

formatListAsString <- function(idList) {
  output_list <- sapply(idList,shQuote)
  output_list <- paste(output_list,collapse = ",")
  return (output_list)
}

getInfoForIDsQuery <- function(table, matchColumn, idList) {
  query_list = formatListAsString(idList = idList)
  return (paste0('select * from ', table, ' where ', matchColumn, " in (", query_list,");"))
}

getInfoFromTableQuery <- function(table) {
  return (paste0('select * from ', table))
}

getMetaboliteCountsForClassesQuery <- function() {
  return ("select class_level_name, class_name, count(1) as pop_hits from metabolite_class
                 group by class_level_name, class_name")
}

getAnalyteCountsForPathwaysQuery <- function(pathwayRampIds) {
  pwIdsStr <- listToQueryString(ids = pathwayRampIds)
  return (paste0("select pathwayRampId, count(distinct(rampId)) as analyte_count from analytehaspathway where pathwayRampId in (", pwIdsStr,") group by pathwayRampId"))
}

getRampIdsForPathwaysQuery <- function(pathwayType) {
  return (paste0("select pathwayRampId from pathway where type = '", pathwayType, "';"))
}

getChemicalClassFromSourceIDsQuery <- function(sourceIds) {
  query_list = formatListAsString(idList = sourceIds)
  return(paste0("select distinct a.ramp_id, b.sourceId, a.class_level_name, a.class_name, a.source from metabolite_class a, source b
          where b.rampId = a.ramp_id and b.sourceId in (",query_list,")"))
}

getSourceInfoFromSourceIDsQuery <- function(sourceIds) {
  query_list = formatListAsString(idList = sourceIds)
  return (paste0("select sourceId,IDtype as analytesource, rampId from source where sourceId in (",query_list,");"))
}

getSynonymsForSynonymQuery <- function(synonymList) {
  query_list = formatListAsString(idList = synonymList)
  return(paste0("select Synonym as origins, rampId from analytesynonym where Synonym COLLATE NOCASE in (",
                  query_list, ");"))
}

getSummaryDataQuery <- function() {
  return ("select data_key, data_blob from ramp_data_object")
}

getAnalyteIntersectsQuery <- function(analyteType='metabolites', scope='mapped-to-pathway') {
  if (analyteType == 'metabolites') {
    if (scope == 'global') {
      column = 'met_intersects_json'
    } else {
      column = 'met_intersects_json_pw_mapped'
    }
  } else {
    if (scope == 'global') {
      column = 'gene_intersects_json'
    } else {
      column = 'gene_intersects_json_pw_mapped'
    }
  }
  return (paste0("select ",column," from db_version where load_timestamp order by load_timestamp desc limit 1"))
}

getCurrentSourceVersionQuery <- function() {
  return ("select * from version_info where status = 'current'")
}

getRaMPVersionQuery <- function(justVersion = TRUE) {
  if(justVersion) {
    query<-"select ramp_version from db_version where load_timestamp order by load_timestamp desc limit 1"
  } else {
    query<-"select ramp_version, load_timestamp, version_notes, db_sql_url  from db_version where load_timestamp order by load_timestamp desc limit 1"
  }
  return(query)
}

getSynonymsForAnalyteQuery <- function(rampIds) {
  ramp_id_list = formatListAsString(idList = rampIds)
  return (paste0("select rampId as rampId, group_concat(distinct Synonym order by Synonym separator '; ')
     as synonyms from analytesynonym
     where rampId in (", ramp_id_list, ") group by rampId"))
}


getAnalytePathwaysWithOntologyQuery <- function(biospecimen) {
  return(paste0(
    "SELECT analytehaspathway.* from analytehasontology, ontology, analytehaspathway where ontology.commonName in ('",
    biospecimen,
    "') and analytehasontology.rampOntologyId = ontology.rampOntologyId and analytehasontology.rampCompoundId = analytehaspathway.rampId"
  ))
}

getSourceInfoForAnalyteIDsQuery <- function(analyteIDs) {
  return (paste("select distinct sourceId, rampId, geneOrCompound from source where sourceId in (",listToQueryString(ids = analyteIDs),")"))
}

getReactionDetailsQuery <- function(reactionIDs) {
  return (
    paste0('select rxn_source_id, direction, is_transport, has_human_prot, ec_num, label, equation, html_equation
           from reaction where rxn_source_id in (', listToQueryString(ids = reactionIDs),");"))
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

getSourceDataForAnalyteNamesQuery <- function(analyteNames) {
  return (paste0("select * from source where rampId in (select * from (select rampId from analytesynonym where Synonym in (", analyteNames, ")) as subquery);"))
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

getPathwayFromSourceIdQuery <- function(pathwaySourceIDs) {
  list_pathways <- sapply(pathwaySourceIDs, shQuote)
  list_pathways <- paste(list_pathways, collapse = ",")
  return (paste0("SELECT pathwayRampId, sourceId from pathway where sourceId in (",
                 list_pathways,
                 ")"))
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
  return (paste0("select c.sourceId as input_analyte, group_concat(distinct c.commonName COLLATE NOCASE) as input_common_name,
  group_concat(distinct g.commonName COLLATE NOCASE) as rxn_partner_common_name,
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
  return (paste0("select g.sourceId as input_analyte, group_concat(distinct g.commonName COLLATE NOCASE) as input_common_name,
  group_concat(distinct c.commonName COLLATE NOCASE) as rxn_partner_common_name,
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
  return (paste0("select s.Synonym as input_analyte, group_concat(distinct c.commonName COLLATE NOCASE) as input_common_name,
  group_concat(distinct g.commonName COLLATE NOCASE) as rxn_partner_common_name,
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
  return (paste0("select s.Synonym as input_analyte, group_concat(distinct g.commonName COLLATE NOCASE) as input_common_name,
  group_concat(distinct c.commonName COLLATE NOCASE) as rxn_partner_common_name,
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

buildSimpleQuery <- function(selectClauses, distinct = FALSE, tables, whereClauses = NULL, groupByClause = NULL,
                             orderByClause = NULL) {
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

getClassesForAnalytesQuery <- function(analytes, inferIdMapping, includeAnalyteName, useCommonName = TRUE) {
  analyteStr = formatListAsString(idList = analytes)
  selectClauses <- c(
    'metabolite_class.ramp_id',
    'source.sourceId'
  )
  tables <- c('metabolite_class', 'source')
  if (inferIdMapping) {
    whereClauses <- c(paste("source.sourceId in (",analyteStr,")"),
                      'source.rampId = metabolite_class.ramp_id')
  } else {
    whereClauses <- c(paste("source.sourceId in (",analyteStr,")"),
                      'source.sourceId = metabolite_class.class_source_id')
  }
  if (includeAnalyteName) {
    if (useCommonName) {
      selectClauses <- c(selectClauses, 'analyte.common_name as common_names')
      tables <- c(tables, 'analyte')
      whereClauses <- c(whereClauses, 'analyte.rampId = source.rampId')
    } else {
      selectClauses <- c(selectClauses, 'group_concat(distinct source.commonName COLLATE NOCASE) as common_names')
    }
  }
  selectClauses <- c(selectClauses,
                     'metabolite_class.class_level_name',
                     'metabolite_class.class_name',
                     'metabolite_class.source',
                     'count(distinct(metabolite_class.class_source_id)) as directIdClassHits')
  groupByClause <- 'metabolite_class.class_name, metabolite_class.class_level_name, source.sourceId, metabolite_class.ramp_Id, metabolite_class.source'
  orderByClause <- 'directIdClassHits desc'
  return(buildSimpleQuery(
    selectClauses = selectClauses,
    distinct = TRUE,
    tables = tables,
    whereClauses = whereClauses,
    groupByClause = groupByClause,
    orderByClause = orderByClause))
}

getPathwaysForAnalytesQuery <- function(analytes, namesOrIds, includeSMPDB, useCommonName = TRUE) {
  analyte_list <- formatListAsString(idList = analytes)
  selectClauses <- c(
    'p.pathwayName',
    'p.type as pathwaySource',
    'p.sourceId as pathwayId'
  )
  tables <- c('source s', 'analytehaspathway ap', 'pathway p')
  whereClauses <- c('ap.rampId = s.rampId',
                    'p.pathwayRampId = ap.pathwayRampId')
  if (!includeSMPDB) {
    whereClauses <- c(whereClauses, 'ap.pathwaySource != "hmdb"')
  }
  if (namesOrIds == 'ids') {
    selectClauses <- c(selectClauses,
                       's.sourceId as inputId')
    whereClauses <- c(whereClauses, paste0("s.sourceId in (", analyte_list, ")"))
    if (useCommonName) {
      tables <- c(tables, 'analyte')
      selectClauses <- c(selectClauses, 'analyte.common_name as commonName')
      whereClauses <- c(whereClauses, 's.rampId = analyte.rampId')
    } else {
      selectClauses <- c(selectClauses,
                         'group_concat(distinct s.commonName COLLATE NOCASE) as commonName')
    }
    groupByClause <- 'inputId, s.rampId, pathwayId, p.pathwayName, p.type, p.pathwayRampId'
  } else {
    if (useCommonName) {
      tables <- c(tables, 'analyte')
      selectClauses <- c(selectClauses, 'analyte.common_name as commonName')
      whereClauses <- c(whereClauses, 's.rampId = analyte.rampId')
    } else {
      selectClauses <- c(selectClauses,
                         'lower(asyn.Synonym) as commonName')
    }
    selectClauses <- c(selectClauses,
                       'group_concat(distinct s.sourceId COLLATE NOCASE) as sourceIds')
    tables <- c(tables, 'analytesynonym asyn')
    whereClauses <- c(whereClauses, 's.rampId = asyn.rampId', paste0("asyn.Synonym in (", analyte_list, ")"))
    groupByClause <- 'commonName, s.rampId, pathwayId, p.pathwayName, p.type, p.pathwayRampId'
  }
  selectClauses <- c(selectClauses,
                     's.rampId',
                     'p.pathwayRampId')

  orderByClause <- 'pathwayName asc'
  return (buildSimpleQuery(selectClauses = selectClauses, tables = tables, whereClauses = whereClauses,
                           groupByClause = groupByClause, orderByClause = orderByClause))
}

getRampIDsAndSourcesForPathwaysQuery <- function(includeSMPDB = FALSE) {
  selectClauses <- c('rampId', 'pathwaySource')
  tables <- c('analytehaspathway')

  if (!includeSMPDB) {
    return(buildSimpleQuery(distinct = TRUE, selectClauses = selectClauses, tables = tables, whereClauses = c("pathwaySource != 'hmdb'")))
  }
  return(buildSimpleQuery(distinct = TRUE, selectClauses = selectClauses, tables = tables))
}

getAnalytesFromPathwaysQuery <- function(pathways, namesOrIds = 'names', match = "exact", useCommonName = TRUE) {

  if (useCommonName) {
    selectClauses <- c('analyte.common_name as analyteName')
    tables <- c('analyte')
    whereClauses <- c('analyte.rampId = s.rampId')
  } else {
    selectClauses <- c('group_concat(distinct s.commonName COLLATE NOCASE) as analyteName')
    tables <- c()
    whereClauses <- c()
  }

  selectClauses <- c(
    selectClauses,
    'group_concat(distinct s.sourceId COLLATE NOCASE) as sourceAnalyteIDs',
    's.geneOrCompound as geneOrCompound',
    'p.pathwayName as pathwayName',
    'p.sourceId as pathwayId',
    'p.pathwayCategory as pathwayCategory',
    'p.type as pathwayType'
  )
  tables <- c( tables,
    'pathway p',
    'analytehaspathway ap',
    'source s'
  )
  whereClauses <- c(whereClauses,
    's.rampId = ap.rampID',
    'ap.pathwayRampId = p.pathwayRampId',
    "(p.pathwayCategory not like 'smpdb%' or p.pathwayCategory is Null)"
  )
  groupByClause = 's.rampId, p.pathwayName, p.sourceId, p.type, s.geneOrCompound'
  orderByClause = 'p.type desc, p.pathwayName asc, s.geneOrCompound asc'

  if (namesOrIds == 'names') {
    matchColumn = 'pathwayName'
  } else {
    matchColumn = 'sourceId'
    match = 'exact'
  }

  if (match == 'exact') {
    whereClauses <- c(whereClauses, paste0('p.', matchColumn, ' in (', pathways, ')'))
  } else {
    whereClauses <- c(whereClauses, paste0('p.', matchColumn, ' like "%', pathways, '%"'))
  }

  return(buildSimpleQuery(
    selectClauses = selectClauses,
    distinct = FALSE,
    tables = tables,
    whereClauses = whereClauses,
    groupByClause = groupByClause,
    orderByClause = orderByClause))
}

getReactionClassStatsQuery <- function(analyteType = 'all', humanProtein = TRUE) {
  selectClauses <- c(
    'rc.rxn_class',
    'rc.rxn_class_hierarchy',
    'rc.ec_level as stat_ec_level',
    'rc.rxn_class_ec as stat_class_ec'
  )
  tables <- c(
    'reaction_ec_class rc',
    'reaction r')
  whereClauses <- c('rc.rxn_source_id = r.rxn_source_id')
  groupByClause <- 'rc.rxn_class, rc.ec_level, rc.rxn_class_ec'

  if (humanProtein) {
    whereClauses <- c(whereClauses, 'r.has_human_prot = 1')
  }

  if (!is.null(analyteType) && analyteType == 'metabolite') {
    selectClauses <- c(selectClauses, 'count(distinct(rm.met_source_id)) as Total_in_RxnClass_Metab')
    tables <- c(tables, 'reaction2met rm')
    whereClauses <- c(whereClauses, 'rm.rxn_source_id = rc.rxn_source_id')
  } else if (!is.null(analyteType) && analyteType == 'gene') {
    selectClauses <- c(selectClauses, 'count(distinct(rp.uniprot)) as Total_in_RxnClass_Protein')
    tables <- c(tables, 'reaction2protein rp')
    whereClauses <- c(whereClauses, 'rp.rxn_source_id = rc.rxn_source_id')
  } else {
    selectClauses <- c(selectClauses, 'count(distinct(rc.rxn_source_id)) as Total_Rxns_in_Class')
  }

  query = buildSimpleQuery(
    selectClauses = selectClauses,
    distinct = FALSE,
    tables = tables,
    whereClauses = whereClauses,
    groupByClause = groupByClause)
  return(query)
}

getReactionsForAnalytesQuery <- function(analytes, analyteType, useIdMapping, keeperRxns, humanProtein) {
  analytesStr <- listToQueryString(ids = analytes)

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
    whereClauses <- c(whereClauses, paste0("r.rxn_source_id in (",listToQueryString(ids = keeperRxns),")"))
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
