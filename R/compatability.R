#' Remove this and all callers for RaMP 4.0
#' @noRd
assertDBparamIsRight <- function(firstParam, dbParam) {
  functionCall <- sys.call(-1)
  firstParamIsDB <- is(firstParam, 'RaMP')
  firstParamName <- substitute(firstParam)
  dbParamIsNotDB <- !is(dbParam, 'RaMP')
  dbParamName <- substitute(dbParam)

  if (firstParamIsDB || dbParamIsNotDB) {

    error_message <- paste0(
      "As of RaMP 3.0, the RaMP object should be the last parameter, or a named parameter.\n",
      "Please update: ", deparse(functionCall), "\n",
      if (firstParamIsDB) paste0(" - ", firstParamName, " shouldn't be a RaMP object\n") else "",
      if (dbParamIsNotDB) paste0(" - ", dbParamName, " should be a RaMP object.\n") else ""
    )

    # Raise the error
    stop(error_message, call. = FALSE)
  }
}

#' For handling renamed parameters
#' @noRd
handleRenamedParameter <- function(argument, oldName, version) {
  function_call <- sys.call(-1)
  function_call_list <- as.list(function_call)
  function_string <- paste0(deparse(function_call), collapse="")
  argumentName <- substitute(argument)
  if (oldName %in% names(function_call_list)) {
    warning(paste0(
      "Please update this function call: '", function_string, "'\n",
      "As of RaMP ", version, ", '", oldName, "' has been renamed to '", argumentName, "'\n",
       "In the future this warning will become an error.\n"), call. = FALSE)
    return (eval.parent(function_call_list[[oldName]]))
  }
  return (argument)
}
#' @rdname package-deprecated
#' @title deprecated functions - to be removed in RaMP 4.0
#' @description Use plotPathwayResults instead of pathwayResultsPlot
#' @param pathwaysSig output of FilterFisherResults
#' @param pval Which p value to plot, choose from Raw, FDR or Holm-adjusted
#' @param perc_analyte_overlap Minimum overlap for pathways to be considered similar
#' (Default = 0.2)
#' @param perc_pathway_overlap Minimum overlap for clusters to merge (Default = 0.2)
#' @param min_pathway_tocluster Minimum number of 'similar' pathways required to start
#' a cluster (medoid) (Default = 3)
#' @param text_size Scales all text in figure (Default=16)
#' @param sig_cutoff Aesthetic, shows pvalue cutoff for significant pathways
#' @param interactive If TRUE, return interactive plotly object instead of ggplot object
#' @export
pathwayResultsPlot <- function(db = RaMP(), pathwaysSig, pval = "FDR", perc_analyte_overlap = 0.5,
                               perc_pathway_overlap = 0.5, min_pathway_tocluster = 3,
                               text_size = 8, sig_cutoff = 0.05, interactive=FALSE)
{
  .Deprecated('plotPathwayResults')
  return(plotPathwayResults(pathwaysSig = pathwaysSig, pVal = pval, percAnalyteOverlap = perc_analyte_overlap,
                            percPathwayOverlap = perc_pathway_overlap, minPathwayToCluster = min_pathway_tocluster,
                            textSize = text_size, sigCutoff = sig_cutoff, interactive = interactive, db = db))
}

#' @rdname package-deprecated
#' @title deprecated functions - to be removed in RaMP 4.0
#' @description Use getChemClass instead of chemicalClassSurvey
#' @param mets a list object of source prepended metaboite ids, representing a metabolite set of interest
#' @param background an optional list of source prepended metaboite ids to be used as the background reference of
#' metabolites for enrichment. The background can be either a list of ids, a file name containing the id list,
#' one id per column (no file header row) or a specificed biospecimen type (available biospecimen types: "Blood",
#' "Adipose tissue", "Heart", "Urine", "Brain", "Liver", "Kidney","Saliva", or "Feces").
#' @param background_type one of 'database' (all analytes in the RaMP Database), 'list' (a list of input ids),
#' or 'file' in which case the background parameter will be a file path, or 'biospecimen' where the specified background parameter is
#' a RaMP HMDB metabolite ontology term (see background parameter, above, for the most common biospecimen background values).
#' @param includeRaMPids include internal RaMP identifiers (default is "FALSE")
#' @param inferIdMapping if FALSE, the survey only reports on class annotations made directly on the input ids.
#' If inferIdMapping is set to TRUE, the ids are cross-referenced or mapped to related ids that contain metabolite class annotations.
#' The default is TRUE.
#' @export
chemicalClassSurvey <- function(db = RaMP(), mets, background = "database",
                                background_type="database", includeRaMPids = FALSE,
                                inferIdMapping = TRUE) {
  .Deprecated('getChemClass')
  return (getChemClass(mets = mets, background = background,
                       backgroundType = background_type, includeRaMPids = includeRaMPids,
                       inferIdMapping = inferIdMapping, db = db))
}

#' @rdname package-deprecated
#' @title deprecated functions - to be removed in RaMP 4.0
#' @description Use runEnrichChemClass instead of chemicalClassEnrichment
#' @param mets a vector of source prepended metabolite ids
#' @param background an optional list of source prepended metaboite ids to be used as the background reference of
#' metabolites for enrichment. The background can be either a list of ids, a file name containing the id list,
#' one id per column (no file header row) or a specificed biospecimen type (available biospecimen types: "Blood",
#' "Adipose tissue", "Heart", "Urine", "Brain", "Liver", "Kidney","Saliva", or "Feces").
#' @param background_type one of 'database' (all analytes in the RaMP Database), 'list' (a list of input ids),
#' or 'file' in which case the background parameter will be a file path, or 'biospecimen' where the specified background parameter is
#' a RaMP HMDB metabolite ontology term (see background parameter, above. for the most common biospecimen background values).
#' @param inferIdMapping if FALSE, the method only reports on class annotations made directly on the input ids.
#' If inferIdMapping is set to TRUE, the input ids are cross-referenced or mapped to other existing ids that contain metabolite class annotations.
#' Following id cross references can expand coverage if the input type is other than HMDB ids or LIPIDMAPS ids.
#' The default value is FALSE.
#' @export
chemicalClassEnrichment <- function(db = RaMP(), mets, background = "database",
                                    background_type = "database", inferIdMapping=F) {
  .Deprecated('runEnrichChemClass')
  return (runEnrichChemClass(mets = mets, background = background,
                             backgroundType = background_type, inferIdMapping = inferIdMapping, db = db))
}

#' @rdname package-deprecated
#' @title deprecated functions - to be removed in RaMP 4.0
#' @description Use filterEnrichResults instead of FilterFishersResults
#' @param fishers_df The data frame generated by runFisherTest
#' @param pval_type Specifies which p-value to use as the filter threshold.
#' Permitted values are 'pval' and 'fdr' for chemical class and pathway enrichment.
#' Pathway enrichment also includes an optional 'holm' value for holm p-value corrections. Default is 'fdr'.
#' @param pval_cutoff return pathways where pval_type p-values are < pval_cutoff
#' @export
FilterFishersResults <- function(fishers_df, pval_type = 'fdr', pval_cutoff = 0.1) {
  .Deprecated('filterEnrichResults')
  return (filterEnrichResults(enrichResults = fishers_df, pValType = pval_type, pValCutoff = pval_cutoff))
}

#' @rdname package-deprecated
#' @title deprecated functions - to be removed in RaMP 4.0
#' @description Use writeFishersResults instead of write_FishersResults
#' @param fishResults a data frame returned by function runCombinedFisherTest()
#' @param outputfile name of output file
#' @param rampid whether or not to include rampId (default is FALSE)
#' @export
write_FishersResults <- function(fishResults = "none", outputfile = "none", rampid = FALSE) {
  .Deprecated('writeFishersResults')
  return (writeFishersResults(fishResults = fishResults, outputFile = outputfile, includeRaMPids = rampid))
}

#' @rdname package-deprecated
#' @title deprecated functions - to be removed in RaMP 4.0
#' @description Use runEnrichPathways instead of runCombinedFisherTest
#' @param db a RaMP databse object
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
#' @param pathway_definitions If "RaMP" (default), use pathway definitions within RaMP-DB. Else, supply path to gmx file containing custom pathway definitions. GMX files are a tab-separated format that contain one analyte set per column, with the name of the set in the first row, and constituent analytes in subsequent rows. Please supply a .xls or .xlsx file. If supplying pathway definitions for genes and metabolites, ensure that metabolite definitions are on tab 1, and gene definitions are on tab2.
#' @param include_smpdb Include pathways from smpdb/hmdb in analysis. Excluded by default since definitions are highly redundant
#' @export
runCombinedFisherTest <- function(
  db = RaMP(),
  analytes,
  NameOrIds = "ids",
  total_genes = 20000,
  min_analyte = 2,
  MCall = F,
  alternative = "less",
  min_path_size = 5,
  max_path_size = 150,
  includeRaMPids = FALSE,
  background_type = "database",
  background = "database",
  pathway_definitions = "RaMP",
  include_smpdb = FALSE) {
  .Deprecated('runEnrichPathways')
  return (runEnrichPathways(
    analytes = analytes,
    namesOrIds = NameOrIds,
    totalGenes = total_genes,
    minAnalyte = min_analyte,
    alternative  = alternative,
    minPathwaySize = min_path_size,
    maxPathwaySize = max_path_size,
    includeRaMPids = includeRaMPids,
    backgroundType = background_type,
    background = background,
    pathwayDefinitions = pathway_definitions,
    includeSmpdb = include_smpdb,
    db = db
  ))
}
