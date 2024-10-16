#' Plots a network based on gene-metabolite relationships
#' @importFrom magrittr %>%
#'
#' @param catalyzeDf a data.frame output by rampFastCata() that contains analytes that are in the same reaction
#' @return  An interactive HTML plot that allows the user to pan/zoom into regions of interest. User genes/
#' metabolites are highlighted in blue, whereas analytes found by the function are orange.
#' @examples
#' \dontrun{
#' inputs.of.interest <- c("kegg:C00186" , "hmdb:HMDB0000148", "kegg:C00780", "hmdb:HMDB0000064",
#'       "ensembl:ENSG00000115850", "uniprot:Q99259")
#'
#' new.transcripts <- rampFastCata(analytes = inputs.of.interest, db = rampDB)
#'
#' plotCataNetwork(head(new.transcripts$HMDB_Analyte_Associations, n=100))
#' }
#' @export
plotCataNetwork <- function(catalyzeDf = "") {

  if(nrow(catalyzeDf) == 0) {
    message("The input data has 0 rows. plotCataNetwork function is returning without generating a plot.")
    return()
  }
  if (length(intersect(c("query_relation", "input_analyte", "input_common_name", "rxn_partner_common_name", "rxn_partner_ids", "Source" ),
                       names(catalyzeDf)))!=6) {
    stop("Please make sure that the input is the resulting data.frame returned by the rampFastCata() function")
  }

  edges_nodes <- buildCataNetworkDataframe(catalyzeDf)


  #ledges <- data.frame(color = unique(myedges$color.color),
  #           label = c("From HMDB", "From Rhea", "From Both"))
  lnodes <- data.frame(icon.color = c("black", "#d2e5f6"),
                       id =  c("Input", "Output"),
                       label = c("Input Analyte", "Output Analyte"),
                       shape = "icon",
                       icon.face = "'Font Awesome 5 Free'",
                       icon.code = "f111")
  lnodes <- rbind(lnodes, data.frame(icon.color = c("#ca1f7b", "#cc5500", "#008080"),
                                     id = c("From HMDB", "From Rhea", "From Both"),
                                     label = c("From HMDB", "From Rhea", "From Both"),
                                     shape = "icon",
                                     icon.face = "'Font Awesome 5 Free'",
                                     icon.code = "f068"))
  lnodes <- rbind(lnodes, data.frame(icon.color = c("#F9e4bc", "#F9e4bc"),
                                     id = c("Metabolite", "Gene/Protein"),
                                     label = c("Metabolite", "Gene/Protein"),
                                     shape = "icon",
                                     icon.face = "'Font Awesome 5 Free'",
                                     icon.code = c("f111", "f0c8")))

  # Now plot
  visNetwork::visNetwork(edges_nodes$mynodes, edges_nodes$myedges) %>%
    visNetwork::visInteraction(dragNodes = TRUE,
                               dragView = TRUE,
                               navigationButtons=TRUE,zoomView = TRUE) %>%
    visNetwork::visLayout(randomSeed = 123) %>%
    visNetwork::visNodes(size = 20, font = list(size = 25, background = "white")) %>%
    visNetwork::visEdges(width = 1.5, dashes = TRUE, selectionWidth = 3) %>%
    visNetwork::visLegend(useGroups = F, addNodes = lnodes, main = "Legend") %>%
    visNetwork::visGroups()%>%
    visNetwork::addFontAwesome(version = "5.13.0") %>%
    visNetwork::visPhysics(
      barnesHut = list(
        gravitationalConstant = -100,
        centralGravity = 0,
        springConstant = 0
      ),
      stabilization = TRUE) %>%   # explicit edge options
    visNetwork::visOptions(highlightNearest = list(enabled = TRUE, degree = 1,
                                                   labelOnly = FALSE, hover = TRUE))
}

#' Cluster and plot significant pathways by FDR-adjusted pVal
#' @param pathwaysSig output of FilterFishersResults
#' @param pVal Which p value to plot, choose from Raw, FDR or Holm-adjusted
#' @param percAnalyteOverlap Minimum overlap for pathways to be considered similar
#' (Default = 0.2)
#' @param percPathwayOverlap Minimum overlap for clusters to merge (Default = 0.2)
#' @param minPathwayToCluster Minimum number of 'similar' pathways required to start
#' a cluster (medoid) (Default = 3)
#' @param textSize Scales all text in figure (Default=16)
#' @param sigCutoff Aesthetic, shows pvalue cutoff for significant pathways
#' @param interactive If TRUE, return interactive plotly object instead of ggplot object
#' @param db a RaMP database object
#' @examples
#' \dontrun{
#' pathwayResultsPlot(pathwaysSig = filtered.fisher.results, textSize = 8, percAnalyteOverlap = 0.2,
#'    minPathwayToCluster = 2, percPathwayOverlap = 0.2, interactive = FALSE, db = rampDB )
#' }
#' @export
pathwayResultsPlot <- function(pathwaysSig, pVal = "FDR", percAnalyteOverlap = 0.5,
                                 percPathwayOverlap = 0.5, minPathwayToCluster = 3,
                               textSize = 8, sigCutoff = 0.05, interactive=FALSE,
                               db = RaMP()) {


  if( !('cluster_assignment' %in% colnames(pathwaysSig$fishresult))) {
    fishClustering <- findCluster(db = db, fishersDf = pathwaysSig,
                                  percAnalyteOverlap = percAnalyteOverlap,
                                  percPathwayOverlap = percPathwayOverlap,
                                  minPathwayToCluster = minPathwayToCluster
    )
    # assign this here if clustering is preformed here...
    fishresult <- fishClustering$fishresults
  } else {
    message("The input pathway result has already been clustered. Defaulting to existing clustering.")
    fishresult <- pathwaysSig$fishresults
  }

  if (pathwaysSig$analyteType == "genes" | pathwaysSig$analyteType == "metabolites") {
    inPath <- fishresult$Num_In_Path
    totPath <- fishresult$Total_In_Path
  } else {
    inPath <- apply(fishresult, 1, function(x) {
      if (is.na(x["Num_In_Path_Metab"])) {
        return(as.numeric(x["Num_In_Path_Gene"]))
      } else if (is.na(x["Num_In_Path_Gene"])) {
        return(as.numeric(x["Num_In_Path_Metab"]))
      } else {
        return(as.numeric(x["Num_In_Path_Metab"]) + as.numeric(x["Num_In_Path_Gene"]))
      }
    })
    totPath <- apply(fishresult, 1, function(x) {
      if (is.na(x["Total_In_Path_Metab"])) {
        return(as.numeric(x["Total_In_Path_Gene"]))
      } else if (is.na(x["Total_In_Path_Gene"])) {
        return(as.numeric(x["Total_In_Path_Metab"]))
      } else {
        return(as.numeric(x["Total_In_Path_Metab"]) + as.numeric(x["Total_In_Path_Gene"]))
      }
    })
  }

  if (pVal == "FDR") {
    clusterDF <- data.frame(
      x = -log10(fishresult[, grepl("FDR", colnames(fishresult))]),
      y = fishresult$pathwayName,
      inPath = inPath,
      totPath = totPath,
      cluster = fishresult$cluster_assignment,
      pathwaysource = fishresult$pathwaySource,
      analytes = fishresult$analytes
    )
    ylab <- "-log10(FDR pVal)"
  } else if (pVal == "Holm") {
    clusterDF <- data.frame(
      x = -log10(fishresult[, grepl("Holm", colnames(fishresult))]),
      y = fishresult$pathwayName, inPath <- fishresult$Num_In_Path,
      totPath <- fishresult$Total_In_Path,
      cluster = fishresult$cluster_assignment,
      pathwaysource = fishresult$pathwaysource,
      analytes = fishresult$analytes
    )
    ylab <- "-log10(Holm pVal)"
  } else if (pVal == "Raw") {
    clusterDF <- data.frame(
      x = -log10(fishresult[
        ,
        grepl(
          "^(?=.*Pval)(?!.*FDR)(?!.*Holm)(?=.*Combined).*",
          colnames(fishresult)
        )
      ]),
      y = fishresult$pathwayName, inPath <- fishresult$Num_In_Path,
      totPath <- fishresult$Total_In_Path,
      cluster = fishresult$cluster_assignment,
      pathwaysource = fishresult$pathwaysource,
      analytes = fishresult$analytes
    )
    ylab <- "-log10(pVal)"
  } else {
    print("Invalid p value selection, choose from 'Raw', 'FDR' or 'Holm'")
    stop()
  }

  clusterDF <- clusterDF[order(clusterDF$y, decreasing = TRUE), ]
  ## browser()
  clusterDF <- clusterDF %>% tidyr::separate_rows("cluster", sep = ", ")

  clusterDF$cluster <- sapply(clusterDF$cluster, function(x) {
    ifelse(x == "Did not cluster", return(x), return(paste0("Cluster ", x)))
  })

  clusterDF$pathway.db <- with(clusterDF, {
    paste0(y," (", pathwaysource,")")
  })
  p <- clusterDF %>%
    dplyr::mutate("pathway.db" =
                    with(clusterDF,{tidytext::reorder_within(pathway.db,
                                                             x,
                                                             cluster)})) %>%
    ggplot2::ggplot(
      ggplot2::aes_string(y = "x", x = "pathway.db")
    ) +
    ggplot2::geom_segment(ggplot2::aes_string(xend = "pathway.db", y = 0, yend = "x")) +
    suppressWarnings(ggplot2::geom_point(
      stat = "identity",
      ggplot2::aes_string(
        colour = "pathwaysource",
        size = "inPath",
        text = "analytes"
      )
    )) +
    ggplot2::geom_hline(yintercept = -log10(sigCutoff), linetype = "dotted") +
    ggplot2::theme_bw(base_size = textSize) +
    ggplot2::coord_flip() +
    tidytext::scale_x_reordered() +
    ggplot2::labs(x = "", y = ylab) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(colour = "black"),
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      strip.text.y = ggplot2::element_text(angle = 0),
      axis.text = ggplot2::element_text(face = "bold")
    ) +
    with(clusterDF, {
      ggplot2::facet_grid(cluster ~ ., space = "free", scales = "free")
    }) +
    ggplot2::guides(colour = ggplot2::guide_legend(
      override.aes =
        list(size = 10),
      title = "Source Database"
    )) +
    ggplot2::scale_size_area(
      breaks = c(2, 4, 6, 8, 10),
      name = "# of Altered Analytes in Pathway"
    )

  if(interactive){
    return(plotly::ggplotly(p, tooltip="text"))
  }else if (!interactive){
    return(p)
  }else{
    stop("'interactive' must be a boolean")
  }
}

#' Plots an interactive sunburst plot of reaction class
#'
#' @param reactionClassesResults output of getReactionClassesForAnalytes()
#' @return  An interactive HTML sunburst plot that allows the user to pan/zoom into reaction classes of interest.
#' @examples
#' \dontrun{
#' analytes.of.interest = c('chebi:58115', 'chebi:456215', 'chebi:58245', 'chebi:58450',
#' 'chebi:17596', 'chebi:16335', 'chebi:16750', 'chebi:172878',
#' 'chebi:62286', 'chebi:77897', 'uniprot:P30566','uniprot:P30520',
#' 'uniprot:P00568', 'uniprot:P23109', 'uniprot:P22102', 'uniprot:P15531')
#'
#' reaction.classes <- getReactionClassesForAnalytes( analytes = analytes.of.interest, db = RaMP())
#'
#' plotReactionClasses(reaction.classes)
#' }
#'
#' @export

plotReactionClasses <- function(reactionClassesResults) {

  if (missing(reactionClassesResults)) {
    stop("Input is missing. Please input the resulting list of dataframes returned by the getReactionClassesForAnalytes() function")
  }

  if(sum(reactionClassesResults$class_ec_level_1$reactionCount) == 0) {
    message("The input dataframe has no reaction results. plotReactionClasses function is returning without generating a plot.")
    return()
  }

  if (length(intersect(c("class_ec_level_1","class_ec_level_2", "class_ec_level_3", "class_ec_level_4"),names(reactionClassesResults)))!=4) {
    stop("Please make sure that the input is the resulting list of dataframes returned by the getReactionClassesForAnalytes() function")
  }



  sunburst_ontology_reactionclass <- buildReactionClassesSunburstDataframe(reactionClassesResults = reactionClassesResults)

  fig <- plotly::plot_ly(
    color = I("black"),
    marker = list(colors = ~ sunburst_ontology_reactionclass$color)
  )
  fig <- fig %>%
    plotly::add_trace(
      ids = sunburst_ontology_reactionclass$ids,
      labels = sunburst_ontology_reactionclass$labels,
      parents = sunburst_ontology_reactionclass$parents,
      hovertemplate = sunburst_ontology_reactionclass$hovertemplate,
      type = 'sunburst',
      maxdepth = 2,
      domain = list(column = 1),
      name = ""
    )
  fig <- fig %>%
    plotly::layout(
      margin = list(
        l = 0,
        r = 0,
        b = 0,
        t = 0
      )
      # ,marker = list(colors = list(sunburst_ontology_reactionclass$color))
    )

  return(fig)

}

#' Plots an interactive upset plot of overlapping input compounds at reaction class level 1
#'
#' @param reactionsResults output of getReactionsForAnalytes()
#' @param includeCofactorMets include metabolites labeled at cofactors within ChEBI (Default = FALSE)
#' @return  An interactive HTML upset plot that allows the user to visualize the overlap in the number of
#' input compounds across level 1 of reaction classes.
#' @examples
#' \dontrun{
#' analytes.of.interest = c('chebi:58115', 'chebi:456215', 'chebi:58245', 'chebi:58450',
#' 'chebi:17596', 'chebi:16335', 'chebi:16750', 'chebi:172878',
#' 'chebi:62286', 'chebi:77897', 'uniprot:P30566','uniprot:P30520',
#' 'uniprot:P00568', 'uniprot:P23109', 'uniprot:P22102', 'uniprot:P15531')
#'
#' reactionsLists <- RaMP::getReactionsForAnalytes(analytes = analytes.of.interest,
#'     includeTransportRxns = F, humanProtein = T, db = rampDB)
#'
#' plotAnalyteOverlapPerRxnLevel(reactionsResults = reactionsLists)
#' }
#'
#' @export
plotAnalyteOverlapPerRxnLevel <- function(reactionsResults, includeCofactorMets = FALSE) {

  if (missing(reactionsResults)) {
    stop("Input is missing. Please input the resulting list of dataframes returned by the getReactionsForAnalytes() function")
  }

  if(nrow(reactionsResults$met2rxn) ==0 && nrow(reactionsResults$prot2rxn) == 0) {
    message("The input has no reaction results. plotAnalyteOverlapPerRxnLevel function is returning without generating a plot.")
    return()
  }

  if (length(intersect(c("met2rxn","prot2rxn", "metProteinCommonReactions"),names(reactionsResults)))!=3) {
    stop("Please make sure that the input is the resulting list of dataframes returned by the getReactionsForAnalytes() function")
  }

  input2reactions_list <- buildAnalyteOverlapPerRxnLevelUpsetDataframe(reactionsResults = reactionsResults, includeCofactorMets = includeCofactorMets)

  fig <- upsetjs::upsetjs() %>%
    upsetjs::fromList(input2reactions_list) %>%
    upsetjs::generateDistinctIntersections() %>% upsetjs::interactiveChart() %>%
    upsetjs::chartLabels(set.name = "Number of Reactions") %>%
    upsetjs::chartLayout(width.ratios = c(0.15,0.2,0.65),
                         height.ratios = c(0.6,.4)) %>%
    upsetjs::chartFontSizes(set.label = "13px", chart.label = "14px")

  return(fig)

}

#' Cluster and plot significant ontologies by FDR-adjusted pVal
#' @param ontologiesSig output of FilterFishersResults
#' @param pVal Which p value to plot, choose from Raw, FDR or Holm-adjusted
#' @param textSize Scales all text in figure (Default=16)
#' @param sigCutoff Aesthetic, shows pvalue cutoff for significant ontologies
#' @param interactive If TRUE, return interactive plotly object instead of ggplot object
#' @param db a RaMP database object
#' @export
ontologyEnrichmentResultsPlot <- function(ontologiesSig, pVal = "FDR",
                                          textSize = 8,
                                          sigCutoff = 0.05, interactive=FALSE,
                                          db = RaMP()) {

  inOntology <- ontologiesSig$Num_In_Ontology
  totOntology <- ontologiesSig$Total_In_Ontology


  if (pVal == "FDR") {
    plotDF <- data.frame(
      x = -log10(ontologiesSig[, grepl("FDR", colnames(ontologiesSig))]),
      y = ontologiesSig$Ontology,
      inOntology = inOntology,
      totOntology = totOntology,
      ontologytype = ontologiesSig$HMDBOntologyType
    )
    ylab <- "-log10(FDR pVal)"
  } else if (pVal == "Holm") {
    plotDF <- data.frame(
      x = -log10(ontologiesSig[, grepl("Holm", colnames(ontologiesSig))]),
      y = ontologiesSig$Ontology, inOntology <- ontologiesSig$Num_In_Ontology,
      totOntology <- ontologiesSig$Total_In_Ontology,
      ontologytype = ontologiesSig$HMDBOntologyType
    )
    ylab <- "-log10(Holm pVal)"
  } else if (pVal == "Raw") {
    plotDF <- data.frame(
      x = -log10(ontologiesSig[
        ,
        grepl(
          "^(?=.*Pval)(?!.*FDR)(?!.*Holm)(?=.*Combined).*",
          colnames(ontologiesSig)
        )
      ]),
      y = ontologiesSig$Ontology, inOntology <- ontologiesSig$Num_In_Ontology,
      totOntology <- ontologiesSig$Total_In_Ontology,
      ontologytype = ontologiesSig$HMDBOntologyType
    )
    ylab <- "-log10(pVal)"
  } else {
    print("Invalid p value selection, choose from 'Raw', 'FDR' or 'Holm'")
    stop()
  }

  ## plotDF <- plotDF[order(plotDF$y, decreasing = TRUE), ]
  plotDF$y <- factor(plotDF$y,
                     levels  = unique(plotDF$y[order(plotDF$x, decreasing = FALSE)]))

  p <- plotDF %>%
    ggplot2::ggplot(
      ggplot2::aes_string(y = "x", x = "y")
    ) +
    ggplot2::geom_segment(ggplot2::aes_string(xend = "y", y = 0, yend = "x")) +
    suppressWarnings(ggplot2::geom_point(
      stat = "identity",
      ggplot2::aes_string(
        colour = "ontologytype",
        size = "inOntology",
      )
    )) +
    ggplot2::geom_hline(yintercept = -log10(sigCutoff), linetype = "dotted") +
    ggplot2::theme_bw(base_size = textSize) +
    ggplot2::coord_flip() +
    tidytext::scale_x_reordered() +
    ggplot2::labs(x = "", y = ylab) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(colour = "black"),
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      strip.text.y = ggplot2::element_text(angle = 0),
      axis.text = ggplot2::element_text(face = "bold")
    ) +
    with(plotDF, {
      ggplot2::facet_grid(ontologytype ~ ., space = "free", scales = "free")
    }) +
    ggplot2::guides(colour = ggplot2::guide_legend(
      override.aes =
        list(size = 10),
      title = "Ontology Type"
    )) +
    ggplot2::scale_size_area(
      breaks = c(2, 4, 6, 8, 10),
      name = "# of Altered Metabolites in Ontology"
    )

  if(interactive){
    return(plotly::ggplotly(p, tooltip="text"))
  }else if (!interactive){
    return(p)
  }else{
    stop("'interactive' must be a boolean")
  }
}


#' Plots an interactive sunburst plot of chemical classes
#'
#' @param chemicalClassResults output of getChemClass()
#' @param plotType choice of 'sunburst' or 'treemap' plot type (default = 'sunburst')
#' @return  An interactive HTML sunburst plot that allows the user to pan/zoom into reaction classes of interest.
#' @export

plotChemicalClass <- function(chemicalClassResults, plotType = "sunburst") {

  if (missing(chemicalClassResults)) {
    stop("Input is missing. Please input the resulting list of dataframes returned by the chemicalClassSurvey() function")
  }

  if(nrow(chemicalClassResults$met_classes) == 0) {
    message("The input dataframe has no metabolites. plotChemicalClass function is returning without generating a plot.")
    return()
  }

  if (length(intersect(c("count_summary","met_classes", "query_report"),names(chemicalClassResults)))!=3) {
    stop("Please make sure that the input is the resulting list of dataframes returned by the chemicalClassSurvey() function")
  }



  sunburst_ontology_chemicalClass <- buildChemicalClassDataframe(chemicalClassResults = chemicalClassResults)

  fig <- plotly::plot_ly(
    color = I("black")
  )
  fig <- fig %>%
    plotly::add_trace(
      ids = sunburst_ontology_chemicalClass$ids,
      labels = sunburst_ontology_chemicalClass$labels,
      parents = sunburst_ontology_chemicalClass$parents,
      type = plotType,
      maxdepth = 2,
      domain = list(column = 1)
    )
  fig <- fig %>%
    plotly::layout(
      margin = list(l = 0, r = 0, b = 0, t = 0),
      colorway = pals::brewer.divdiv()
    )

  return(fig)

}
