#' Plots a network based on gene-metabolite relationships
#' @importFrom magrittr %>%
#'
#' @param catalyzedf a data.frame output by rampFastCata() that contains analytes that are in the same reaction
#' @return  An interactive HTML plot that allows the user to pan/zoom into regions of interest. User genes/metabolites are highlighted in blue, whereas analytes found by the function are orange.
#' @examples
#' \dontrun{
#' pkg.globals <- setConnectionToRaMP(dbname="ramp2",username="root",conpass="",host = "localhost")
#' catalyzedf <- rampFastCata(analytes="creatine",NamesOrIds="names")
#' plotCataNetwork(catalyzedf)
#' }
#' @export
plotCataNetwork <- function(catalyzedf = "") {

        if(nrow(catalyzedf) == 0) {
          message("The input dataframe has 0 rows. plotCataNetwork function is returning without generating a plot.")
          return()
        }

        if (length(intersect(c("input_analyte","rxn_partner_common_name", "rxn_partner_ids"),colnames(catalyzedf)))!=3) {
                stop("Please make sure that the input is the resulting data.frame returned by the rampFastCata() function")
        }

        toplot = catalyzedf[,c("input_analyte","rxn_partner_common_name")]
        colnames(toplot)<-c("from","to")

        #colopts <- brewer.pal(12,"Set3")
        mycol=rep("black")
        mysize<-rep(8,nrow(toplot))
        mynames<-rep(NA,nrow(toplot))

        myedges <- cbind(toplot,mycol,mysize,mynames)

        # Now set nodes
        mynodes=c(unique(toplot$from),unique(toplot$to))
        mycol=c(rep("blue",length(unique(toplot$from))),
                rep("orange",length(unique(toplot$to))))
        mysize <- rep(8,length(mynodes))
        mynames <- mynodes

        mynodes=data.frame(color=mycol,size=mysize,id=mynames,label=mynames)

        # Now plot
        visNetwork::visNetwork(mynodes, myedges, width = "100%",height="1000px") %>%
		visNetwork::visInteraction(dragNodes = FALSE,
                 dragView = TRUE,hideEdgesOnDrag=TRUE,hideNodesOnDrag=TRUE,
                 navigationButtons=TRUE,zoomView = TRUE) %>%
  		visNetwork::visLayout(randomSeed = 123) %>%
		visNetwork::visPhysics(
		  barnesHut = list(
		    gravitationalConstant = -100,
		    centralGravity = 0,
		    springConstant = 0
		  ),
		  stabilization = TRUE)
        #return(NULL) #return(list(nodes=mynodes,edges=myedges))
}

#' Cluster and plot significant pathways by FDR-adjusted pval
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
pathwayResultsPlot <- function(db = RaMP(), pathwaysSig, pval = "FDR", perc_analyte_overlap = 0.2,
                                 perc_pathway_overlap = 0.2, min_pathway_tocluster = 3,
                                 text_size = 8, sig_cutoff = 0.05, interactive=FALSE) {

  if( !('cluster_assignment' %in% colnames(pathwaysSig$fishresult))) {
    fishClustering <- findCluster(db = db, pathwaysSig,
                                  perc_analyte_overlap = perc_analyte_overlap,
                                  perc_pathway_overlap = perc_pathway_overlap,
                                  min_pathway_tocluster = min_pathway_tocluster
    )
  } else {
    message("The input pathway result has already been clustered. Defaulting to existing clustering.")
    fishresult <- pathwaysSig$fishresults
  }

  fishresult <- fishClustering$fishresults
  if (pathwaysSig$analyte_type == "genes" | pathwaysSig$analyte_type == "metabolites") {
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

  if (pval == "FDR") {
    clusterDF <- data.frame(
      x = -log10(fishresult[, grepl("FDR", colnames(fishresult))]),
      y = fishresult$pathwayName,
      inPath = inPath,
      totPath = totPath,
      cluster = fishresult$cluster_assignment,
      pathwaysource = fishresult$pathwaySource,
      analytes = fishresult$analytes
    )
    ylab <- "-log10(FDR pval)"
  } else if (pval == "Holm") {
    clusterDF <- data.frame(
      x = -log10(fishresult[, grepl("Holm", colnames(fishresult))]),
      y = fishresult$pathwayName, inPath <- fishresult$Num_In_Path,
      totPath <- fishresult$Total_In_Path,
      cluster = fishresult$cluster_assignment,
      pathwaysource = fishresult$pathwaysource,
      analytes = fishresult$analytes
    )
    ylab <- "-log10(Holm pval)"
  } else if (pval == "Raw") {
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
    ylab <- "-log10(pval)"
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
    ggplot2::geom_hline(yintercept = -log10(sig_cutoff), linetype = "dotted") +
    ggplot2::theme_bw(base_size = text_size) +
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
