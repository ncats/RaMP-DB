#' Function that writes output from getPathwayFromAnalyte() to a CSV file
#'
#' @param myPathways data frame returned by function getPathwayFromAnalyte()
#' @param outputFile name of output file
#' @param ... Internal Use - for handling deprecated parameter names
#' @export
writePathwaysToCSV <- function(myPathways = "none", outputFile = "none", ...) {
  myPathways <- handleRenamedParameter(argument = myPathways, oldName = "mypathways", version = "3.0")
  outputFile <- handleRenamedParameter(argument = outputFile, oldName = "outputfile", version = "3.0")

  if(length(myPathways) == 1){
    if (myPathways == "") {
      stop("Be sure to specify the output of the function getPathwayFromAnalyte() and an output file")
    }}
  if(length(outputFile) == 1){
    if (outputFile == "") {
      stop("Be sure to specify the output of the function getPathwayFromAnalyte() and an output file")
    }}
  if (!all(c(
    "pathwayName", "pathwaySource",
    "pathwayId", "inputId", "commonName"
  ) %in% colnames(myPathways))) {
    stop("Make sure that your input data is the output of the function getPathwayFromAnalyte()")
  }
  out <- myPathways[, c(
    "pathwayName", "pathwaySource",
    "pathwayId", "inputId", "commonName"
  )]
  utils::write.csv(out, file = outputFile, row.names = F, quote=TRUE)
}


#' Function that writes Fishers Test results, after clustering to a CSV file
#'
#' @param fishResults a data frame returned by function runCombinedFisherTest()
#' @param outputFile name of output file
#' @param includeRaMPids whether or not to include rampIds (default is FALSE)
#' @export
writeFishersResults <- function(fishResults = "none", outputFile = "none", includeRaMPids = FALSE) {
  if(length(fishResults) == 1){
    if (fishResults == "") {
      stop("Be sure to specify the output of the function findCluster()")
    }}
  clusters <- fishResults$cluster_list
  if (is.null(clusters)) {
    out <- fishResults$fishresults
    mycols <- setdiff(colnames(out), c("pathwayRampId", "pathwayName"))
    mycols <- c("pathwayName", mycols)
    utils::write.csv(out[, mycols], file = outputFile, row.names = F)
  } else {
    cluster_list <- fishResults$cluster_list
    out <- fishResults
    rampOut <- out$fishresults
    if (!is.null(rampOut)) {
      if (out$analyte_type == "both") {
        rampOut <- rampOut[, c(
          "pathwayName", "Pval.Metab", "Num_In_Path.Metab", "Total_In_Path.Metab",
          "Pval.Gene", "Num_In_Path.Gene", "Total_In_Path.Gene", "Pval_combined",
          "Pval_combined_FDR", "Pval_combined_Holm", "pathwaysourceId", "pathwaysource",
          "cluster_assignment", "rampIds"
        )]

        colnames(rampOut) <- c(
          "Pathway Name", "Raw Fisher's P Value (Metabolites)", "User Metabolites in Pathway",
          "Total Metabolites in Pathway", "Raw Fisher's P Value (Genes)", "User Genes in Pathway",
          "Total Genes in Pathway", "Raw Fisher's P Value (Combined)", "FDR Adjusted P Value (Combined)",
          "Holm Adjusted P Value (Combined)", "Source ID", "Source DB", "Cluster", "rampIds"
        )
        rampOut <- rampOut[order(rampOut[, "Holm Adjusted P Value (Combined)"]), ]
      } else {
        results_fisher <- rampOut[, c(
          "pathwayName", "Pval", "Pval_FDR", "Pval_Holm", "pathwaysourceId", "pathwaysource",
          "Num_In_Path", "Total_In_Path", "cluster_assignment", "rampIds"
        )]
        colnames(rampOut) <- c(
          "Pathway Name", "Raw Fisher's P Value", "FDR Adjusted P Value", "Holm Adjusted P Value",
          "Source ID", "Source DB", "User Analytes in Pathway", "Total Analytes in Pathway",
          "Cluster", "rampIds"
        )
        rampOut <- rampOut[order(rampOut[, "Holm Adjusted P Value"]), ]
      }
      # 	      utils::write.csv(rampOut,outputFile,row.names = FALSE)
    } else {
      rampOut <- "No significant results"
      # 	utils::write.csv(c("No significant results"),outputFile,row.names = FALSE)
    }
  }

  if (!includeRaMPids) {
    rampOut <- rampOut[, -ncol(rampOut)]
  }
  if (outputFile=="none") {
    return(rampOut[order(rampOut[, "Cluster"]), ])
  } else {
    utils::write.csv(rampOut[order(rampOut[, "Cluster"]), ], outputFile, row.names = FALSE)
  }
}
