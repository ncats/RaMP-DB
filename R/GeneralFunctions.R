#' Load pathway overlap matrices for find_clusters function
#'
#' @return A list of pathway overlap matrices for clustering
loadOverlapMatrices<- function() {
  gene_result <- metabolite_result <- analyte_result <- c()
  load(system.file(package = "RaMP",... = "extdata/genes_overlap_matrix.RData"))
  load(system.file(package = "RaMP",... = "extdata/metabolites_overlap_matrix.RData"))
  load(system.file(package = "RaMP",... = "extdata/analytes_overlap_matrix.RData"))
  return(list(gene=gene_result,metab=metabolite_result,
	analyte=analyte_result))
}
