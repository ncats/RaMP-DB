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

#' Update and save the overlap matrices based on current version of RaMP
#' 
#' @param method a string that specifies algorithm to compute overlap matrix
#' should be 'balanced' or 'weighted'
#' @param all a string that specifies which matrices to compute, should be in
#' 'all','analyte'
#' @param username a string that specifies name of MySQL database
#' @param dbname a string that specifies database name of MySQL database
#' @param conpass a string that specifies password for database connection
#' @param host a string that specifes host for database connection
#' @export
updateOverlapMatrices <- function(method,all,
                                  conpass,
                                  host = 'localhost',dbname = 'ramp',
                                  username = 'root'){
  if(!(method %in% c('balanced','weighted'))){
    stop('Wrong input for argument method')
  }
  if(!(all %in%c('all','metabolite','gene','analyte'))){
    stop('Wrong input for argument all')
  }
  
  if(all == 'all'){
    result <- updateOverlapMatrix(min_analyte = 5,overlapmethod = 'balanced',together = T)
    save(result[[1]],system.file(package = "RaMP",... = "extdata/metabolites_overlap_matrix.RData"))
    save(result[[2]],system.file(package = "RaMP",... = "extdata/genes_overlap_matrix.RData"))
  } else if(all == 'analyte'){
    result <- updateOverlapMatrix(min_analyte = 5,overlapmethod = 'balanced',together = F)
    save(result,system.file(package = "RaMP",... = "extdata/analytes_overlap_matrix.RData"))
   
  }
}
