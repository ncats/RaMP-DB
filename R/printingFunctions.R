##' Generic function for removing RaMP identifiers and making query results more pleasant to look at
##' @param data A dataframe from the following functions: runCominedFishersTests, getAnalyteFromPathway, getPathwayFromAnalyte, chemicalClassSurvey, and getChemicalProperties
##' @param showNRows Only print first n rows of input (Shows all if NA)
##' @return Printed out and prettified version of input
##' @author Andrew Christopher Patt
## @import is from methods
#' @noRd
cleanup<- function(data, showNRows = 6) {
  if (!inherits(data,"data.frame") & (!inherits(data,"list") & length(data) != 1)) {
    stop("Input should be a dataframe resulting from runCombinedFishersTest, getAnalyteFromPathway, getPathwayFromAnalyte, chemicalClassSurvey, or getChemicalProperties")
  }
#  if (class(data) == "list") {
   if( is( data, "list")){
    data <- data[[1]]
  }
  rownames(data) <- NULL
  data <- data %>%
    dplyr::select(!dplyr::contains("RaMP"))
  ## if (!is.na(showNRows)) {
  ##   norows <- nrow(data)
  ##   if (norows > showNRows) {
  ##     data <- data[1:showNRows, ]
  ##   }
  ## }
  ## if (!is.na(truncate_cells_at)) {
  ##   data <- apply(data, 2, function(x) {
  ##     if (class(x) == "character") {
  ##       x <- sapply(x, function(y) {
  ##         if (is.na(y)) {
  ##           return(y)
  ##         } else if (nchar(y) > truncate_cells_at) {
  ##           return(paste0(substring(
  ##             y, 1,
  ##             truncate_cells_at
  ##           ), "..."))
  ##         } else {
  ##           return(y)
  ##         }
  ##       })
  ##     } else {
  ##       return(x)
  ##     }
  ##   })
  ##   return(data.frame(data))
  ##   ## if(!is.na(showNRows)){
  ##   ##     print(paste0("Dataframe has ",norows," row(s)"))
  ##   ## }
  ## } else {
    return(data.frame(data))
    ##   if(!is.na(showNRows)){
    ##     print(paste0("Dataframe has ",norows," row(s)"))
    ## }
  ## }
}

