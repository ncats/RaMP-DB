##' Generic function for removing RaMP identifiers and making query results more pleasant to look at
##' @param data A dataframe from the following functions: runCominedFishersTests, getAnalyteFromPathway, getPathwayFromAnalyte, chemicalClassSurvey, and getChemicalProperties
##' @param truncate_cells_at String length to truncate cells at (No truncation if NA)
##' @param show_n_rows Only print first n rows of input (Shows all if NA)
##' @return Printed out and prettified version of input
##' @author Andrew Christopher Patt
##' @export
display_results <- function(data, show_n_rows = 6, truncate_cells_at = 20) {
  if (class(data) != "data.frame" & (class(data) != "list" & length(data) != 1)) {
    stop("Input should be a dataframe resulting from runCombinedFishersTest, getAnalyteFromPathway, getPathwayFromAnalyte, chemicalClassSurvey, or getChemicalProperties")
  }
  if (class(data) == "list") {
    data <- data[[1]]
  }
  rownames(data) <- NULL
  data <- data %>%
    dplyr::select(!dplyr::contains("RaMP"))
  if (!is.na(show_n_rows)) {
      norows <- nrow(data)
      if(norows > show_n_rows){
          data <- data[1:show_n_rows, ]
      }
  }
  if (!is.na(truncate_cells_at)) {
    data <- apply(data, 2, function(x) {
        if (class(x) == "character") {
            x <- sapply(x, function(y){
                if(is.na(y)){
                    return(y)
                }else if (nchar(y) > truncate_cells_at) {
                    return(paste0(substring(y, 1,
                                            truncate_cells_at), "..."))
                } else {
                    return(y)
                }
            })
        } else {
            return(x)
        }
    })
    print(data.frame(data))
    if(!is.na(show_n_rows)){
        print(paste0("Dataframe has ",norows," row(s)"))
    }
  } else {
      print(data.frame(data))
      if(!is.na(show_n_rows)){
        print(paste0("Dataframe has ",norows," row(s)"))
    }
  }
}
