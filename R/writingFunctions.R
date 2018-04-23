#' The function write a data.frame to a csv files.
#' 
#' @param df a data frame returned by functions that queires database.
#' @param f.name a string that represents output file name.
#' @export
write_to_csv <- function(df,f.name){
  write.csv(df,row.names = F)
}