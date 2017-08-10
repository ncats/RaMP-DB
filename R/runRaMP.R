#' run shiny app
#' 
#' if has database locally, user must provide database name
#' username password
#' @param dbname a string that is database name including all tables
#' @param username a string that is username for database
#' @param password a string this is password for database
#' @param host a string that stand for host 
#' @param update a logic value determines if RaMP update local data
#' @export
runRaMPapp <- function(dbname = "mathelabramp",
                       username = "root",
                       password = "Ehe131224",
                       host = NULL,
                       update = FALSE){
  if(update){
    # Next version...
  }
  appDir <- system.file("shinyApp",package = "RaMP")
  if (appDir == ""){
    stop("Could not find example directory. Try re-installing 'ramp'.")
  }
  con <<- dbConnect(
    drv = RMySQL::MySQL(),
    dbname = dbname,
    username = username,
    password = password
  )
  shiny::runApp(appDir,display.mode = "normal")
  
}
