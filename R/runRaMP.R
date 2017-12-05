#' Run Shiny App
#' 
#' This function connects to the RaMP database (provided within the package), and launches the RaMP web application.  The password for the mysql connection must be provided.
#' 
#' @param dbname the name of the dadtabase (by default a string that is database name including all tables
#' @param username a string that is username for database
#' @param password a string this is password for database
#' @param host a string that stand for host 
#' @param update a logic value determines if RaMP update local data
#' @export
runRaMPapp <- function(dbname = "ramp",
                       username = "root",
                       password = NULL,
                       host = NULL,
                       update = FALSE){
  
  if(is.null(password)) {
	stop("Please define the password for the mysql connection")
  }

  con <<- DBI::dbConnect(
    drv = RMySQL::MySQL(),
    dbname = dbname,
    username = username,
    password = password
  )
  appDir <- system.file("shinyApp",package = "RaMP")
  if (appDir == ""){
    dbDisconnect(con)
    stop("Could not find example directory. Try re-installing 'ramp'.")
  }
 
  shiny::runApp(appDir,display.mode = "normal")
  
}

