#' Run Shiny App
#' 
#' This function launches the RShiny app.  It requires a connection to the RaMP database as input, which requires running the function connectoToRaMP() and providing the MySQL password.  
#' 
#' @param con a connection object returned from the function connectToRaMP()
#' @examples
#' \dontrun{
#' con <- connectToRaMP(dbname="ramp",username="root",password="mypassword")
#' runRaMPapp(con=con)
#' }
#' @export
runRaMPapp <- function(con = NULL) {
 
  if(is.null(con)) {
	stop("Please connect to the database first using the funciton connectToRaMP")
  }

  appDir <- system.file("shinyApp",package = "RaMP")
  if (appDir == ""){
    DBI::dbDisconnect(con)
    stop("Could not find example directory. Try re-installing 'ramp'.")
  }
 
  shiny::runApp(appDir,display.mode = "normal")
  
}

