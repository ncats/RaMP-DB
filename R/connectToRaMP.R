#' Connect to RaMP database (requires mysql password when running locally)
#' 
#' @param dbname the name of the database (by default a string that is database name including all tables
#' @param username a string that is username for database (default: root)
#' @param conpass password for database (string)
#' @param host a string that stand for host 
#' @examples
#' \dontrun{
#' con <- connectToRaMP(dbname="ramp",username="root",conpass="mypassword")
#' }
#' @export
connectToRaMP <- function(dbname = "ramp",
                       username = "root",
                       conpass = NULL,
                       host = NULL){
#                       update = FALSE){
  
  if(is.null(conpass)) {
	stop("Please define the password for the mysql connection")
  }

  con <- DBI::dbConnect(
    drv = RMySQL::MySQL(),
    dbname = dbname,
    username = username,
    password = conpass
  )
 return(con) 
}

