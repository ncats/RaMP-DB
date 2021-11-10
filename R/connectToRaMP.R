#' Set Connection Parameters or RaMP
#'
#' @param dbname the name of the database (by default a string that is database name including all tables
#' @param username a string that is username for database (default: root)
#' @param conpass password for database (string)
#' @param host a string that stand for host
#'
#' @examples
#' \dontrun{
#' pkg.globals <- setConnectionToRaMP(dbname="ramp2",username="root",conpass="",host = "localhost")
#' }
#' @export
setConnectionToRaMP <- function(dbname="ramp",
                       username = "root",
                       conpass = NULL,
                       host ="localhost"){
  pkg.globals <- new.env()
  pkg.globals$dbname=dbname
  pkg.globals$username=username
  pkg.globals$conpass=conpass
  pkg.globals$host=host
  return(pkg.globals)
}

#' Connect to RaMP database (requires mysql password when running locally)
#'
#' @examples
#' \dontrun{
#' pkg.globals <- setConnectionToRaMP(dbname="ramp2",username="root",conpass="",host = "localhost")
#' con <- connectToRaMP()
#' }
#' @return MySQL connection based on given dbname,username,password, and host
#' @export
connectToRaMP <- function() {
  if(is.null(get("conpass",pkg.globals))) {
        stop("Please define the password for the mysql connection using the setConnectionToRaMP() function")
  }

  con <- DBI::dbConnect(
    drv = RMariaDB::MariaDB(),
    dbname = get("dbname",pkg.globals),
    username = get("username",pkg.globals),
    password = get("conpass",pkg.globals),
    host = get("host",pkg.globals)
  )
 return(con)
}

#' Connect to RaMP database (requires mysql password when running locally)
#' 
#' @param dbname the name of the database (by default a string that is database name including all tables
#' @param username a string that is username for database (default: root)
#' @param conpass password for database (string)
#' @param host a string that stand for host 
#' 
#' @examples
#' \dontrun{
#' con <- connectToRaMP(dbname="ramp",username="root",conpass="mypassword",host = "localhost")
#' }
#' @return MySQL connection based on given dbname,username,password, and host
#' @export
OLDconnectToRaMP <- function(dbname = "ramp",
                       username = "root",
                       conpass = NULL,
                       host ="localhost"){

  
  if(is.null(conpass)) {
	stop("Please define the password for the mysql connection")
  }

  con <- DBI::dbConnect(
    drv = RMariaDB::MariaDB(),
    dbname = dbname,
    username = username,
    password = conpass,
    host = host
  )
 return(con) 
}

