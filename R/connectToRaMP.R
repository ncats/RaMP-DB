#' Set Connection Parameters or RaMP
#'
#' @param dbname the name of the database (by default a string that is database name including all tables
#' @param username a string that is username for database (default: root)
#' @param conpass password for database (string)
#' @param host a string that stand for host
#' @param socket optional, location of mySQL.sock file (useful when running RaMP on remote clusters)
#'
#' @examples
#' \dontrun{
#' pkg.globals <- setConnectionToRaMP(dbname="ramp2",username="root",conpass="",host = "localhost")
#' }
#' @export
setConnectionToRaMP <- function(dbname="ramp",
                       username = "root",
                       conpass = "",
                       host ="localhost",
		       socket = "", is_sqlite = F, sqlite_file_path = ""){
  pkg.globals <- new.env()
  pkg.globals$dbname=dbname
  pkg.globals$username=username
  pkg.globals$conpass=conpass
  pkg.globals$host=host
  if(socket ==""){
      pkg.globals$socket=NULL
  }else{
      pkg.globals$socket = socket
  }

  if(is_sqlite) {
    pkg.globals$is_sqlite = T
    pkg.globals$sqlite_file_path = sqlite_file_path
  }

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
  if(!exists("pkg.globals")) {
    stop("Be sure the run the setConnectionToRaMP() and assign it to pkg.globals");
  }

  if(get("is_sqlite",pkg.globals)) {
    db = RSQLite::SQLite()
    con = RSQLite::dbConnect(db, get("sqlite_file_path", pkg.globals))
  } else {
    if(!is.null(get("socket",pkg.globals))) {
      con <- RMariaDB::dbConnect(
        drv = RMariaDB::MariaDB(),
        dbname = get("dbname",pkg.globals),
        username = get("username",pkg.globals),
        password = get("conpass",pkg.globals),
        host = get("host",pkg.globals),
        unix.socket = get("socket",pkg.globals)
      )
    } else {
      con <- RMariaDB::dbConnect(
        drv = RMariaDB::MariaDB(),
        dbname = get("dbname",pkg.globals),
        username = get("username",pkg.globals),
        password = get("conpass",pkg.globals),
        host = get("host",pkg.globals)
      )
    }
  }
  return(con)
}


