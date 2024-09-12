library(properties)

## if (Sys.getenv("MYSQL_TEST") == "true") {
##    dbpass <- properties::read.properties('../../tests/local_mysql.dbprops.txt')
## } else {
##}

rampDB <- NULL

## if(grepl("SQLite", dbpass$dbname)) {
##  print("Testing SQLite DB")
 rampDB <- RaMP()
## } else {
 ## print("Testing MySQL DB")
 ## rampDB <- RaMP:::.RaMP(
 ##   driver = RMariaDB::MariaDB(), dbname = dbpass$dbname, username = dbpass$username,
 ##   conpass = dbpass$conpass, host = dbpass$hostname, port = as.integer(dbpass$port))
#}

sqliteRampDB <- RaMP()
dbf <- sqliteRampDB@dbname
