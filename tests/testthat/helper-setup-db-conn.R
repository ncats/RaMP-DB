#library(properties)
#
#if (Sys.getenv("MYSQL_TEST") == "true") {
#    dbpass <- properties::read.properties('../../tests/local_mysql.dbprops.txt')
#} else {
#    dbpass <- properties::read.properties('../../dbprops.txt')
#}
#
# get db props
#rampDB <- NULL
#
#if(grepl("SQLite", dbpass$dbname)) {
#  print("Testing SQLite DB")
#  v <- gsub("SQLite_v", "", dbpass$dbname)
#  rampDB <- RaMP(v)
#} else {
#  print("Testing MySQL DB")
#  rampDB <- RaMP:::.RaMP(
#    driver = RMariaDB::MariaDB(), dbname = dbpass$dbname, username = dbpass$username,
#    conpass = dbpass$conpass, host = dbpass$hostname, port = as.integer(dbpass$port))
#}

# set dbf database file for sqlite-specific tests
#sqliteRampDB <- RaMP()
#dbf <- sqliteRampDB@dbname
