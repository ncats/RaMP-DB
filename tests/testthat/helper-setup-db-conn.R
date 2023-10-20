library(properties)

# get db props
dbpass <- properties::read.properties('../../dbprops.txt')
rampDB <- NULL

if(grepl("SQLite", dbpass$dbname)) {
  rampDB <- RaMP()
} else {
  print("MySQL")
  rampDB <- RaMP:::.RaMP(
    driver = RMariaDB::MariaDB(), dbname = dbpass$dbname, username = dbpass$username,
    conpass = dbpass$conpass, host = dbpass$hostname, port = as.integer(dbpass$port))
}

# set dbf database file for sqlite-specific tests
sqliteRampDB <- RaMP()
dbf <- sqliteRampDB@dbname
