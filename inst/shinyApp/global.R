# The global variables that are important to keywords search and fisher test.

# read db connection from db.properties file and store them in GlobalEnv
dbProp <- properties::read.properties('db.properties')

.GlobalEnv$.host <- dbProp$host
.GlobalEnv$.dbname <- dbProp$dbname
.GlobalEnv$.username <- dbProp$username
.GlobalEnv$.conpass <- dbProp$conpass

FisherTestData <- list(
  metabolites = read.csv("FisherTestDataMetabolites.csv"),
  genes = read.csv("FisherTestDataGenes.csv")
)

# Unlist automatically select first column, the first column should
# be the data containing keywords
con <- DBI::dbConnect(RMySQL::MySQL(),
                      username = .username,
                      dbname = .dbname,
                      password = .conpass,
                      host = .host)
# kw_biofluid <- unique(unlist(read.csv("biofluid.csv",header = F,stringsAsFactors = F),use.names = F))
# kw_pathway <- unique(unlist(read.csv("pathway.csv",header = F,stringsAsFactors = F),use.names = F))
# kw_analyte <- unique(unlist(read.csv("analyte.csv",header = F,stringsAsFactors = F),use.names = F))
# kw_source <- unique(unlist(read.csv("source.csv",header = F,stringsAsFactors = F),use.names = F))
kw_biofluid <- unname(unlist(unique(as.vector(DBI::dbGetQuery(con,'select commonName from ontology;')))))
kw_pathway <- unname(unlist(unique(as.vector(DBI::dbGetQuery(con,'select pathwayName from pathway;')))))
kw_analyte <- unname(unlist(unique(as.vector(DBI::dbGetQuery(con,'select Synonym from analytesynonym;')))))
kw_source <- unname(unlist(unique(as.vector(DBI::dbGetQuery(con,'select sourceId from source;')))))

load("/inst/extdata/FT_data.Rdata")

DBI::dbDisconnect(con)
