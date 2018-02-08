# The global variables that are important to keywords search and fisher test.

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
kw_biofluid <- unname(unlist(unique(as.vector(dbGetQuery(con,'select commonName from ontology;')))))
kw_pathway <- unname(unlist(unique(as.vector(dbGetQuery(con,'select pathwayName from pathway;')))))
kw_analyte <- unname(unlist(unique(as.vector(dbGetQuery(con,'select Synonym from analytesynonym;')))))
kw_source <- unname(unlist(unique(as.vector(dbGetQuery(con,'select sourceId from source;')))))

DBI::dbDisconnect(con)
