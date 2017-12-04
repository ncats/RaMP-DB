# con <- dbPool(
#   drv = RMySQL::MySQL(),
#   dbname = "mathelabramp",
#   username = "root",
#   password = "Ramp340!",
#   minSize = 2,
#   maxSize = 100,
#   idleTimeout = 3600000
# )
# con <- dbConnect(
#   drv = RMySQL::MySQL(),
#   dbname = "mathelabramp",
#   username = "root",
#   password = "Ehe131224"
# )

FisherPathwayTable <- read.csv("dataFisherTest.csv")

kw_biofluid <- unique(unlist(read.csv("biofluid.csv",header = F,stringsAsFactors = F),use.names = F))
kw_pathway <- unique(unlist(read.csv("pathway.csv",header = F,stringsAsFactors = F),use.names = F))
kw_analyte <- unique(unlist(read.csv("analyte.csv",header = F,stringsAsFactors = F),use.names = F))
