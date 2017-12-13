# The global variables that are important to keywords search and fisher test.
FisherPathwayTable <- read.csv("dataFisherTest.csv")
FisherPathwayTable_m <- read.csv("FisherTestDataMetabolites.csv")
FisherPathwayTable_g <- read.csv("FisherTestDataGenes.csv")
# Unlist automatically select first column, the first column should
# be the data containing keywords
kw_biofluid <- unique(unlist(read.csv("biofluid.csv",header = F,stringsAsFactors = F),use.names = F))
kw_pathway <- unique(unlist(read.csv("pathway.csv",header = F,stringsAsFactors = F),use.names = F))
kw_analyte <- unique(unlist(read.csv("analyte.csv",header = F,stringsAsFactors = F),use.names = F))
kw_source <- unique(unlist(read.csv("source.csv",header = F,stringsAsFactors = F),use.names = F))