library(RMySQL)
library(VennDiagram)
now <- proc.time()
con <- dbConnect(MySQL(), user = "root", password = "Ehe131224", dbname = "mathelabramp")

query <- "select pathwayName from pathway;"

pathway_list <- dbGetQuery(con,query)

pathway_list <- unlist(pathway_list)
dbGetQuery(con,"select * from analyte limit 5;")
tot_metabolites <- dbGetQuery(con,"select count(*) from analytesynonym;")
tb <- data.frame(inpathway=NULL,outpathway=NULL)

for (pathway in pathway_list){
  tot_in_pathway <- dbGetQuery(con,paste0("select count(*) from analyte 
                                            where rampId in (select rampId from 
                                          analytehaspathway where 
                                          pathwayRampId in (select pathwayRampId 
                                          from pathway where pathwayName = \"",
                                          pathway,"\"));"))
  tot_out_pathways <- tot_metabolites - tot_in_pathway
  item <- data.frame(inpathway = tot_in_pathway,outpathway = tot_out_pathways)
  row.names(item) <- pathway
  tb <- rbind(tb,item)
  print(tot_in_pathway)
}
print("timing ...")
print(proc.time() - now)
colnames(tb) <- c("inPathway","outPathway")
write.csv(tb,"dataFisherTest.csv")

df1 <- read.table("dataFisherTest.csv")
df2 <- read.csv("dataFisherTest.csv")

df2[grepl("Metabolism",df2),]

# generate Keyword search
query1 <- "select synonym from analytesynonym;"
query2 <- "select pathwayName from pathway;"
query3 <- "select commonName from ontology;"
query4 <- "select distinct * from analyte;"
analyte <- dbGetQuery(con,query1)
pathway <- dbGetQuery(con,query2)
ontology <- dbGetQuery(con,query3)
nrow(unique(analyte))
length(unique(analyte$rampId))

analyte2 <- dbGetQuery(con,query4)
id2 <- analyte2$rampId
length(id2[grepl("RAMP_C",id2)])
write.csv(analyte,"analyte.csv",row.names = F)
write.csv(pathway,"pathway.csv",row.names = F)
write.csv(ontology,"biofluid.csv",row.names = F)

# try Keyword search.
kw_analyte[grepl("VitaminE",kw_analyte)]
kw <- kw_analyte[grep("VitaminA",kw_analyte)]

kw
kw <- kw[order(nchar(kw),kw)]
kw
kw2 <- kw_pathway[grepl("a",kw_pathway)]
kw2
