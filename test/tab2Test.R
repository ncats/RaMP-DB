library(RMySQL)
con <- dbConnect(MySQL(), user = "root", password = "Ehe131224", dbname = "mathelabramp")
dbDisconnect(con)
query <- "select pathwayName from pathway where rand()<0.01;"

random_list <- dbGetQuery(con,query)
write.table(random_list,"pathwayslist2.txt",quote = F,row.names = F)
string <- "Glucose-Alanine Cycle"
string <- random_list
df <- rampFastMetaFromPath(string)
