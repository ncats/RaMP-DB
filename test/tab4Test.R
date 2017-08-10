library(RMySQL)

con <- dbConnect(MySQL(), user = "root", password = "Ramp340!", dbname = "mathelabramp")
dbDisconnect(con)


string <- c("VitaminA","VitaminB1","VitaminC","VitaminD","VitaminE","GTG",
            "Glucose","GST3","GABARAPL2","ATG10","GABARAP")
string <- "triglyceride"
string <- c("GABARAPL2","ATG10","GABARAP")

mdf <- rampFastCata(string,T)
mdf[,1:4]
df3 <- rampFastCata(string)
nrow(unique(df3))
nrow(unique(mdf))

# manually 
string <- sapply(string,shQuote)
string
string <- paste0(string,collapse = ",")
now <- proc.time()
query <- paste0("select Synonym,rampId from analyteSynonym where Synonym in (",string,");")
query


df1 <- dbGetQuery(con,query)

c_id <- unique(df1$rampId)
g_id <- c_id[grepl("RAMP_G",c_id)]
c_id <- c_id[grepl("RAMP_C",c_id)]
c_id <- sapply(c_id,shQuote)
c_id <- paste(c_id,collapse = ",")


query2 <- paste0("select rampCompoundId as rampId,rampGeneId as rampId2 from catalyzed where 
                 rampCompoundId in (",c_id,");")
query2
df2 <- dbGetQuery(con,query2)

id2 <-  unique(df2$rampId2)
id2 <- sapply(id2,shQuote)
id2 <- paste0(id2,collapse = ",")
query_c <- paste0("select Synonym,rampId from analyteSynonym where rampId in (",id2,");")
df3 <- dbGetQuery(con,query_c)
