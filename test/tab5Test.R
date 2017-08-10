library(RMySQL)

con <- dbConnect(MySQL(),
                 dbname = "mathelabramp",
                 user = "root",
                 password = "Ramp340!")
string <- c("VitaminA","Blood")
biof <- vector()
analyte <- vector()

# function run

df <- rampFastBiofluid(string)




# manually run
for(i in 1:length(string)){
  if(is.element(string[i],kw_analyte)){
    analyte <- c(analyte,string[i])
  } else if (is.element(string[i],kw_biofluid)) {
    biof <- c(biof,string[i])
  }
}
if(length(analyte) > 0){
  analyte <- sapply(analyte,shQuote)
  analyte <- paste(analyte,collapse = ",")
  query_a1 <- paste("select synonym,rampId from analytesynonym where synonym in (",analyte,");")
  dfa1 <- dbGetQuery(con,query_a1)
  c_id <- dfa1$rampId
  c_id <- unique(c_id)
  c_id <- sapply(c_id,shQuote)
  c_id <- paste(c_id,collapse = ",")
  query_a2 <- paste("select rampCompoundId as rampId,rampOntologyIdLocation as rampId2 from analytehasontology where rampCompoundId in (",
                    c_id,");")
  query_a2
  dfa2 <- dbGetQuery(con,query_a2)
  oid <- unique(dfa2$rampId2)
  oid <- sapply(oid,shQuote)
  oid <- paste(oid,collapse = ",")
  query_a3 <- paste0("select commonName as biofluidLocation,biofluidORcellular 
                     as type,rampOntologyIdLocation as rampId2 from ontology where rampOntologyIdLocation in (",
                     oid,");")
  dfa3 <- dbGetQuery(con,query_a3)
  mdfa <- dplyr::left_join(dfa3,dfa2)
  mdfa <- dplyr::left_join(mdf,dfa1)
  mdfa <- unique(mdf[,c(1,2,5)])
}

if(length(biof) > 0){
  biof <- sapply(biof,shQuote)
  biof <- paste(biof,collapse = ",")
  query_b1 <- paste0("select commonName as biofluidLocation,rampOntologyIdLocation as rampId from ontology 
                     where commonName in (",biof,");")
  dfb1 <- dbGetQuery(con,query_b1)
  oid <- unique(dfb1$rampId)
  oid <- sapply(oid,shQuote)
  oid <- paste(oid,collapse = ",")
  query_b2 <- paste0("select rampOntologyIdLocation as rampId,rampCompoundId as 
                     rampId2 from analytehasontology where rampOntologyIdLocation 
                     in (",oid,");")
  dfb2 <- dbGetQuery(con,query_b2)
  cid <- unique(dfb2$rampId2)
  cid <- sapply(cid,shQuote)
  cid <- paste(cid,collapse = ",")
  query_b3 <- paste0("select synonym,rampId as rampId2 from analytesynonym
                     where rampId in (",cid,");")
  query_b4 <- paste0("select rampId as rampId2,sourceId,IDtype from source where
                     rampId in (",cid,");")
  dfb3 <- dbGetQuery(con,query_b3)
  dfb4 <- dbGetQuery(con,query_b4)
  mdfb <- dplyr::left_join(dfb3,dfb4)
  mdfb <- dplyr::left_join(mdfb,dfb2)
  mdfb <- dplyr::left_join(mdfb,dfb1)
  mdfb <- mdfb[,c(1,3,4,6)] 
  mdfb <- unique(mdfb)
}
source <- dbGetQuery(con,"select * from source;")
anahasonto <- dbGetQuery(con,"select * from analytehasontology;")

nrow(source[grepl("hmdb",source$IDtype),])
query <- "select"
compound <- mdfb$synonym
length(unique(compound))
