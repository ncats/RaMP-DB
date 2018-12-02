processData <- function(conpass=NULL,
                        dbname="ramp",username="root",
                        host = "localhost") {
  
  if(is.null(conpass)) {
    stop("Please define the password for the mysql connection")
  }
  
  # get all rows form analytehaspathway
  query <- "select * from analytehaspathway"
  con <- DBI::dbConnect(RMySQL::MySQL(), user = "root",
                        password = "password",
                        dbname = "ramp",
                        host = "localhost")
  
  allRampIds <- DBI::dbGetQuery(con,query)
  
  
  if(is.null(allRampIds)) {
    
    stop("Data doesn't exist")
    
  } else {
    
    
    # length(unique(allRampIds$pathwayRampId)) 
    # total unique pathwayRampIDs - 51,526
    
    # data frames for metabolites with pathawayRampID, Freq based  on Source(kegg, reactome, wiki, hmdb)
    
    allRampId_C <- allRampIds[grep("RAMP_C", allRampIds$rampId), ]
    unique_allRampId_C <- unique(allRampId_C[,c("rampId", "pathwayRampId")])
    unique_pathwayRampId_source <- unique(allRampId_C[,c("pathwayRampId", "pathwaySource")])
    
    # length(unique_pathwayRampId_source$pathwayRampId) - 36,039
    
    freq_unique_allRampId_C <- as.data.frame(table(unique_allRampId_C[,"pathwayRampId"]))
    
    # length of total RAMP_C pathwayRampIDs - 36,039
    
    names(freq_unique_allRampId_C)[1] = 'pathwayRampId'
    merge_Pathwayfreq_source <- merge(freq_unique_allRampId_C, unique_pathwayRampId_source, by="pathwayRampId")
    
    # subset data based on source -  kegg, reactome, wiki, hmdb
    
    kegg_metab <- subset(merge_Pathwayfreq_source, merge_Pathwayfreq_source$pathwaySource == "kegg")
    #length(kegg_metab$pathwayRampId) - 264
    reactome_metab <- subset(merge_Pathwayfreq_source, merge_Pathwayfreq_source$pathwaySource == "reactome")
    #length(reactome_metab$pathwayRampId) - 1866
    wiki_metab <- subset(merge_Pathwayfreq_source, merge_Pathwayfreq_source$pathwaySource == "wiki")
    #length(wiki_metab$pathwayRampId) - 213
    hmdb_metab <- subset(merge_Pathwayfreq_source, merge_Pathwayfreq_source$pathwaySource == "hmdb")
    #length(hmdb_metab$pathwayRampId) - 33,696
    
    
    
    # data frames for Genes with pathawayRampID, Freq based  on Source(kegg, reactome, wiki, hmdb)
    
    allRampId_G <- allRampIds[grep("RAMP_G", allRampIds$rampId), ]
    unique_allRampId_G <- unique(allRampId_G[,c("rampId", "pathwayRampId")])
    unique_pathwayG_source <- unique(allRampId_G[,c("pathwayRampId", "pathwaySource")])
    
    # length(unique(unique_pathwayG_source$pathwayRampId)) - 51,461
    
    freq_unique_allRampId_G <- as.data.frame(table(unique_allRampId_G[,"pathwayRampId"]))
    
    # length of total RAMP_G pathwayRampIDs - 51,461
    
    names(freq_unique_allRampId_G)[1] = 'pathwayRampId'
    merge_PathwayG_source <- merge(freq_unique_allRampId_G, unique_pathwayG_source, by="pathwayRampId")
    
    # subset data based on source -  kegg, reactome, wiki, hmdb
    
    kegg_gene <- subset(merge_PathwayG_source, merge_PathwayG_source$pathwaySource == "kegg")
    #length(kegg_gene$pathwayRampId) - 323
    reactome_gene <- subset(merge_PathwayG_source, merge_PathwayG_source$pathwaySource == "reactome")
    #length(reactome_gene$pathwayRampId) - 2155
    wiki_gene <- subset(merge_PathwayG_source, merge_PathwayG_source$pathwaySource == "wiki")
    #length(wiki_gene$pathwayRampId) - 408
    hmdb_gene <- subset(merge_PathwayG_source, merge_PathwayG_source$pathwaySource == "hmdb")
    #length(hmdb_gene$pathwayRampId) - 48575
    
    save(kegg_metab,
         reactome_metab,
         wiki_metab,
         hmdb_metab,
         kegg_gene,
         reactome_gene,
         wiki_gene,
         hmdb_gene, 
         file = "FT_data.Rdata" )
    #head(load(file="FT_data.Rdata"))
  }
}
