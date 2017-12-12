#' Use fast search algorithm to find all pathways from given analytes
#' 
#' @param pathway a string or a data.fram that contains all pathways
#' @return a data.frame that contains all search results
rampFastMetaFromPath <- function(pathway){
  now <- proc.time()
  print("fired")
  # con <- dbConnect(MySQL(), 
  #                  user = "root", 
  #                  password = "Ramp340!", 
  #                  dbname = "mathelabramp")
  # on.exit(dbDisconnect(con))
  if(is.character(pathway)){
    if(grepl("\n",pathway)[1]){
      list_pathway <- strsplit(pathway,"\n")
      list_pathway <- unlist(list_pathway)
    } else if(grepl(",",pathway)[1]){
      list_pathway <- strsplit(pathway,",")
      list_pathway <- unlist(list_pathway)
    } else {
      list_pathway <- pathway
    }
  } else if(is.data.frame(pathway)){
    list_pathway <- unlist(pathway)
  } else {
    return("Wrong Data Format")
  }
  list_pathway <- sapply(list_pathway,shQuote)
  list_pathway <- paste(list_pathway,collapse = ",")
  query1 <- paste0("select pathwayName,pathwayRampId from pathway where pathwayName
                   in (",list_pathway,");")
  df1 <- DBI::dbGetQuery(con,query1)
  query2 <- paste0("select pathwayRampId,rampId from analytehaspathway where 
                   pathwayRampId in (select pathwayRampId from pathway where 
                   pathwayName in (",list_pathway,"));")
 
  df2 <- DBI::dbGetQuery(con,query2)
  cid_list <- unlist(df2[,2])
  cid_list <- sapply(cid_list,shQuote)
  cid_list <- paste(cid_list,collapse = ",")
  query3 <- paste0("select synonym,geneOrCompound,rampId from analytesynonym where rampId in (",
                   cid_list,");")
  df3 <- DBI::dbGetQuery(con,query3)
  query4 <- paste0("select rampId,sourceId,IDtype from source where rampId in (",
                   cid_list,");")
  df4 <- DBI::dbGetQuery(con,query4)
  mdf1 <- merge(df3,df4,all.x = T)
  mdf1 <- merge(mdf1,df2,all.x = T)
  mdf1 <- merge(mdf1,df1,all.x = T)
  mdf1 <- unique(mdf1)
  print("Timing ..")
  print(proc.time() - now)
  mdf1[,3:7]
}


#' Generate data.frame that has all search results from given files
#' 
#' identifing the file type, then it returns table output to 
#' shiny renderTable function as preview of searching data
#' 
#' @param infile a file object given from files 
#' 
#' @return a data.frame either from multiple csv file
#' or search through by a txt file.
rampFileOfPathways_tab2 <- function(infile){
  name <- infile[[1,'name']]
  summary <- data.frame(pathway  = character(0),id = character(0),
                        source = character(0),metabolite = character(0))
  rampOut <- list()
  for (i in 1:length(infile[,1])){
      rampOut <- readLines(infile[[i,'datapath']])
      summary <- rampFastMetaFromPath(rampOut)
  }
  return(summary)
}
