#' Find all synonym from a given metabolite's name
#' @param full bool if return whole data.frame
#' @param find_synonym bool if find all synonyms or just return same synonym
#' as input (there are some common synonyms that will mess up whole searching)
#' @export
rampFindSynonymFromSynonym <- function(synonym,full = F,find_synonym = F){
  if(is.character(synonym)){
    if(grepl("\n",synonym)[1]){
      list_metabolite <- strsplit(synonym,"\n")
      list_metabolite <- unlist(list_metabolite)
    } else if(grepl(",",synonym)[1]){
      list_metabolite <- strsplit(synonym,",")
      list_metabolite <- unlist(list_metabolite)
    } else {
      list_metabolite <- synonym
    }
  } else if(is.data.frame(synonym)){
    list_metabolite <- unlist(synonym)
  } else{
    message("Wrong Format of argument")
    return(NULL)
  }
  if(!find_synonym){
    message("Dont Find synonym due to common synonyms (Triglyceride?).")
    list_metabolite <- unique(list_metabolite)
    list_metabolite <- sapply(list_metabolite,shQuote)
    list_metabolite <- paste(list_metabolite,collapse = ",")
    query <- paste0("select synonym as origins,rampId from analyteSynonym where Synonym in(",
                    list_metabolite,
                    ");")
    df1 <- DBI::dbGetQuery(con,query)
    return(df1)
  }
  list_metabolite <- unique(list_metabolite)
  list_metabolite <- sapply(list_metabolite,shQuote)
  list_metabolite <- paste(list_metabolite,collapse = ",")
  query <- paste0("select synonym as origins,rampId from analyteSynonym where Synonym in(",
                  list_metabolite,
                  ");")
  df1 <- DBI::dbGetQuery(con,query)
  
  rampid <- df1$rampId
  rampid <- sapply(rampid,shQuote)
  rampid <- paste(rampid,collapse = ",")
  query <- paste0("select * from analyteSynonym where rampId in(",rampid,");")
  df2 <- DBI::dbGetQuery(con,query)
  df2 <- merge(df1,df2)
  if(full){
    return(df2)
  }
  synonym <- df2$Synonym
  print('Hello world')
  synonym
}
#' Find all source from given list of RaMP Ids
#' @param rampId could be a data frame return by rampFindSynonymFromSynonym
#'  containing all information related to synonym. Or can be a list of 
#'  rampId
#'  @param full return whole searching result or not.
#' @export
rampFindSourceFromId <- function(rampId=NULL,full = T){
  if(is.data.frame(rampId)){
    list_id <- rampId$rampId
  } else if(is.character(rampId)){
    if(grepl("\n",rampId)[1]){
      list_id <- strsplit(rampId,"\n")
      list_id <- unlist(list_id)
    } else if(grepl(",",rampId)[1]){
      list_id <- strsplit(rampId,",")
      list_id <- unlist(list_id)
    } else {
      list_id <- rampId
    }
  } else{
    message("Wrong format of input")
    return(NULL)
  }
  list_id <- unique(list_id)
  list_id <- sapply(list_id,shQuote)
  list_id <- paste(list_id,collapse = ",")
  query <- paste0("select * from source where rampId in (",list_id,");")
  df <- DBI::dbGetQuery(con,query)
  if(full){
    return(df)
  } else{
    return(df[,1])
  }
}


#' Fast search given a list of metabolites
#' @param sourceid a vector of synonym that need to be searched
#' @return a list contains all metabolits as name and pathway inside.
#' 
#' Apply famil function...
#' 
#' @export
rampFastPathFromSource<- function(sourceid,find_synonym = FALSE){
  # progress<- shiny::Progress$new()
  # progress$set(message = "Querying databases ...",value = 0)
  now <- proc.time()
  # con <- dbConnect(MySQL(), user = "root", password = "Ramp340!", dbname = "mathelabramp")
  # on.exit(dbDisconnect(con))
  # find synonym
  
  #synonym <- rampFindSynonymFromSynonym(synonym,find_synonym=find_synonym)
  
  list_metabolite <- unique(sourceid)
  list_metabolite <- sapply(list_metabolite,shQuote)
  list_metabolite <- paste(list_metabolite,collapse = ",")
  query1 <- paste0("select * from source where sourceid in (",
                   list_metabolite,");")
  df1<- DBI::dbGetQuery(con,query1)
  colnames(df1)[1] <-"sourceId2"
  #return(df1)
  rampid <- df1$rampId
  rampid <- sapply(rampid,shQuote)
  rampid <- paste(rampid,collapse = ",")
  query2 <- paste0("select * from analytehaspathway where 
                   rampId in (",rampid,");")
  df2 <- DBI::dbGetQuery(con,query2)
  #return(df2)
  id_list <- unique(df2$pathwayRampId)
  id_list <- sapply(id_list,shQuote)
  id_list <- paste(id_list,collapse = ",")
  print(id_list)
  query3 <- paste0("select * from pathway where pathwayRampId in (",
                   id_list,");")
  df3 <- DBI::dbGetQuery(con,query3)
  #return(df3)
  mdf <- merge(df3,df2,all.x=T)
  mdf <- merge(mdf,df1,all.x = T)
  mdf <- unique(mdf)
  print("timing ...")
  print(proc.time()- now)
  return(unique(mdf[,c(6,4,3,7)]))
}
