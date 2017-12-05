#' Function that sends query with tryCatch
#' Return the customized error messages
#' @param query
#' @export
rampReadQuery <- function(query){
  out <- tryCatch({
    DBI::dbGetQuery(con,query)
  },
  error = function(cond){
    message("No searching result in one of the query ...")
    message(cond)
    return(NULL)
  },
  finally = {
    message("Done for searching ...")
    stop("No Searching Result")
  }
  )
  return(out)
}

#' Use fast algorithm to search through database to find analytes
#' 
#' @param synonym a string,or a data.frame that has all analytes' name user
#' want to search for
#' @param options a logic value used to define output (deprecate)
#' @return a data.frame that contains all search results
#' @export
rampFastMetaInPath <- function(synonym,options = TRUE){
  now <- proc.time()
  # con <- dbConnect(MySQL(), user = "root", password = "Ramp340!", dbname = "mathelabramp")
  # on.exit(dbDisconnect(con))
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
  } else {
    return("Wrong Data Format")
  }

    list_metabolite <- sapply(list_metabolite,shQuote)
    list_metabolite <- paste(list_metabolite,collapse = ",")
    
    query1 <- paste0("select distinct Synonym,rampId from analytesynonym where Synonym in (",
                     list_metabolite,");")
    df1<- DBI::dbGetQuery(con,query1)
    colnames(df1) <- c("Synonym1","rampId1")
    query2 <- paste0("select pathwayRampId,rampId from analytehaspathway where 
                     rampId in (select rampId from analytesynonym where Synonym in (",
                     list_metabolite,"));")
    df2 <- DBI::dbGetQuery(con,query2)
    df2<- unique(df2)
    colnames(df2) <- c("pathwayRampId1","rampId1")
    pid_list <- df2[,1]
    pid_list <- sapply(pid_list,shQuote)
    pid_list <- paste(pid_list,collapse = ",")
    query3 <- paste0("select rampId,pathwayRampId from analytehaspathway where pathwayRampId in (",
                     pid_list,");")
    df3 <- DBI::dbGetQuery(con,query1)
    df3 <- unique(df3)
    colnames(df3) <- c("rampId2","pathwayRampId1")
    cid_list <- df3[,1]
    cid_list <- sapply(cid_list,shQuote)
    cid_list <- paste(cid_list,collapse = ",")
    query4 <- paste0("SELECT Synonym,geneOrCompound,rampId FROM analytesynonym 
                     WHERE rampID IN (",cid_list,");")
    df4 <- dbGetQuery::dbGetQuery(con,query4)
    df4<- unique(df4)
    colnames(df4) <- c("Synonym2","geneOrCompound","rampId2")
    query5 <- paste0("select rampId,sourceId,IDtype from source where rampId in(",
                     cid_list,");")
    df5 <- DBI::dbGetQuery(con,query5)
    if(is.null(df5)){
      return("No searching result ...")
    }
    df5 <- unique(df5)
    colnames(df5) <- c("rampId2","sourceId","IDtype")
    mdf <- merge(df4,df5,all = T)
    mdf <- mdf[!is.na(mdf[,'Synonym2']),]
    mdf2 <- merge(mdf,df3,all.x = T)
    mdf3 <- merge(mdf2,df2,all.x = T)
    mdf4 <- merge(mdf3,df1,all.x=T)
    mdf4 <- unique(mdf4)
    print("Timing ..")
    print(proc.time() - now)
    mdf4 <- mdf4[,4:8]
    colnames(mdf4) <- c("metabolite from same pathway","type","sourceId","source",
                        "given metabolites")
    return(mdf4[,1:4])
}

#' Query: given a synonym, returns other genes or metabolites from the same pathway
#' @param synonym synonym (character string)
rampFastMetaInPath2 <- function(synonym){
  now <- proc.time()
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
  } else {
    return("Wrong Data Format")
  }
  
  list_metabolite <- sapply(list_metabolite,shQuote)
  list_metabolite <- paste(list_metabolite,collapse = ",")
  
  query1 <- paste0("select distinct Synonym,rampId from analytesynonym where Synonym in (",
                   list_metabolite,");")
  df1<- DBI::dbGetQuery(con,query1)
  colnames(df1) <- c("Synonym1","rampId1")
  query2 <- paste0("select pathwayRampId,rampId from analytehaspathway where 
                   rampId in (select rampId from analytesynonym where Synonym in (",
                   list_metabolite,"));")
  df2 <- DBI::dbGetQuery(con,query2)
  df2<- unique(df2)
  colnames(df2) <- c("pathwayRampId1","rampId1")
  pid_list <- df2[,1]
  pid_list <- sapply(pid_list,shQuote)
  pid_list <- paste(pid_list,collapse = ",")
  query3 <- paste0("select rampId,pathwayRampId from analytehaspathway where pathwayRampId in (",
                   pid_list,");")
  
  df3 <- DBI::dbGetQuery(con,query3)
  df3 <- unique(df3)
  colnames(df3) <- c("rampId2","pathwayRampId1")
  cid_list <- df3[,1]
  cid_list <- sapply(cid_list,shQuote)
  cid_list <- paste(cid_list,collapse = ",")
  query4 <- paste0("SELECT Synonym,geneOrCompound,rampId FROM analytesynonym 
                   WHERE rampID IN (",cid_list,");")
  df4 <- DBI::dbGetQuery(con,query4)
  df4<- unique(df4)
  colnames(df4) <- c("Synonym2","geneOrCompound","rampId2")
  query5 <- paste0("select rampId,sourceId,IDtype from source where rampId in(",
                   cid_list,");")
  df5 <- DBI::dbGetQuery(con,query5)
  if(is.null(df5)){
    return("No searching result ...")
  }
  df5 <- unique(df5)
  colnames(df5) <- c("rampId2","sourceId","IDtype")
  mdf <- merge(df4,df5,all = T)
  mdf <- mdf[!is.na(mdf[,'Synonym2']),]
  mdf2 <- merge(mdf,df3,all.x = T)
  mdf3 <- merge(mdf2,df2,all.x = T)
  mdf4 <- merge(mdf3,df1,all.x=T)
  mdf4 <- unique(mdf4)
  print("Timing ..")
  print(proc.time() - now)
  mdf4 <- mdf4[,4:8]
  colnames(mdf4) <- c("metabolite from same pathway","type","sourceId","source",
                      "given metabolites")
  
  mdf5 <- mdf4[,c(1,4)]
}
