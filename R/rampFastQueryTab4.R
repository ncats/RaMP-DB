#' Return analytes that has catalyzation relation with given analyte
#' 
#' @param synonym string or data.frame contains given analyte
#' @param con a connection object returned from the function connectToRaMP()
#' @return a data.frame that contains search results
rampFastOneCata <- function(synonym,con) {
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
  }
  list_metabolite <- unique(list_metabolite)
  list_metabolite <- sapply(list_metabolite,shQuote)
  list_metabolite <- paste(list_metabolite,collapse = ",")
  print(list_metabolite)
  query1 <- paste0("select Synonym as analyte1,rampId,geneOrCompound as type1 from analytesynonym where Synonym in (",
                   list_metabolite,");")
  df1<- DBI::dbGetQuery(con,query1)
  print(df1$rampId)
  df_c <- df_g <- NULL
  mdf_c <- mdf_g <- NULL
  if(length(grep("RAMP_C",df1$rampId)!=0)){
    df_c <- df1[grep("RAMP_C",unique(df1$rampId)),]
    print("Get Compound ...")
    c_id <- unique(df_c$rampId)
    if(length(c_id) == 0){
      message("No searching result")
      return(NULL)
    }
    print(length(c_id))
    c_id <- sapply(c_id,shQuote)
    c_id <- paste(c_id,collapse = ",")
    query_c <- paste0("select rampCompoundId as rampId,rampGeneId as rampId2 from catalyzed where rampCompoundId in (",
                      c_id,");")
    print("Geting gene Id from Compound Id ...")
    df_c2 <- DBI::dbGetQuery(con,query_c)
    if(nrow(df_c2) == 0){
      message("No searching result")
      return(NULL)
    }
    print("Getting analyte from gene Id ...")
    analyte2_list <- unique(df_c2$rampId2)
    analyte2_list <- sapply(analyte2_list,shQuote)
    analyte2_list <- paste(analyte2_list,collapse = ",")
    query2 <- paste0("select Synonym as analyte2,rampId as rampId2,geneOrCompound as type2 from analyteSynonym 
                     where rampId in (",analyte2_list,");")
    df_c3 <- DBI::dbGetQuery(con,query2)
    query3 <- paste0("select sourceId,rampId as rampId2,IDtype from source where rampId in (",
                     analyte2_list,");")
    print("Get source ...")
    df_c4 <- DBI::dbGetQuery(con,query3)
    df_c4 <- unique(df_c4)
    
    print("merge 1")
    # mdf_c <- merge(df_c3,df_c4,all = T)
    mdf_c <- dplyr::left_join(df_c3,df_c4)
    print("merge 2")
    # mdf_c <- merge(mdf_c,df_c2,all.x = T)
    mdf_c <- dplyr::left_join(mdf_c,df_c2)
    print("merge 3")
    # mdf_c <- merge(mdf_c,df_c,all.x = T)
    mdf_c <- dplyr::left_join(mdf_c,df_c)
    mdf_c <- unique(mdf_c[,c(1,3,4,5,7,8)])
  }
  
  if(length(grep("RAMP_G",df1$rampId))!=0){
    print("Also find gene inside")
    df_g <- df1[grep("RAMP_G",unique(df1$rampId)),]
    print("Get gene ...")
    g_id <- df_g$rampId
    g_id <- sapply(unique(df_g),shQuote)
    g_id <- paste(g_id,collapse = ",")
    query_g <- paste0("select rampGeneId as rampId,rampCompoundId as rampId2 from catalyzed where rampGeneId in (",
                      g_id,");")
    df_g2 <- DBI::dbGetQuery(con,query_g)
    if(nrow(df_g2) == 0){
      message("No searching result for gene")
      return(NULL)
    }
    analyte2_list <- df_g2$rampId2
    analyte2_list <- sapply(analyte2_list,shQuote)
    analyte2_list <- paste(analyte2_list,collapse = ",")
    query2 <- paste0("select Synonym as analyte2,rampId as rampId2,geneOrCompound as type2 from analyteSynonym 
                     where rampId in (",analyte2_list,");")
    print(query2)
    df_g3 <-DBI::dbGetQuery(con,query2)
    if(nrow(df_g3))
    query3 <- paste0("select sourceId,rampId as rampId2,IDtype from source where rampId in (",
                     analyte2_list,");")
    df_g4 <- DBI::dbGetQuery(con,query3)
    
    # mdf_g <- merge(df_g3,df_g4,all = T)
    mdf_g <- dplyr::left_join(df_g3,df_g4)
    # mdf_g <- merge(mdf_g,df_g2,all.x = T)
    mdf_g <- dplyr::left_join(mdf_g,df_g2)
    # mdf_g <- merge(mdf_g,df_g,all.x = T)
    mdf_g <- dplyr::left_join(mdf_g,df_g)
    mdf_g <- unique(mdf_g[,c(1,3,4,5,7,8)])
  }
  mdf <- rbind(mdf_c,mdf_g)
  print("Done ...")
  print("timing ...")
  print(proc.time() - now)
  return(mdf)
}

#' Return analytes that has catalyzation relation with given list of analytes
#' 
#' @param synonym string or data.frame contains given analyte
#' @param con a connection object returned from the function connectToRaMP()
#' @return a data.frame that contains search results
rampFastMulCata <- function(synonym,con) {
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
  }
  list_metabolite <- unique(list_metabolite)
  list_metabolite <- sapply(list_metabolite,shQuote)
  list_metabolite <- paste(list_metabolite,collapse = ",")
  print(list_metabolite)
  query1 <- paste0("select Synonym as analyte1,rampId,geneOrCompound as type1 from analytesynonym where Synonym in (",
                   list_metabolite,");")
  df1<- DBI::dbGetQuery(con,query1)
  print(df1$rampId)
  df_c <- df_g <- NULL
  mdf_c <- mdf_g <- NULL
  if(length(grep("RAMP_C",df1$rampId)!=0)){
    df_c <- df1[grep("RAMP_C",unique(df1$rampId)),]
    print("Get Compound ...")
    c_id <- unique(df_c$rampId)
    if(length(c_id) == 0){
      message("No searching result")
      return(NULL)
    }
    print(length(c_id))
    c_id <- sapply(c_id,shQuote)
    c_id <- paste(c_id,collapse = ",")
    query_c <- paste0("select rampCompoundId as rampId,rampGeneId as rampId2 from catalyzed where rampCompoundId in (",
                      c_id,");")
    print("Geting gene Id from Compound Id ...")
    df_c2 <- DBI::dbGetQuery(con,query_c)
    if(nrow(df_c2) == 0){
      message("No searching result")
      df_c2 <- data.frame(rampId2 = "Not Found")
    }
    print("Getting analyte from gene Id ...")
    analyte2_list <- unique(df_c2$rampId2)
    analyte2_list <- sapply(analyte2_list,shQuote)
    analyte2_list <- paste(analyte2_list,collapse = ",")
    query2 <- paste0("select Synonym as analyte2,rampId as rampId2,geneOrCompound as type2 from analyteSynonym 
                     where rampId in (",analyte2_list,");")
    df_c3 <- DBI::dbGetQuery(con,query2)
    query3 <- paste0("select sourceId,rampId as rampId2,IDtype from source where rampId in (",
                     analyte2_list,");")
    print("Get source ...")
    df_c4 <- DBI::dbGetQuery(con,query3)
    df_c4 <- unique(df_c4)
    
    print("merge 1")
    # mdf_c <- merge(df_c3,df_c4,all = T)
    mdf_c <- dplyr::left_join(df_c3,df_c4)
    print("merge 2")
    # mdf_c <- merge(mdf_c,df_c2,all.x = T)
    mdf_c <- dplyr::left_join(mdf_c,df_c2)
    print("merge 3")
    # mdf_c <- merge(mdf_c,df_c,all.x = T)
    mdf_c <- dplyr::left_join(mdf_c,df_c)
    mdf_c <- unique(mdf_c[,c(1,3,4,5,7,8)])
  }
  
  if(length(grep("RAMP_G",df1$rampId))!=0){
    print("Also find gene inside")
    df_g <- df1[grep("RAMP_G",unique(df1$rampId)),]
    print("Get gene ...")
    g_id <- df_g$rampId
    g_id <- sapply(unique(df_g),shQuote)
    g_id <- paste(g_id,collapse = ",")
    query_g <- paste0("select rampGeneId as rampId,rampCompoundId as rampId2 from catalyzed where rampGeneId in (",
                      g_id,");")
    df_g2 <- DBI::dbGetQuery(con,query_g)
    if(nrow(df_g2) == 0){
      message("No searching result for gene")
      
    }
    analyte2_list <- df_g2$rampId2
    analyte2_list <- sapply(analyte2_list,shQuote)
    analyte2_list <- paste(analyte2_list,collapse = ",")
    query2 <- paste0("select Synonym as analyte2,rampId as rampId2,geneOrCompound as type2 from analyteSynonym 
                     where rampId in (",analyte2_list,");")
    print(query2)
    df_g3 <-DBI::dbGetQuery(con,query2)
    query3 <- paste0("select sourceId,rampId as rampId2,IDtype from source where rampId in (",
                     analyte2_list,");")
    df_g4 <- DBI::dbGetQuery(con,query3)
    
    # mdf_g <- merge(df_g3,df_g4,all = T)
    mdf_g <- dplyr::left_join(df_g3,df_g4)
    # mdf_g <- merge(mdf_g,df_g2,all.x = T)
    mdf_g <- dplyr::left_join(mdf_g,df_g2)
    # mdf_g <- merge(mdf_g,df_g,all.x = T)
    mdf_g <- dplyr::left_join(mdf_g,df_g)
    mdf_g <- unique(mdf_g[,c(1,3,4,5,7,8)])
  }
  mdf <- rbind(mdf_c,mdf_g)
  print("Done ...")
  print("timing ...")
  print(proc.time() - now)
  return(mdf)
}


#' Generate data.frame from given files
#' 
#' identifing the file type, then it returns table output to 
#' shiny renderTable function as preview of searching data
#' 
#' @param infile a file object given from files 
#' 
#' @return a data.frame either from multiple csv file
#' or search through by a txt file.
rampFileOfAnalytes_tab4 <- function(infile){
  name <- infile[[1,'name']]
  summary <- data.frame(pathway  = character(0),id = character(0),
                        source = character(0),metabolite = character(0))
  rampOut <- list()
  for (i in 1:length(infile[,1])){
    rampOut <- readLines(infile[[i,'datapath']])
    summary <- rampFastOneCata(rampOut,T)
  }
  return(summary)
}
