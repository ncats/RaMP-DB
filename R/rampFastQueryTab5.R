#' Search by fast algorithm to find analyte or biofluid location
#' 
#' @param string a string or data.frame that contains all given analyte or
#' biofluid location
#' @param analyteOrBiofluid a string or NULL value defined output format
#' @return if analyteOrBiofluid == NULL, it returns a list, otherwise return
#' a data.frame that has analyte or biofluid
#' @export 
rampFastBiofluid <- function(string,analyteOrBiofluid = NULL){
  biof <- vector()
  analyte <- vector()
  # check and distribute individual item to either biofluid or analyte
  if(is.character(string)){
    if(grepl("\n",string)[1]){
      list_metabolite <- strsplit(string,"\n")
      list_metabolite <- unlist(list_metabolite)
    } else if(grepl(",",string)[1]){
      list_metabolite <- strsplit(string,",")
      list_metabolite <- unlist(list_metabolite)
    } else {
      list_metabolite <- string
    }
  } else if(is.data.frame(string)){
    list_metabolite <- unlist(string)
  }
  for(i in 1:length(list_metabolite)){
    if(is.element(list_metabolite[i],kw_analyte)){
      analyte <- c(analyte,list_metabolite[i])
    } else if (is.element(list_metabolite[i],kw_biofluid)) {
      biof <- c(biof,list_metabolite[i])
    }
  }
  mdf <- list()
  if(length(analyte) > 0){
    analyte <- sapply(analyte,shQuote)
    analyte <- paste(analyte,collapse = ",")
    query_a1 <- paste("select synonym,rampId from analytesynonym where synonym in (",analyte,");")
    dfa1 <- DBI::dbGetQuery(con,query_a1)
    c_id <- dfa1$rampId
    c_id <- unique(c_id)
    c_id <- sapply(c_id,shQuote)
    c_id <- paste(c_id,collapse = ",")
    query_a2 <- paste("select rampCompoundId as rampId,rampOntologyIdLocation as rampId2 from analytehasontology where rampCompoundId in (",
                      c_id,");")
    query_a2
    dfa2 <- DBI::dbGetQuery(con,query_a2)
    oid <- unique(dfa2$rampId2)
    oid <- sapply(oid,shQuote)
    oid <- paste(oid,collapse = ",")
    query_a3 <- paste0("select commonName as biofluidLocation,biofluidORcellular 
                       as type,rampOntologyIdLocation as rampId2 from ontology where rampOntologyIdLocation in (",
                       oid,");")
    dfa3 <- DBI::dbGetQuery(con,query_a3)
    mdfa <- dplyr::left_join(dfa3,dfa2)
    mdfa <- dplyr::left_join(mdfa,dfa1)
    mdfa <- unique(mdfa[,c(1,2,5)])
    mdf[['analyte']] <- mdfa
  }
  
  if(length(biof) > 0){
    biof <- sapply(biof,shQuote)
    biof <- paste(biof,collapse = ",")
    query_b1 <- paste0("select commonName as biofluidLocation,rampOntologyIdLocation as rampId from ontology 
                       where commonName in (",biof,");")
    dfb1 <- DBI::dbGetQuery(con,query_b1)
    oid <- unique(dfb1$rampId)
    oid <- sapply(oid,shQuote)
    oid <- paste(oid,collapse = ",")
    query_b2 <- paste0("select rampOntologyIdLocation as rampId,rampCompoundId as 
                       rampId2 from analytehasontology where rampOntologyIdLocation 
                       in (",oid,");")
    dfb2 <- DBI::dbGetQuery(con,query_b2)
    cid <- unique(dfb2$rampId2)
    cid <- sapply(cid,shQuote)
    cid <- paste(cid,collapse = ",")
    query_b3 <- paste0("select synonym,rampId as rampId2 from analytesynonym
                       where rampId in (",cid,");")
    query_b4 <- paste0("select rampId as rampId2,sourceId,IDtype from source where
                       rampId in (",cid,");")
    dfb3 <- DBI::dbGetQuery(con,query_b3)
    dfb4 <- DBI::dbGetQuery(con,query_b4)
    mdfb <- dplyr::left_join(dfb3,dfb4)
    mdfb <- dplyr::left_join(mdfb,dfb2)
    mdfb <- dplyr::left_join(mdfb,dfb1)
    mdfb <- mdfb[,c(1,3,4,6)] 
    mdfb <- unique(mdfb)
    mdf[['biofluid']] <- mdfb
  }
  print(" done for merging ... ")
  if(!is.null(analyteOrBiofluid) && analyteOrBiofluid =="analyteSynonym"){
    message("Ignore provided biofluid location")
    return(mdf[['analyte']])  
  } else if ( !is.null(analyteOrBiofluid) &&analyteOrBiofluid == "ontology"){
    message("Ignore provided metabolite or gene synonyms")
    return(mdf[['biofluid']])
  } else if (is.null(analyteOrBiofluid)){
    return(mdf)
  }
  
}
#' Generate search results based on given file
#' 
#' @param infile a file object given by shiny
#' @return search results of biofluid or analytes
#' @export
rampFileOfBiofluid <- function(infile) {
  name <- infile[[1,'name']]
  
  for (i in 1:length(infile[,1])){
    rampOut <- readLines(infile[[i,'datapath']])
    summary <- rampFastBiofluid(rampOut)
  }
  return(summary)
}
