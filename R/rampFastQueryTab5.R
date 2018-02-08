

rampFastOnto <- function(analytes,conpass = NULL,
                         dbname = 'mathelabramp',
                         host = 'localhost',
                         username = 'root',
                         analytesOrOntology = 'analytes'){
  if(!(analytesOrOntology %in% c('analytes','biofluids')))
  if(is.null(conpass)) {
    stop("Please define the password for the mysql connection")
  }
  now <- proc.time()
  if(is.character(analytes)){
    if(grepl("\n",analytes)[1]){
      list_metabolite <- strsplit(analytes,"\n")
      list_metabolite <- unlist(list_metabolite)
    } else if(grepl(",",analytes)[1]){
      list_metabolite <- strsplit(analytes,"\n")
      list_metabolite <- unlist(list_metabolite)
    } else {
      list_metabolite <- analytes
    }
  } else if(is.data.frame(analytes)){
    list_metabolite <- unlist(analytes)
  }
  list_metabolite <- unique(list_metabolite)
  list_metabolite <- sapply(list_metabolite,shQuote)
  list_metabolite <- paste(list_metabolite,collapse = ",")
  
  con <- DBI::dbConnect(RMySQL::MySQL(),
                        dbname = dbname,
                        username = username,
                        host = host,
                        password = conpass)
  sql <- paste0('select * from source where sourceId in (',list_metabolite,');')
  df <- DBI::dbGetQuery(con,sql)
  print(colnames(df))
  DBI::dbDisconnect(con)
  if(nrow(df) == 0) {
    message('No searching result since this source id 
            is not existed in source table')
    return(NULL)
  }
  rampid <- unique(df$rampId)
  rampid <- sapply(rampid,shQuote)
  rampid <- paste(rampid,collapse = ',')
  con <- DBI::dbConnect(RMySQL::MySQL(),
                        dbname = dbname,
                        username = username,
                        host = host,
                        password = conpass)
  
  sql <- paste0('select * from analytehasontology where rampCompoundId in (',
                rampid,');')
  df2 <- DBI::dbGetQuery(con,sql)
  if(nrow(df2) == 0) {
    message('No searching result because these metabolites are not linked to ontology')
    return(NULL)
  }
  DBI::dbDisconnect(con)
  rampontoid <- unique(df2$rampOntologyIdLocation)
  rampontoid <- sapply(rampontoid,shQuote)
  rampontoid <- paste(rampontoid,collapse = ',')
  sql <- paste0('select * from ontology where rampOntologyIdLocation in (',
                rampontoid,');')
  con <- DBI::dbConnect(RMySQL::MySQL(),
                        dbname = dbname,
                        username = username,
                        host = host,
                        password = conpass)
  df3 <- DBI::dbGetQuery(con,sql)
  print(colnames(df3))
  DBI::dbDisconnect(con)
  mdf <- unique(merge(df3,df2,all.x=T))
  mdf <- unique(merge(mdf,df,all.x = T,by.x = 'rampCompoundId',
                      by.y= 'rampId'))
  mdf <- mdf[c('sourceId','commonName.x','biofluidORcellular')]
  colnames(mdf) <- c('Metabolites','Ontology','Ontology_Type')
  return(mdf)
}

rampFastMetaFromOnto <- function(ontology,conpass = NULL,
                                 dbname = 'mathelabramp',
                                 host = 'localhost',
                                 username = 'root',){
  
}