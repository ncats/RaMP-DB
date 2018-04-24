#' Find table of analyte has pathway from given pathway IDs
#' Aggregate ramp Id to ramp pathway Id
#' GC is C or G
#' @param pathwayRampId a vector of ramp Pathway ID
#' @param GC the analytes type that is either "C" for compound or "G" for gene
#' @param n minimum analytes of which pathway to considered computing overlap
#' @param username a string that specifies name of MySQL database
#' @param dbname a string that specifies database name of MySQL database
#' @param conpass a string that specifies password for database connection
#' @param host a string that specifes host for database connection
#' @return A list with pathway rampID as name, a vector of analytes from this pathway as content.
findAnalyteHasPathway <- function(pathwayRampId,GC = "C",n = 10,
                                  username ='root',
                                  dbname = 'ramp',
                                  conpass,
                                  host = 'localhost'){
  con <- DBI::dbConnect(RMySQL::MySQL(),
                   user = username,
                   dbname=dbname,
                   password = conpass,
                   host = host)
  on.exit(DBI::dbDisconnect(con))
  p_id <- unique(pathwayRampId)
  p_id <- sapply(p_id,shQuote)
  p_id <- paste(p_id,collapse = ",")
  query <-paste0("select * from analytehaspathway where pathwayRampId in (",
                 p_id,
                 ");")
  df <- DBI::dbGetQuery(con,
                   query)
  if(GC == 'both'){
    df2 <- aggregate(df$rampId,list(df$pathwayRampId),FUN = function(x){
      if(length(x) >= n){
        paste(x,collapse = ',')
      } else{
        x <- 0
      }
    })
  }
  else if (GC %in% c('G','C')){
    df2 <- aggregate(df$rampId,list(df$pathwayRampId),FUN = function(x){
      x <- x[grepl(paste0("RAMP_",GC,"_"),x)]
      if(length(x) >= n ){
        paste(x,collapse = ",")
      } else {
        x <- 0
      }
    })
  }
  fdf <- df2[df2$x!=0,]
  fdf2 <- data.frame(fdf[,-1],row.names = fdf[,1],stringsAsFactors = F)
  df.list <- setNames(split(fdf2, seq(nrow(fdf2))), rownames(fdf2))
  df.list <- lapply(df.list,FUN = function(x){
    text <- x[[1]]
    text <- strsplit(text,split = ",")
  })
  df.list <- lapply(df.list,unlist)
}
#'Compute overlaping matrix based on given list return by findAnalyteHasPathway()
#'
#'@param pathwayid a vector that has all ramp pathway id in the pathwaysWithAnalytes
#'@param pathwaysWithAnalytes a list that has pathways ramp Id as name, analytes (ramp compound id only
#' or gene id only or both) as content. Return by findAnalyteHasPathway
#' @param methods must be in c('balanced','weighted') to determine which way to calculate this matrix
#' @return the overlap matrix that has the overlap
compute_overlap_matrix <- function(pathwayid,
                                   pathwaysWithAnalytes,
                                   methods){
  if(!(methods %in% c('balanced','weighted')))
    stop('Wrong option for the input')
  analyte_result <- matrix(NA,nrow = length(pathwayid),ncol = length(pathwayid))
  colnames(analyte_result) <- pathwayid
  rownames(analyte_result) <- pathwayid
  # First method compute intersection over the union
  if(methods == 'balanced'){
    for(i in 1:length(pathwayid)){
      id <- pathwayid[i]
      cid <- pathwaysWithAnalytes[[i]]
      for (j in 1:length(pathwayid)) {
        if(is.na(analyte_result[i,j])){
          if(i==j){
            analyte_result[i,j] <- 1
          }else{
            cid2 <- pathwaysWithAnalytes[[j]]
            shared_metabolite <- unique(intersect(cid,cid2))
            total <- unique(union(cid,cid2))
            analyte_result[i,j] <- length(shared_metabolite)/length(total)
            print(analyte_result[i,j])
            if(is.na(analyte_result[j,i])){
              analyte_result[j,i] <- analyte_result[i,j]
            }
          }
        }
        print(paste("Compute for ",i,",",j))
      }
    }
  }else if (methods == 'weighted'){
    # second method
    for(i in 1:length(pathwayid)){
      id <- pathwayid[i]
      cid <- pathwaysWithAnalytes[[i]]
      for (j in 1:length(pathwayid)) {
        if(is.na(analyte_result[i,j])){
          if(i==j){
            analyte_result[i,j] <- 1
          }else{
            cid2 <- pathwaysWithAnalytes[[j]]
            shared_metabolite <- unique(intersect(cid,cid2))
            total <- unique(union(cid,cid2))
            analyte_result[i,j] <- length(shared_metabolite)/length(unique(cid2))
            print(analyte_result[i,j])
            if(is.na(analyte_result[j,i])){
              analyte_result[j,i] <- length(shared_metabolite)/length(unique(cid))
            }
          }
        }
        print(paste("Compute for ",i,",",j))
      }
    }
  }

  return(analyte_result)
}
#' Update the overlap matrix for store in the shiny app directory
#'
#' @param min_analyte a int that specifies the minimum of analytes the
#' pathway should have to be considered compute for overlap matrix
#' @param method a string that specifies the way to compute overlap matrix,
#' must be 'balanced' or 'weighted'
#' @param together a boolean value to compute overlap matrix for
#' gene/metabolites separatly or together
#' @param username a string that specifies name of MySQL database
#' @param dbname a string that specifies database name of MySQL database
#' @param conpass a string that specifies password for database connection
#' @param host a string that specifes host for database connection
updateOverlapMatrix <- function(min_analyte,method,together,conpass,
                                host = 'localhost',dbname = 'ramp',
                                username = 'root'){
  if(!together){
    con <- DBI::dbConnect(RMySQL::MySQL(),
                     user = username,
                     dbname= dbname,
                     password = conpass,
                     host = host)
    on.exit(DBI::dbDisconnect(con))
    pathways<- DBI::dbGetQuery(con,'select * from pathway;')
    source <- DBI::dbGetQuery(con,'select * from source;')


    dbname <- unique(pathways$type)

    # pathwayInHmdb <- pathways[pathways$type == 'hmdb',]
    pathwayInKegg <- pathways[pathways$type == 'kegg',]
    pathwayInWiki <- pathways[pathways$type == 'wiki',]
    pathwayInReac <- pathways[pathways$type == 'reactome',]

    # define the minimum metabolites/genes


    # Store Compound Ids in List
    # listOfHmdbC <- findAnalyteHasPathway(pathwayInHmdb$pathwayRampId)
    listOfKeggC <- RaMP:::findAnalyteHasPathway(pathwayInKegg$pathwayRampId,n = min_analyte,host = host,
                                         conpass = conpass, dbname = dbname,host = host)
                                         
    listOfWikiC <- RaMP:::findAnalyteHasPathway(pathwayInWiki$pathwayRampId,n = min_analyte,host = host,
                                         conpass = conpass, dbname = dbname,host = host)
    listOfReacC <- RaMP:::findAnalyteHasPathway(pathwayInReac$pathwayRampId,n = min_analyte,host = host,
                                         conpass = conpass, dbname = dbname,host = host)
    # Store Gene Ids in List
    # listOfHmdbG <- findAnalyteHasPathway(pathwayInHmdb$pathwayRampId,GC="G")
    listOfKeggG <- RaMP:::findAnalyteHasPathway(pathwayInKegg$pathwayRampId,GC="G",n = min_analyte,host = host,
                                         conpass = conpass, dbname = dbname,host = host)
    listOfWikiG <- RaMP:::findAnalyteHasPathway(pathwayInWiki$pathwayRampId,GC="G",n = min_analyte,host = host,
                                         conpass = conpass, dbname = dbname,host = host)
    listOfReacG <- RaMP:::findAnalyteHasPathway(pathwayInReac$pathwayRampId,GC="G",n = min_analyte,host = host,
                                         conpass = conpass, dbname = dbname,host = host)
    # Setup minimum number of analytes that will be considered

    # May need to filter out that pathway that has less than 5 metabolites
    # Output to a matrix
    # In order of HMDB Kegg Wiki Reac

    pathwayid <- c(#names(listOfHmdbC),
      names(listOfKeggC),
      names(listOfWikiC),
      names(listOfReacC))

    pathToanalC <- do.call(c,list(#listOfHmdbC,
      listOfKeggC,
      listOfWikiC,
      listOfReacC))
    if(methods == 'balanced' ){
      metabolite_result <- RaMP:::compute_overlap_matrix(pathwayid = pathwayid,
                                                  pathwaysWithAnalytes =  pathToanalC,
                                                  methods = 'balanced')
    }else if(methods == 'weighted'){
      metabolite_result <- RaMP:::compute_overlap_matrix(pathwayid = pathwayid,
                                                   pathwaysWithAnalytes = pathToanalC,
                                                   methods = 'weighted')
    }

    # Output to a matrix
    ### Part for genes ...
    pathwayidG <- c(#names(listOfHmdbG),
      names(listOfKeggG),
      names(listOfWikiG),
      names(listOfReacG))



    pathToanalG <- do.call(c,list(#listOfHmdbG,
      listOfKeggG,
      listOfWikiG,
      listOfReacG))
    # compute for matrix
    if(methods == 'balanced' ){
      gene_result <- RaMP:::compute_overlap_matrix(pathwayid = pathwayidG,pathToanalG,methods = 'balanced')
    } else if(methods == 'weighted' ){
      gene_result <- RaMP:::compute_overlap_matrix(pathwayid = pathwayidG,pathToanalG,methods = 'weighted')
    }
    return(list(
      metabolite = metabolite_result,
      gene = gene_result
    ))
  } else if(together){
    con <- DBI::dbConnect(MySQL(),
                     user = username,
                     dbname=dbname,
                     password = conpass,
                     host = host)

    pathways<- DBI::dbGetQuery(con,'select * from pathway;')



    dbname <- unique(pathways$type)

    # pathwayInHmdb <- pathways[pathways$type == 'hmdb',]
    pathwayInKegg <- pathways[pathways$type == 'kegg',]
    pathwayInWiki <- pathways[pathways$type == 'wiki',]
    pathwayInReac <- pathways[pathways$type == 'reactome',]

    # double minimum number of analytes in the pathway
    min_analyte <- min_analyte * 2

    # Store Compound Ids in List
    # listOfHmdbC <- findAnalyteHasPathway(pathwayInHmdb$pathwayRampId)
    # use both to save metabolites/genes in the list
    listOfKegg <- RaMP:::findAnalyteHasPathway(pathwayInKegg$pathwayRampId,GC = 'both',n = min_analyte)
    listOfWiki <- RaMP:::findAnalyteHasPathway(pathwayInWiki$pathwayRampId,GC = 'both',n = min_analyte)
    listOfReac <- RaMP:::findAnalyteHasPathway(pathwayInReac$pathwayRampId,GC = 'both',n = min_analyte)
    # Append all pathways id together

    pathwayid <- c(#names(listOfHmdbC),
      names(listOfKegg),
      names(listOfWiki),
      names(listOfReac))

    pathToanal <- do.call(c,list(#listOfHmdbC,
      listOfKegg,
      listOfWiki,
      listOfReac))
    if(methods == 'balanced'){
      analyte_result <- RaMP:::compute_overlap_matrix(pathwayid = pathwayid,
                                               pathwaysWithAnalytes =  pathToanal,
                                               methods = 'balanced')
    } else if(methods == 'weighted'){
      analyte_result <- RaMP:::compute_overlap_matrix(pathwayid = pathwayid,
                                               pathwaysWithAnalytes =  pathToanal,
                                               methods = 'weighted')
    }
    return(analyte_result)
  }
}
