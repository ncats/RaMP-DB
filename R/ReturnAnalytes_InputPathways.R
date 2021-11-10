#' Use fast search algorithm to find all pathways from given analytes
#'
#' @param pathway a string or a vector of strings that contains pathways of interest
#' @return a data.frame that contains all search results
#' @examples
#' \dontrun{
#' # To query one pathway:
#' myanalytes <- getAnalyteFromPathway(pathway="sphingolipid metabolism")
#'
#' # To query multiple pathways:
#' pkg.globals <- setConnectionToRaMP(dbname="ramp2",username="root",conpass="",host = "localhost")
#' myanalytes <- getAnalyteFromPathway(pathway=c("De Novo Triacylglycerol Biosynthesis", 
#'	"sphingolipid metabolism"))
#' }
#' @export
getAnalyteFromPathway <- function(pathway) {
  
  now <- proc.time()
  print("fired")
  if(is.character(pathway)){
    if(grepl("\n",pathway)[1]){
      list_pathway <- strsplit(pathway,"\n")
      list_pathway <- unlist(list_pathway)
    } else if(grepl(",",pathway)[1]){
      list_pathway <- strsplit(pathway,"\n")
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
  # Retrieve pathway RaMP id
  con <- connectToRaMP()
  query1 <- paste0("select * from pathway where pathwayName
                   in (",list_pathway,");")
  
  df1 <- DBI::dbGetQuery(con,query1)
  DBI::dbDisconnect(con)
  
  if(nrow(df1)==0) {
    stop("None of the input pathway(s) could be found")}
  
  # Retrieve compound id from RaMP pathway id (query1)
  query2 <- paste0("select pathwayRampId,rampId from analytehaspathway where
                   pathwayRampId in (select pathwayRampId from pathway where
                   pathwayName in (",list_pathway,"));")
  con <- connectToRaMP()
  df2 <- DBI::dbGetQuery(con,query2)
  DBI::dbDisconnect(con)
  cid_list <- unlist(df2[,2])
  cid_list <- sapply(cid_list,shQuote)
  cid_list <- paste(cid_list,collapse = ",")
  
  # Retrieve all common name from compounds associated with RaMP compound ids (query2)
  query3 <- paste0("select * from source where rampId in (",cid_list,");")
  con <- connectToRaMP()
  df3 <- DBI::dbGetQuery(con,query3)
  DBI::dbDisconnect(con)
  
  # Merge all of this together
  mdf1 <- merge(df3,df2,all.x=T)
  mdf1 <- merge(mdf1,df1,all.x=T,by="pathwayRampId")
  mdf1 <- mdf1[!duplicated(mdf1),]
  #mdf1[,-which(colnames(mdf1) %in% "sourceId.y")]
  colnames(mdf1)[which(colnames(mdf1)=="sourceId.x")]="sourceId"
  mdf1$temp<-paste(mdf1$pathwayRampId,mdf1$rampId,sep=";")
  # use regular expression to remove substring before ':'
  #mdf1$sourceId <- gsub('.*:','',mdf1$sourceId)
  out=data.frame(pathwayName=NA,pathwayCategory=NA,pathwayType=NA,compoundName=NA,
                 sourceCompoundIDs=NA,geneOrCompound=NA)
  # Reformat so that you have one metabolite, with all synonyms, in one line:
  count=1
  allout=c()
  for (i in unique(mdf1$temp)) {
    temp <- mdf1[which(mdf1$temp==i),]
    out$sourceCompoundIDs <- paste(paste(temp$IDtype,temp$sourceId,sep=": "),collapse="; ")
    out$pathwayName <- temp[1,"pathwayName"]
    #out$pathwayCategory <- temp[1,"pathwayCategory"]
    out$pathwayType <- temp[1,"type"]
    out$compoundName <- temp[1,"commonName"]
    out$geneOrCompound <- temp[1,"geneOrCompound"]
    count=count+1
    allout=rbind(allout,out)
  }
  
  print("Timing ..")
  print(proc.time() - now)
  
  return(as.data.frame(allout))
}

