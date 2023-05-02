
#' Use fast search algorithm to find all pathways from given analytes
#'
#' @param pathway a string or a vector of strings that contains pathways of interest
#' @param analyte_type a string denoting the type of analyte to return ("gene", "metabolite", "both")
#' @param match type of matching to use, options are "exact" or "fuzzy".  The default is "exact".
#' @param max_pathway_size (default Inf), trims returned results to pathways that have fewer than this number
#' @param names_or_ids are the input pathways input as pathway names or as pathway ids
#' of genes and metabolites
#' @return a data.frame that contains all search results
#' @examples
#' \dontrun{
#' # To query one pathway:
#' myanalytes <- getAnalyteFromPathway2(pathway="sphingolipid metabolism")
#'
#' # To query multiple pathways:
#' pkg.globals <- setConnectionToRaMP(dbname="ramp2",username="root",conpass="",host = "localhost")
#' myanalytes <- getAnalyteFromPathway2(pathway=c("De Novo Triacylglycerol Biosynthesis",
#'	"sphingolipid metabolism"))
#' }
#' @export
getAnalyteFromPathway <- function(pathway, match="exact", analyte_type="both", max_pathway_size = Inf, names_or_ids="names") {
  now <- proc.time()
  print("fired!")
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

  pathwayMatchCol = 'pathwayName'
  if(names_or_ids == 'ids') {
    pathwayMatchCol = 'sourceId'
    match = 'exact'
  }

  # Retrieve pathway RaMP ids
  if (match=='exact') {

    # return pathway name, pathway type, analyte name, source analyte ids, analyte type/class
    sql = paste0("select
    group_concat(distinct s.commonName order by s.commonName asc separator '; ') as analyteName,
    group_concat(distinct s.sourceId order by s.sourceId asc separator '; ') as sourceAnalyteIDs,
    s.geneOrCompound as geneOrCompound,
    p.pathwayName as pathwayName,
    p.sourceId as pathwayId,
    p.pathwayCategory as pathwayCategory,
    p.type as pathwayType
    from pathway p, analytehaspathway ap, source s
    where s.rampId = ap.rampID
    and ap.pathwayRampId = p.pathwayRampId
    and (p.pathwayCategory not like 'smpdb%' or p.pathwayCategory is Null)
    and p.",pathwayMatchCol," in (",list_pathway,") ",
    "group by s.rampId, p.pathwayName, p.sourceId, p.type, s.geneOrCompound
    order by p.type desc, p.pathwayName asc, s.geneOrCompound asc;")
    con <- connectToRaMP()
    df <- RMariaDB::dbGetQuery(con,sql)
    RMariaDB::dbDisconnect(con)
  } else if(match == 'fuzzy') {
    df = data.frame(matrix(nrow=0, ncol=6))
    sql = paste0("select
    group_concat(distinct s.commonName order by s.commonName asc separator '; ') as analyteName,
    group_concat(distinct s.sourceId order by s.sourceId asc separator '; ') as sourceAnalyteIDs,
    s.geneOrCompound as geneOrCompound,
    p.pathwayName as pathwayName,
    p.sourceId as pathwayId,
    p.pathwayCategory as pathwayCategory,
    p.type as pathwayType
    from pathway p, analytehaspathway ap, source s
    where s.rampId = ap.rampID
    and ap.pathwayRampId = p.pathwayRampId
    and (p.pathwayCategory not like 'smpdb%' or p.pathwayCategory is Null)
    and p.",pathwayMatchCol," like '%[SOME_PW_NAME]%' group by s.rampId, p.pathwayName, p.sourceId, p.type, s.geneOrCompound
    order by p.type desc, p.pathwayName asc, s.geneOrCompound asc;")

    con <- connectToRaMP()
    for(p in pathway) {
      currSQL = gsub(pattern = '[SOME_PW_NAME]', replacement = p, x= sql, fixed = T )
      subdf <- RMariaDB::dbGetQuery(con,currSQL)
      df <- rbind(df, subdf)
    }
    RMariaDB::dbDisconnect(con)
  }

  # if we have a result and max_pathway size is not Infinite, filter pathway results by pathway size
  if(nrow(df) > 0 && max_pathway_size != Inf) {
    pwAnalyteCounts <- data.frame(table(df$`pathwayName`))
    pwAnalyteCounts <- pwAnalyteCounts[pwAnalyteCounts$Freq <= max_pathway_size,]
    df <- df[df$`pathwayName` %in% unlist(pwAnalyteCounts$Var1),]
  }

  if(analyte_type=="gene") {
    print("gene return...")
    allout <- df[which(df$`geneOrCompound`=="gene"),]
  } else if (analyte_type=="metabolite") {
    print("met return...")
    allout <- df[which(df$`geneOrCompound`=="compound"),]
  } else {
    allout <- df
  }

  print("Timing ..")
  print(proc.time() - now)

  return(allout)
}


# (Original Implementation)
# getAnalyteFromPathway <- function(pathway, match="exact", analyte_type="both") {
#   now <- proc.time()
#   print("fired")
#   if(is.character(pathway)){
#     if(grepl("\n",pathway)[1]){
#       list_pathway <- strsplit(pathway,"\n")
#       list_pathway <- unlist(list_pathway)
#     } else if(grepl(",",pathway)[1]){
#       list_pathway <- strsplit(pathway,"\n")
#       list_pathway <- unlist(list_pathway)
#     } else {
#       list_pathway <- pathway
#     }
#   } else if(is.data.frame(pathway)){
#     list_pathway <- unlist(pathway)
#   } else {
#     return("Wrong Data Format")
#   }
#   list_pathway <- sapply(list_pathway,shQuote)
#   list_pathway <- paste(list_pathway,collapse = ",")
#
#   # Retrieve pathway RaMP ids
#   if (match=='exact') {
#   	query1 <- paste0("select * from pathway where pathwayName
#                    in (",list_pathway,");")
# 	con <- connectToRaMP()
# 	df1 <- RMariaDB::dbGetQuery(con,query1)
# 	RMariaDB::dbDisconnect(con)
#   } else if (match=='fuzzy') {
# 	print("running fuzzy")
# 	df1=c()
# 	for (i in 1:length(pathway)) {
# 		# note here that we are using pathway, not list_pathway which
# 		# formats for 'exact' but not 'fuzzy'
# 		con <- connectToRaMP()
# 		query1 <- paste0('select * from pathway where pathwayName
#                    like "%',pathway[i],'%";')
# 		df1 <- rbind(df1,RMariaDB::dbGetQuery(con,query1))
# 		RMariaDB::dbDisconnect(con)
# 	}
#   } else {
# 	stop("Please be sure to set the match parameter to 'exact' or 'fuzzy'.")
#   }
#
#   if(nrow(df1)==0) {
#     stop("None of the input pathway(s) could be found")}
#
#   # Retrieve compound id from RaMP pathway id (query1)
#   #query2 <- paste0("select pathwayRampId,rampId from analytehaspathway where
#   #                 pathwayRampId in (select pathwayRampId from pathway where
#   #                 pathwayName in (",list_pathway,"));")
#   pidlist <- sapply(df1$pathwayRampId,shQuote)
#   pidlist <- paste(pidlist,collapse = ",")
#
#   query2 <- paste0("select pathwayRampId, rampId from analytehaspathway
# 	where pathwayRampId in (",pidlist,");")
#   con <- connectToRaMP()
#   df2 <- RMariaDB::dbGetQuery(con,query2)
#   RMariaDB::dbDisconnect(con)
#   cid_list <- unlist(df2[,2])
#   cid_list <- sapply(cid_list,shQuote)
#   cid_list <- paste(cid_list,collapse = ",")
#
#   # Retrieve all common name from compounds associated with RaMP compound ids (query2)
#   query3 <- paste0("select * from source where rampId in (",cid_list,");")
#   con <- connectToRaMP()
#   df3 <- RMariaDB::dbGetQuery(con,query3)
#   RMariaDB::dbDisconnect(con)
#
#   # Merge all of this together
#   mdf1 <- merge(df3,df2,all.x=T)
#   mdf1 <- merge(mdf1,df1,all.x=T,by="pathwayRampId")
#   mdf1 <- mdf1[!duplicated(mdf1),]
#   #mdf1[,-which(colnames(mdf1) %in% "sourceId.y")]
#   colnames(mdf1)[which(colnames(mdf1)=="sourceId.x")]="sourceId"
#   mdf1$temp<-paste(mdf1$pathwayRampId,mdf1$rampId,sep=";")
#   # use regular expression to remove substring before ':'
#   #mdf1$sourceId <- gsub('.*:','',mdf1$sourceId)
#   out=data.frame(pathwayName=NA,pathwayType=NA,analyteName=NA,
#                  sourceAnalyteIDs=NA,geneOrCompound=NA)
#   # Reformat so that you have one metabolite, with all synonyms, in one line:
#   count=1
#   allout=c()
#   for (i in unique(mdf1$temp)) {
#     temp <- mdf1[which(mdf1$temp==i),]
#     #out$sourceCompoundIDs <- paste(paste(temp$IDtype,temp$sourceId,sep=": "),collapse="; ")
#     out$sourceAnalyteIDs <- paste(unique(temp$sourceId),collapse="; ")
#     out$pathwayName <- temp[1,"pathwayName"]
#     #out$pathwayCategory <- temp[1,"pathwayCategory"]
#     out$pathwayType <- temp[1,"type"]
#     out$analyteName <- temp[1,"commonName"]
#     out$geneOrCompound <- temp[1,"geneOrCompound"]
#     count=count+1
#     allout=rbind(allout,out)
#   }
#
#   print("Timing ..")
#   print(proc.time() - now)
#
#   if(analyte_type=="gene") {
# 	allout <- allout[which(allout$geneOrCompound=="gene"),]
#   } else if (analyte_type=="metabolite") {
#         allout <- allout[which(allout$geneOrCompound=="compound"),]
#   } else {
#     allout <- allout
#   }
#  return(allout)
# }
