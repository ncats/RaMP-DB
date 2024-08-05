
#' Use fast search algorithm to find all pathways from given analytes
#'
#' @param pathway a string or a vector of strings that contains pathways of interest
#' @param analyteType a string denoting the type of analyte to return ("gene", "metabolite", "both")
#' @param match type of matching to use, options are "exact" or "fuzzy".  The default is "exact".
#' @param maxPathwaySize (default Inf), trims returned results to pathways that have fewer than this number
#' @param namesOrIds are the input pathways input as pathway names or as pathway ids
#' of genes and metabolites
#' @param db a RaMP database object, if not specified a new one is created with RaMP::RaMP()
#' @return a data.frame that contains all search results
#' @examples
#' \dontrun{
#'
#' getAnalyteFromPathway( pathway=c("Wnt Signaling Pathway", "Sphingolipid Metabolism"), db = rampDB )
#'
#' }
#' @export
getAnalyteFromPathway <- function(pathway, match="exact", analyteType="both", maxPathwaySize = Inf, namesOrIds="names", db = RaMP() ) {
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
  if(namesOrIds == 'ids') {
    pathwayMatchCol = 'sourceId'
    match = 'exact'
  }

  isSQLite = .is_sqlite(db)

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

    if(isSQLite) {
      sql = paste0("select
    group_concat(distinct s.commonName COLLATE NOCASE) as analyteName,
    group_concat(distinct s.sourceId COLLATE NOCASE) as sourceAnalyteIDs,
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
    }

    df <- RaMP::runQuery(sql, db)

  } else if(match == 'fuzzy') {
    df = data.frame(matrix(nrow=0, ncol=7))
    colnames(df) <- c('analyteName', 'sourceAnalyteIDs', 'geneOrCompound',
                      'pathwayName', 'pathwayId', 'pathwayCategory', 'pathwayType')
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

    if(isSQLite) {
      sql = paste0("select
    group_concat(distinct s.commonName COLLATE NOCASE) as analyteName,
    group_concat(distinct s.sourceId COLLATE NOCASE) as sourceAnalyteIDs,
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
    }

    for(p in pathway) {
      if(nchar(p)>2) {
        currSQL = gsub(pattern = '[SOME_PW_NAME]', replacement = p, x= sql, fixed = T )
        subdf <- RaMP::runQuery(currSQL, db)
        df <- rbind(df, subdf)
      }
    }

  }

  # if we have a result and max_pathway size is not Infinite, filter pathway results by pathway size
  if(nrow(df) > 0 && maxPathwaySize != Inf) {
    pwAnalyteCounts <- data.frame(table(df$`pathwayName`))
    pwAnalyteCounts <- pwAnalyteCounts[pwAnalyteCounts$Freq <= maxPathwaySize,]
    df <- df[df$`pathwayName` %in% unlist(pwAnalyteCounts$Var1),]
  }

  if(analyteType=="gene") {
    print("gene return...")
    allout <- df[which(df$`geneOrCompound`=="gene"),]
  } else if (analyteType=="metabolite") {
    print("met return...")
    allout <- df[which(df$`geneOrCompound`=="compound"),]
  } else {
    allout <- df
  }

  print("Timing ..")
  print(proc.time() - now)

  return(allout)
}
