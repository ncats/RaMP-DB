
#' Use fast search algorithm to find all pathways from given analytes
#'
#' @param pathway a string or a vector of strings that contains pathways of interest
#' @param analyteType a string denoting the type of analyte to return ("gene", "metabolite", "both")
#' @param match type of matching to use, options are "exact" or "fuzzy".  The default is "exact".
#' @param maxPathwaySize (default Inf), trims returned results to pathways that have fewer than this number
#' @param namesOrIds are the input pathways input as pathway names or as pathway ids
#' of genes and metabolites
#' @param db a RaMP database object, if not specified a new one is created with RaMP::RaMP()
#' @param ... Internal Use - for handling deprecated parameter names
#' @return a data.frame that contains all search results
#' @examples
#' \dontrun{
#'
#' getAnalyteFromPathway( pathway=c("Wnt Signaling Pathway", "Sphingolipid Metabolism"), db = rampDB )
#'
#' }
#' @export
getAnalyteFromPathway <- function(pathway, match="exact", analyteType="both", maxPathwaySize = Inf,
                                  namesOrIds="names", db = RaMP(), ...) {
  analyteType <- handleRenamedParameter(argument = analyteType, oldName = "analyte_type", version = "3.0")
  maxPathwaySize <- handleRenamedParameter(argument = maxPathwaySize, oldName = "max_pathway_size", version = "3.0")
  namesOrIds <- handleRenamedParameter(argument = namesOrIds, oldName = "names_or_ids", version = "3.0")
  assertDBparamIsRight(firstParam = pathway, dbParam = db)

  now <- proc.time()
  print("fired!")

  # Retrieve pathway RaMP ids
  df <- db@api$getAnalytesFromPathways(pathways = pathway, namesOrIds = namesOrIds, match = match)

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
