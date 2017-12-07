#' Generate HTML table output from given data frame
#' 
#' The HTML output is accepted by shiny HTML output, and 
#' form well-format table on webpage
#' @param data The data converted to table format (data.frame)
#' @param num_row The rows of table converted from raw data (integer)
#' @return A string that contains all html output
#' @export
rampTablize <- function(data,num_row =5){
  if(is.character(data))
    return(NULL)
  tables <- list()
  tb_title <- unique(data[,ncol(data)])
  for(item in tb_title){
    table <- data[data[,ncol(data)] %in% item,]
    colnames(table)[1]<-paste(nrow(table),"items from",item)
    if (nrow(table)>num_row){
      tables[[item]] <- print(xtable::xtable(table[1:num_row,]),
                              type = "html",
                              html.table.attributes ="class = 'data table table-bordered table-condensed'",
                              caption.placement = "top",
                              include.rownames = F)
      
  
    } else {
      tables[[item]] <- print(xtable::xtable(table[1:nrow(table),]),
                              type = "html",
                              html.table.attributes ="class = 'data table table-bordered table-condensed'",
                              caption.placement = "top",
                              include.rownames = F)
    }
  }
  return(lapply(tables,paste))
  
}
#' Send query to databases to get metabolites 
#'
#' From user supplied synonym metabolites, it search through whole databases 
#' , and returns all metabolites which are in the same pathways.
#'
#' @param names the string that user defined for given synonym of metablite
#' @param maxItems the integer given from slider bar to limit maximum elements
#' returned by query
#' @param geneOrcompound string defined from radio input that describes type
#' of returned metablites from query.
#' @examples
#' \dontrun{
#' ramp <- rampGenesFromComp("VitamineE",20,"compound")
#' @return a dataframe that contain given synonym as column name}
#' @export
rampGenesFromComp <- function(names, maxItems, geneOrcompound = NULL) {
    # con <- dbConnect(MySQL(), user = "root", password = "Ramp340!", dbname = "mathelabramp")
    result <- DBI::dbGetQuery(con, paste0("SELECT analytesynonym.Synonym,analytesynonym.geneOrCompound,source.sourceId,source.IDtype
                            FROM analytesynonym,source WHERE analytesynonym.rampID IN
                                     (SELECT rampId FROM analytehaspathway WHERE pathwayRampId IN
                                     (SELECT pathwayRampId FROM analytehaspathway WHERE rampId IN
                                     (SELECT rampID FROM analytesynonym WHERE synonym = \"",
                                     names, "\"))) AND source.rampId = analytesynonym.rampId LIMIT ", maxItems, ";"))

    colnames(result) <- c(paste0("Search result for ", names),"metabolites_type","sourceID","source_type")
    # dbDisconnect(con)
    return(result)
}

#' Send query to databases to get metabolites
#' 
#' From given pathway name, it searches through whole databases to find
#' all metabolites in that pathway.
#'
#' @param names the string value that describes pathway name.
#' @param maxItems the integer value given to define maximum items
#' return by this query.
#' @param geneOrCompound the string value that define types (gene or
#' compound) of returned metabolites.
#' @return a dataframe that has given pathwayname as column names
#' ramp <- rampNameFromPath("Glucose-Alanine Cycle",20,compound)
#' @export
rampNameFromPath <- function(names, maxItems, geneOrCompound) {
  # con <- dbConnect(MySQL(), user = "root", password = "Ramp340!", dbname = "mathelabramp")
  if (geneOrCompound != "both"){
      result <- DBI::dbGetQuery(con, paste0("select analytesynonym.Synonym,analytesynonym.geneOrCompound,source.sourceId,source.IDtype 
                                        from analytesynonym,source where analytesynonym.rampId in 
                                       (select rampId from analytehaspathway where pathwayRampId in 
                                       (select pathwayRampId from pathway where pathwayName = \"",
                                        names,"\")) and source.rampId = analytesynonym.rampId and analytesynonym.geneOrCompound = \"",
                                        geneOrCompound,"\" limit ", maxItems,";"))
  } else {
      result <- DBI::dbGetQuery(con, paste0("select analytesynonym.Synonym,analytesynonym.geneOrCompound,source.sourceId,source.IDtype 
                                        from analytesynonym,source where analytesynonym.rampId in 
                                       (select rampId from analytehaspathway where pathwayRampId in 
                                       (select pathwayRampId from pathway where pathwayName = \"",
                                       names,"\")) and source.rampId = analytesynonym.rampId limit ", maxItems,";"))
  }
  if (nrow(result)>0){
    colnames(result) <- c(paste0("Result_for_", names),"metabolites_type","source_ID", "source_type")
  } else {
    return ("Not Found")
  }
  # dbDisconnect(con)
  return(result)
}

#' Send query to database that find pathway
#' 
#' From given metabolite's synonym, it search through all databases to find
#' all pathways that have the metabolite.
#'
#' @param synonym the string that describe given synonym name.
#' @param maxItems the integer that limit maximum returned by
#' query
#' @return a dataframe that has given synonym as column name
#' @export
rampPathFromMeta <- function(synonym, maxItems) {
    # con <- dbConnect(MySQL(), user = "root", password = "Ramp340!", dbname = "mathelabramp")
    sql <- paste0("select pathwayRampId from analytehaspathway where rampId in
                  (select rampID from analytesynonym where Synonym = \"",synonym, "\");")
    cat(sql)
    query <- DBI::dbGetQuery(con,sql)
    container <- data.frame(matrix(nrow = 0,ncol = 3))
    colnames(container) <- c(paste0("Search result for ", synonym), "source_id", "source_type")
    for(id in unlist(query)){
      sql <- paste0("select pathwayName,sourceId,type from pathway where pathwayRampId =\"",id,"\";")
      query2 <- DBI::dbGetQuery(con,sql)
      colnames(query2) <- c(paste0("Search result for ", synonym), "source_id", "source_type")
      container <- rbind(container,query2)
    }
    
    # dbDisconnect(con)

    return(container)
}

#' The function does key word search from given synonym on within databases
#'
#' @param word string value that describe the keyword user want to search
#' @param database string value that describe where the user want to search
#' for key word
#' @return If there is at least one itmes in database vaguely matching key
#' word, it will return a data frame that contains all search result. Otherwise
#' it will return a string value to inform user.
#' @export
rampKWsearch <- function(word, database) {
    # con <- dbConnect(MySQL(), user = "root", password = "Ramp340!", dbname = "mathelabramp")
    item <- "*"
    if (database == "analytesynonym") {
        item <- "Synonym"
    } else if (database == "pathway") {
        item <- "pathwayName"
    } else if (database == "ontology") {
        item <- "commonName"
    } else {
        return("Not Found")
    }
  
    query <- paste0("select distinct ", item, " from ", database,
                    " where (", item, " like \"%", 
                    word, "%\") order by ",
                    "char_length(",
                    item,"), ",
                    item," ASC limit 20;")

    rampOut <- DBI::dbGetQuery(con, query)
    # dbDisconnect(con)
    if (length(unlist(rampOut)) == 0)
        rampOut <- "Not Found"
    return(rampOut)
}

#' The function send query to databases to find synonym 
#' 
#' based on catalyzation relationship. 
#' This function also automate identification about if the given
#' metabolites is gene or compound
#'
#' @param synonym string value that describe synonym of metabolites
#' @param maxItems Integer value that describe maximum returned by
#' query
#' @return a dataframe that has all searching result with given synonym
#' as column name. If \code{synonym} is compound, it returns all gene that
#' catalyze that compound. If \code{synonym} is gene, it returns all compound
#' catalyzed by this gene.
#' @export
rampCataOut <- function(synonym, maxItems = 1000) {
    # con <- dbConnect(MySQL(), user = "root", password = "Ramp340!", dbname = "mathelabramp")
    rampOut <- DBI::dbGetQuery(con, paste0("select synonym from analytesynonym where rampId in
                        (select if((select geneOrCompound from analyteSynonym where synonym = \"",
        synonym, "\" limit 1) =\"compound\",rampGeneId,rampCompoundId) from catalyzed where
                        if((select geneOrCompound from analyteSynonym where synonym = \"",
        synonym, "\" limit 1) =
                        \"compound\",rampCompoundId,rampGeneId) in
                        (select rampID from analytesynonym where synonym = \"",
        synonym, "\")) limit ", maxItems, ";"))
    # dbDisconnect(con)
    return(rampOut)
}

#' The function send query to databases to find ontology or metabolites 
#' 
#' Based on ontology relationship, it search throught databases to find 
#' metabolites or ontology from given name which could be either ontology 
#' location or metabolite's synonym.
#'
#' @param synonym string value that describe given metabolites or ontology
#' @param maxItems integer value that describe maximum items returned by query
#' @return a dataframe that has given synonym as column name.
#' @export
rampOntoOut <- function(synonym, maxItems) {
    con <- dbConnect(MySQL(), user = "root", password = "Ramp340!", dbname = "mathelabramp")
    out <- DBI::dbGetQuery(con, paste0("select * from analyteSynonym where Synonym = \"", synonym, "\";"))
    if (nrow(out) > 0) {
        rampOut <- DBI::dbGetQuery(con, paste0("select commonName,biofluidORcellular from ontology where rampOntologyIdLocation in
      (select rampOntologyidlocation from analyteHasOntology where rampCompoundId in
       (select rampID from analyteSynonym where synonym = \"",
            synonym, "\")) limit ", maxItems, ";"))

    } else {
        rampOut <- DBI::dbGetQuery(con, paste0("select synonym,geneOrCompound from analyteSynonym where rampId in
      (select rampCompoundId from analyteHasOntology where rampOntologyIdLocation in
       (select rampOntologyIdLocation from ontology where commonName=\"",
            synonym, "\")) limit ", maxItems, ";"))
    }
    colnames(rampOut) <- c(paste("Search result for", synonym), "type")

    # dbDisconnect(con)
    return(rampOut)
}

#' Kill all databases connection.
#' 
#' If the connections is forgotten to be disconnected, it will cause potential trouble
#' for further connection.
#' 
#' @export
killDbConnections <- function() {
    all_cons <- dbListConnections(MySQL())

    print(all_cons)

    for (con in all_cons) +dbDisconnect(con)
    print(paste(length(all_cons), " connections killed."))
}

