#' run shiny app
#' 
#' if has database locally, user must provide database name
#' username password
#' @param dbname a string that is database name including all tables
#' @param username a string that is username for database
#' @param password a string this is password for database
#' @param host a string that stand for host 
#' @param update a logic value determines if RaMP update local data
#' @export
runRaMPapp <- function(dbname = "mathelabramp",
                       username = "root",
                       password = "Ehe131224",
                       host = NULL,
                       update = FALSE){
  
  con <<- DBI::dbConnect(
    drv = RMySQL::MySQL(),
    dbname = dbname,
    username = username,
    password = password
  )
  appDir <- system.file("shinyApp",package = "RaMP")
  if (appDir == ""){
    dbDisconnect(con)
    stop("Could not find example directory. Try re-installing 'ramp'.")
  }
  #if(update){
    # message("Updating local data usually takes 10 mins ...")
    # query <- "select pathwayName from pathway;"
    # 
    # pathway_list <- dbGetQuery(con,query)
    # 
    # pathway_list <- unlist(pathway_list)
    # tot_metabolites <- dbGetQuery(con,"select count(*) from analytesynonym;")
    # tot_metabolites <- tot_metabolites[1,1]
    # tb <- data.frame(inpathway=NULL,outpathway=NULL)
    # 
    # for (pathway in pathway_list){
    #   tot_in_pathway <- dbGetQuery(con,paste0("select count(*) from analyte 
    #                                           where rampId in (select rampId from 
    #                                           analytehaspathway where 
    #                                           pathwayRampId in (select pathwayRampId 
    #                                           from pathway where pathwayName = \"",
    #                                           pathway,"\"));"))
    #   tot_out_pathways <- tot_metabolites - tot_in_pathway
    #   item <- data.frame(inpathway = tot_in_pathway,outpathway = tot_out_pathways)
    #   row.names(item) <- pathway
    #   tb <- rbind(tb,item)
    # }
    # colnames(tb) <- c("inPathway","outPathway")
    # write.csv(tb,paste0(appDir,"/dataFisherTest.csv"))
    # query1 <- "select synonym from analytesynonym;"
    # query2 <- "select pathwayName from pathway;"
    # query3 <- "select commonName from ontology;"
    # 
    # analyte <- dbGetQuery(con,query1)
    # pathway <- dbGetQuery(con,query2)
    # ontology <- dbGetQuery(con,query3)
    # 
    # write.csv(analyte,paste0(appDir,"/analyte.csv"),row.names = F)
    # write.csv(pathway,paste0(appDir,"/pathway.csv"),row.names = F)
    # write.csv(ontology,paste0(appDir,"/biofluid.csv"),row.names = F)
  #}
  
 
  shiny::runApp(appDir,display.mode = "normal")
  
}
