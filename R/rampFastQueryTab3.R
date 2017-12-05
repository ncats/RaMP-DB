#' Do fisher test for only one pathway from search result
#' clicked on highchart
#' @param pathway_list a vector of all pathways searched from given metabolites
#' @param pathway single pathway name that is applied by Fisher test
#' @param num_of_meta integer that describe how many initial metabolites
#' @export
rampOneFisherTest <- function(pathway_list,pathway,num_of_meta){
  # con <- dbConnect(MySQL(), user = "root", password = "Ramp340!", dbname = "mathelabramp")
  # on.exit(dbDisconnect(con))
  contingencyTb <- matrix(0,nrow = 2,ncol = 2)
  colnames(contingencyTb) <- c("In Pathway","Not In Pathway")
  rownames(contingencyTb) <- c("All Metabolites","User's Metabolites")
  tot_in_pathway <- DBI::dbGetQuery(con,paste0("select count(*) from analyte 
                                            where rampId in (select rampId from 
                                            analytehaspathway where 
                                            pathwayRampId in (select pathwayRampId 
                                            from pathway where pathwayName = \"",
                                          pathway,"\"));"))
  tot_in_pathway <- tot_in_pathway[[1]]
  
  tot_metabolites <- DBI::dbGetQuery(con,"select count(*) from analyte;")
  tot_metabolites <- tot_metabolites[[1]]
  tot_out_pathway <- tot_metabolites - tot_in_pathway
  user_in_pathway <- nrow(pathway_list[[pathway]])
  user_out_pathway <- num_of_meta - user_in_pathway
  contingencyTb[1,1] <- tot_in_pathway[[1]]
  contingencyTb[1,2] <- tot_out_pathway[[1]]
  contingencyTb[2,1] <- user_in_pathway
  contingencyTb[2,2] <- user_out_pathway
  
  result <- fisher.test(contingencyTb)
  result$p.value <- round(result$p.value,4)
  return(result)
}
#' Return a list named by pathway
#' 
#' each pathway contains all metabolites inside that pathway
#' @param df A dataframe that has information for bar plot
#' @export
rampGenerateBarPlot <- function(df){
  path_meta_list <- list()
  for (i in 1:nrow(df)){
    if (!is.element(df[i,]$pathway,names(path_meta_list))){
      path_meta_list[[df[i,]$pathway]] <- data.frame(metabolite = df[i,]$metabolite,stringsAsFactors = F)
    } else {
      path_meta_list[[df[i,]$pathway]] <- rbind(path_meta_list[[df[i,]$pathway]],df[i,]$metabolite)
      path_meta_list[[df[i,]$pathway]] <- unique(path_meta_list[[df[i,]$pathway]])
    }
  }
  return(path_meta_list)
}
#' Fisher test for given list of pathways
#' 
#' From user input, the function accept a list of pathways, and number of metabolites 
#' from which the pathways are
#' 
#' @param pathway_list The list contains pathway as names and metabolites in that pathway 
#' under each names (list)
#' @param num_user_metabolites number of metabolites given by user when they want to 
#' search for pathways.
#' @return a data.frame contains all fisher test result with pathway name
#' as column name
#' @export
rampFisherTest <- function(pathway_meta_list,num_user_metabolites,FisherPathwayTable){
  # con <- dbConnect(MySQL(), user = "root", password = "Ehe131224", dbname = "mathelabramp")
  # on.exit(dbDisconnect(con))
  # progress<- shiny::Progress$new()
  # progress$set(message = "Fisher testing ...",value = 0)
  # on.exit(progress$close())
  tot_metabolites <- DBI::dbGetQuery(con,"select count(*) from analyte;")
  tot_metabolites <- tot_metabolites[[1]]
  cumulative_fisher_results <- list()
  
  inc <- length(pathway_meta_list)
  for (pathway in names(pathway_meta_list)){
    # progress$inc(1/inc,"computing ...")
    contingencyTb <- matrix(0,nrow = 2,ncol = 2)
    colnames(contingencyTb) <- c("In Pathway","Not In Pathway")
    rownames(contingencyTb) <- c("All Metabolites","User's Metabolites")
    print(pathway)
    df <- FisherPathwayTable[FisherPathwayTable[,1] == pathway,]
    tot_in_pathway <- df[,2]
    print(tot_in_pathway)
    tot_out_pathway <- df[,3]
    print(tot_out_pathway)
    user_in_pathway <- nrow(pathway_meta_list[[pathway]])
    user_out_pathway <- num_user_metabolites - user_in_pathway
    contingencyTb[1,1] <- tot_in_pathway
    contingencyTb[1,2] <- tot_out_pathway
    contingencyTb[2,1] <- user_in_pathway
    
    contingencyTb[2,2] <- user_out_pathway
    cumulative_fisher_results[[pathway]] <- contingencyTb
  
  }
  
  cumulative_fisher_results <- lapply(cumulative_fisher_results,fisher.test)
  
  return(cumulative_fisher_results)
}
#' Fast search given a list of metabolites
#' @param sysnonym a vector of synonym that need to be searched
#' @return a list contains all metabolits as name and pathway inside.
#' 
#' Apply famil function...
#' 
#' @export
rampFastPathFromMeta<- function(synonym,find_synonym = FALSE){
  # progress<- shiny::Progress$new()
  # progress$set(message = "Querying databases ...",value = 0)
  now <- proc.time()
  # con <- dbConnect(MySQL(), user = "root", password = "Ramp340!", dbname = "mathelabramp")
  # on.exit(dbDisconnect(con))
  # find synonym
  
  synonym <- rampFindSynonymFromSynonym(synonym,find_synonym=find_synonym)
  
  list_metabolite <- unique(synonym)
  list_metabolite <- sapply(list_metabolite,shQuote)
  list_metabolite <- paste(list_metabolite,collapse = ",")
  query1 <- paste0("select distinct Synonym,rampId from analytesynonym where Synonym in (",
                    list_metabolite,");")
  df1<- DBI::dbGetQuery(con,query1)
  query2 <- paste0("select pathwayRampId,rampId from analytehaspathway where 
                    rampId in (select rampId from analytesynonym where Synonym in (",
                    list_metabolite,"));")
  df2 <- DBI::dbGetQuery(con,query2)
  id_list <- df2[,1]
  id_list <- sapply(id_list,shQuote)
  id_list <- paste(id_list,collapse = ",")
  query3 <- paste0("select pathwayName,sourceId,type,pathwayRampId from pathway where pathwayRampId in (",
                    id_list,");")
  df3 <- DBI::dbGetQuery(con,query3)
  mdf <- merge(df1,df2,all.x=T)
  mdf <- mdf[!is.na(mdf[,3]),]
  mdf <- merge(mdf,df3,all.x = T)
  # mdf <- mdf[!is.na(mdf[,3]),]
  colnames(mdf)[3] <- "metabolite"
  mdf <- mdf[,c(4,5,6,3)]
  mdf <- unique(mdf)
  print("timing ...")
  print(proc.time()- now)
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
#' @export
rampFileOfPathways <- function(infile){
  name <- infile[[1,'name']]
    summary <- data.frame(pathway  = character(0),id = character(0),
                          source = character(0),metabolite = character(0))
    rampOut <- list()
    for (i in 1:length(infile[,1])){
      if(infile[[i,'type']]!="text/plain"){
        rampOut[[i]] <- read.csv(infile[[i,'datapath']])
        name <- infile[[i,'name']]
        print(infile[[i,'type']])
        rampOut[[i]]$new.col <- substr(name,1,nchar(name) - 4)
        colnames(rampOut[[i]]) <- c("pathway","id","source","metabolite")
        summary <- rbind(summary,rampOut[[i]])
      } else {
        rampOut <- readLines(infile[[i,'datapath']])
        summary <- rampFastPathFromMeta(rampOut)
      }
    }
    return(summary)
}

#' highchart output for RaMP and fisher test.
#' 
#' Based on given x,y data, type and click event, it returns a highcharter object
#' to highchartOutput function to display bar plot.
#' 
#' @param x_data vector contains data for categorical x-axis
#' @param y_data vector contains frequency of each pathway
#' @param type plot's type of this highcharter object
#' @param event_func Javascript code that define the click event
#' @return highcharter object
#' @export
rampHcOutput <- function(x_data,y_data,type = 'column',event_func){
  fomatterFunc <- JS("function(){
                        html = '<strong> Pathway ' + this.x +' has frequency: ' + this.y +'. metabolites are :'
                        '+this.y.detail+'</strong>;'
                        return html;
                     }")
  hc<-highchart() %>%
    hc_chart(type = type,
             # options3d = list(enabled = TRUE, beta = 15, alpha = 15),
             borderColor = '#ceddff',
             borderRadius = 10,
             borderWidth = 2,
             zoomType = "x",
             backgroundColor = list(
               linearGradient = c(0, 0, 500, 500),
               stops = list(
                 list(0, 'rgb(255, 255, 255)'),
                 list(1, 'rgb(219, 228, 252)')
               ))) %>%
    hc_title(text = "<strong>Search result from given metabolites</strong>",
             margin = 20,align = "left",
             style = list(color = "black",useHTML = TRUE)) %>%
    hc_xAxis(categories = x_data) %>%
    hc_yAxis(allowDecimals = FALSE,
             title = list(
               text = "Frequency"
             )) %>%
    hc_add_series(data = y_data,
                  name = "pathway"
                  ) %>%
    hc_plotOptions(
      series = list(stacking = FALSE,
                    events = list(
                      click = event_func
                    ))) %>%
    hc_tooltip(headerFormat = "<span>{point.key} has frequency {point.y}</span>
               <div style ='margin:0px;ma-xwidth:300px;overflow-y:hidden;'>",
               pointFormat = "<p class = 'tab3-tooltip-hc'>Metabolites: {point.detail}</p>",
               footerFormat = "</div>",
              # formatter = fomatterFunc,
               shared = TRUE,
               useHTML = TRUE) %>%
    hc_exporting(enabled = TRUE)
  return(hc)
}
