#' Do fisher test for only one pathway from search result
#' clicked on highchart
#' @param pathway_list a vector of all pathways searched from given metabolites
#' @param pathway single pathway name that is applied by Fisher test
#' @param num_of_meta integer that describe how many initial metabolites
#' @param conpass password for database access (string)
#' @param dbname name of the mysql database (default is "ramp")
#' @param username username for database access (default is "root")
rampOneFisherTest <- function(pathway_list,pathway,num_of_meta,conpass=NULL,
	dbname="ramp",username="root"){

  if(is.null(conpass)) {
        stop("Please define the password for the mysql connection")
  }

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

  con <- DBI::dbConnect(RMySQL::MySQL(), user = username,
	password = conpass, 
	dbname = dbname)
  
  tot_metabolites <- DBI::dbGetQuery(con,"select count(*) from analyte;")
  DBI::dbDisconnect(con)
  tot_metabolites <- tot_metabolites[[1]]
  tot_out_pathway <- tot_metabolites - tot_in_pathway
  user_in_pathway <- nrow(pathway_list[[pathway]])
  user_out_pathway <- num_of_meta - user_in_pathway
  contingencyTb[1,1] <- tot_in_pathway[[1]]
  contingencyTb[1,2] <- tot_out_pathway[[1]]
  contingencyTb[2,1] <- user_in_pathway
  contingencyTb[2,2] <- user_out_pathway
  
  result <- stats::fisher.test(contingencyTb)
  result$p.value <- round(result$p.value,4)
  return(result)

}
#' Return a list named by pathway
#' 
#' each pathway contains all metabolites inside that pathway
#' @param df A dataframe that has information for bar plot
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
#' @param pathway_meta_list The list contains pathway as names and metabolites in that pathway 
#' under each names (list)
#' @param num_user_metabolites number of metabolites given by user when they want to 
#' search for pathways.
#' @param FisherPathwayTable Fisher Pathway Table
#' @param con a connection object returned from the function connectToRaMP()
#' @return a data.frame contains all fisher test result with pathway name
#' as column name
rampFisherTest <- function(pathway_meta_list,
			   num_user_metabolites,
			   FisherPathwayTable,con){
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
  
  cumulative_fisher_results <- lapply(cumulative_fisher_results,stats::fisher.test)
  
  return(cumulative_fisher_results)
}
#' Fast search given a list of metabolites
#' @param synonym a vector of synonym that need to be searched
#' @param find_synonym bool if find all synonyms or just return same synonym
#' @param conpass password for database access (string)
#' @param dbname name of the mysql database (default is "ramp")
#' @param username username for database access (default is "root")
#' @return a list contains all metabolits as name and pathway inside.
#' 
#' Apply famil function...
<<<<<<< HEAD
rampFastPathFromMeta<- function(input,find_synonym = FALSE,synonymOrIdS = "ids"){
=======
rampFastPathFromMeta<- function(synonym,find_synonym = FALSE,
	conpass=NULL,dbname="ramp",username="root"){

  if(is.null(conpass)) {
        stop("Please define the password for the mysql connection")
  }


>>>>>>> 3a5564b02de570058f5719dca0fede3fb525f7be
  # progress<- shiny::Progress$new()
  # progress$set(message = "Querying databases ...",value = 0)
  now <- proc.time()
  # con <- dbConnect(MySQL(), user = "root", password = "Ramp340!", dbname = "mathelabramp")
  # on.exit(dbDisconnect(con))
  # find synonym
<<<<<<< HEAD
  if(synonymOrIdS == "synonyms"){
    synonym <- rampFindSynonymFromSynonym(input,find_synonym=find_synonym)
    list_metabolite <- unique(synonym)
    
    list_metabolite <- sapply(list_metabolite,shQuote)
    list_metabolite <- paste(list_metabolite,collapse = ",")
  } else if (synonymOrIdS == "ids"){
    source <- rampFindSourceRampId(input)
    list_metabolite <- unique(source)
    sourceIDTable <- list_metabolite
    list_metabolite <- list_metabolite$rampId
    list_metabolite <- sapply(list_metabolite,shQuote)
    list_metabolite <- paste(list_metabolite,collapse = ",")
    
  }
  # Parse data to fit mysql
  # Can be simplified here
  
  if(synonymOrIdS == "synonyms"){
    query1 <- paste0("select distinct Synonym,rampId from analytesynonym where Synonym in (",
                      list_metabolite,");")
    df1<- DBI::dbGetQuery(con,query1)
    query2 <- paste0("select pathwayRampId,rampId from analytehaspathway where 
                      rampId in (select rampId from analytesynonym where Synonym in (",
                      list_metabolite,"));")
    df2 <- DBI::dbGetQuery(con,query2)
  }else if (synonymOrIdS == "ids"){
    query2 <- paste0("select pathwayRampId,rampId from analytehaspathway where 
                      rampId in (",
                     list_metabolite,");")
    df2 <- DBI::dbGetQuery(con,query2)
    
  }
=======

  synonym <- rampFindSynonymFromSynonym(synonym,find_synonym=find_synonym,
	conpass=conpass)
    
  list_metabolite <- unique(synonym)
  list_metabolite <- sapply(list_metabolite,shQuote)
  list_metabolite <- paste(list_metabolite,collapse = ",")
  query1 <- paste0("select distinct Synonym,rampId from analytesynonym where Synonym in (",
                    list_metabolite,");")


  con <- DBI::dbConnect(RMySQL::MySQL(), user = username,
        password = conpass,
        dbname = dbname)

  tot_metabolites <- DBI::dbGetQuery(con,"select count(*) from analyte;")
  df1<- DBI::dbGetQuery(con,query1)
  DBI::dbDisconnect(con)

  query2 <- paste0("select pathwayRampId,rampId from analytehaspathway where 
                    rampId in (select rampId from analytesynonym where Synonym in (",
                    list_metabolite,"));")
  con <- DBI::dbConnect(RMySQL::MySQL(), user = username,
        password = conpass,
        dbname = dbname)
  df2 <- DBI::dbGetQuery(con,query2)
  DBI::dbDisconnect(con)

>>>>>>> 3a5564b02de570058f5719dca0fede3fb525f7be
  id_list <- df2[,1]
  id_list <- sapply(id_list,shQuote)
  id_list <- paste(id_list,collapse = ",")
  query3 <- paste0("select pathwayName,sourceId as pathwaysourceId,type as pathwaysource,pathwayRampId from pathway where pathwayRampId in (",
                    id_list,");")
  con <- DBI::dbConnect(RMySQL::MySQL(), user = username,
        password = conpass,
        dbname = dbname)
  df3 <- DBI::dbGetQuery(con,query3)
<<<<<<< HEAD
  if(synonymOrIdS == "synonyms"){
    mdf <- merge(df1,df2,all.x=T)
    mdf <- mdf[!is.na(mdf[,3]),]
    mdf <- merge(mdf,df3,all.x = T)
    # mdf <- mdf[!is.na(mdf[,3]),]
    colnames(mdf)[3] <- "metabolite"
    print("timing ...")
    print(proc.time()- now)
  } else if (synonymOrIdS == "ids"){
    mdf <- merge(df3,df2,all.x = T)
    mdf <- merge(mdf,source,all.x = T,by.y = "rampId")
  }
=======
  DBI::dbDisconnect(con)

  mdf <- merge(df1,df2,all.x=T)
  mdf <- mdf[!is.na(mdf[,3]),]
  mdf <- merge(mdf,df3,all.x = T)
  # mdf <- mdf[!is.na(mdf[,3]),]
  colnames(mdf)[3] <- "metabolite"
  print("timing ...")
  print(proc.time()- now)
  return(mdf)
  mdf <- mdf[,c(4,5,6,3)]
  mdf <- unique(mdf)
  
>>>>>>> 3a5564b02de570058f5719dca0fede3fb525f7be
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
rampFileOfPathways <- function(infile){
  name <- infile[[1,'name']]
    summary <- data.frame(pathway  = character(0),id = character(0),
                          source = character(0),metabolite = character(0))
    rampOut <- list()
    for (i in 1:length(infile[,1])){
      if(infile[[i,'type']]!="text/plain"){
        rampOut[[i]] <- utils::read.csv(infile[[i,'datapath']])
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
rampHcOutput <- function(x_data,y_data,type = 'column',event_func){
  fomatterFunc <- highcharter::JS("function(){
                        html = '<strong> Pathway ' + this.x +' has frequency: ' + this.y +'. metabolites are :'
                        '+this.y.detail+'</strong>;'
                        return html;
                     }")
  hc<-highcharter::highchart() %>%
    highcharter::hc_chart(type = type,
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
    highcharter::hc_title(text = "<strong>Search result from given metabolites</strong>",
             margin = 20,align = "left",
             style = list(color = "black",useHTML = TRUE)) %>%
    highcharter::hc_xAxis(categories = x_data) %>%
    highcharter::hc_yAxis(allowDecimals = FALSE,
             title = list(
               text = "Frequency"
             )) %>%
    highcharter::hc_add_series(data = y_data,
                  name = "pathway"
                  ) %>%
    highcharter::hc_plotOptions(
      series = list(stacking = FALSE,
                    events = list(
                      click = event_func
                    ))) %>%
    highcharter::hc_tooltip(headerFormat = "<span>{point.key} has frequency {point.y}</span>
               <div style ='margin:0px;ma-xwidth:300px;overflow-y:hidden;'>",
               pointFormat = "<p class = 'tab3-tooltip-hc'>Metabolites: {point.detail}</p>",
               footerFormat = "</div>",
              # formatter = fomatterFunc,
               shared = TRUE,
               useHTML = TRUE) %>%
    highcharter::hc_exporting(enabled = TRUE)
  return(hc)
}
#' Generate raw data for fisher test based on the given output
#' Out put is from rampFastMetaFromPathway
#' Required format 
<<<<<<< HEAD
#' @param rampOut The data frame generated by rampFastPathFromMeta
#' @metaboliteOrGene A bool value that specifiy if the fisher test data for
#' metabolites or genes is the output
#' @return a data frame with three columns pathwayID, Number in pathway,
#' number out of pathway
rampFisherTestData <- function(rampOut,metaboliteOrGene){
  if(metaboliteOrGene == "metabolites"){
    rampOut2 <- rampOut[grepl("RAMP_C_",rampOut$rampId),]
  } else if (metaboliteOrGene == "genes") {
    rampOut2 <- rampOut[grepl("RAMP_G_",rampOut$rampId),]
  } else {
    warning("Please input either 'metabolites' or 'genes' for metaboliteOrGene")
    return(NULL)
  }
  df <- aggregate(rampOut2$rampId,
                  list(rampOut2$pathwayRampId),FUN = function(x){
                    x <- unique(x)
                    length(x)
                  }
                  )
  total_metabolites <- length(unique(rampOut$rampId[grepl("RAMP_C",rampOut$rampId)]))
  total_genes <- length(unique(rampOut$rampId[grepl("RAMP_G",rampOut$rampId)]))
  if(metaboliteOrGene == "metabolites"){
    df$newCol <- total_metabolites - df$x
  } else if(metaboliteOrGene == "genes"){
    df$newCol <- total_genes - df$x
  }
  print("Hello")
  colnames(df) <- c("pathwayRampId","inPathway","outPathway")
  df
}
=======
>>>>>>> 3a5564b02de570058f5719dca0fede3fb525f7be
