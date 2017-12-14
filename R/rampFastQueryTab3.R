#' Do fisher test for only one pathway from search result
#' clicked on highchart
#' @param pathwaydf a data frame resulting from rampFastPathFromMeta
#' @param total_analytes number of total genes or metabolites analyzed in the experiment (e.g. background) (default is 500, with assumption that analyte_type is "metabolite")
#' @param analyte_type "metabolites" or "genes" (default is "metabolites")
#' @param conpass password for database access (string)
#' @param dbname name of the mysql database (default is "ramp")
#' @param username username for database access (default is "root")
rampOneFisherTest <- function(pathwaydf,total_analytes=500,
	analyte_type="metabolites",conpass=NULL,
	dbname="ramp",username="root"){

  if(is.null(conpass)) {
        stop("Please define the password for the mysql connection")
  }

  contingencyTb <- matrix(0,nrow = 2,ncol = 2)
  colnames(contingencyTb) <- c("In Pathway","Not In Pathway")
  rownames(contingencyTb) <- c("All Metabolites","User's Metabolites")

  # Get the total number of analytes in the input pathway:
  pid <- unique(pathwaydf$pathwayRampId);
  if(length(pid)>1) {
	stop("This function is meant to do a Fisher's test on one pathway only (only input one info on one pathway")
  }
  # Retrieve the Ramp compound ids associated with the ramp pathway id and count them:
  query1 <- paste0("select rampId from analytehaspathway where pathwayRampId in (\"",
	pid,"\")")  

  con <- DBI::dbConnect(RMySQL::MySQL(), user = username,
        password = conpass,
        dbname = dbname)
  cids <- DBI::dbGetQuery(con,query1)[[1]]
  DBI::dbDisconnect(con)

  if(analyte_type=="metabolites") {
	tot_in_pathway <- length(grep("RAMP_C",cids))
  } else if (analyte_type=="genes") {
	tot_in_pathway <- length(grep("RAMP_G",cids))
  }
  else {stop("Please define analyte_type as 'metabolites' or 'genes'")}

  # Get the number of analytes
#  tot_in_pathway <- DBI::dbGetQuery(con,paste0("select count(*) from analyte 
#                                            where rampId in (select rampId from 
#                                            analytehaspathway where 
#                                            pathwayRampId in (select pathwayRampId 
#                                            from pathway where sourceId = \"",
#                                            unique(pathwaydf$pathwaysourceId),"\"));"))
#  DBI::dbDisconnect(con)
#  tot_in_pathway <- tot_in_pathway[[1]]

 # Now get the total number of metabolites
 con <- DBI::dbConnect(RMySQL::MySQL(), user = username,
	password = conpass, 
	dbname = dbname)
  
  tot_analytes <- DBI::dbGetQuery(con,"select count(*) from analyte;")[[1]]
  DBI::dbDisconnect(con)
  tot_out_pathway <- tot_analytes - tot_in_pathway

  # fill the rest of the table out 
  user_in_pathway <- nrow(pathwaydf)
  user_out_pathway <- total_analytes - user_in_pathway
  contingencyTb[1,1] <- tot_in_pathway
  contingencyTb[1,2] <- tot_out_pathway
  contingencyTb[2,1] <- user_in_pathway
  contingencyTb[2,2] <- user_out_pathway
  
  result <- stats::fisher.test(contingencyTb)
  pval <- round(result$p.value,4)
  return(pval)
}

#' Reformat the result of query (get pathways from analyte(s)) for input into barplot
#' function
#' 
#' each pathway contains all metabolites inside that pathway
#' @param df A dataframe that has information for bar plot
#' @return a list with the analyte names for each pathway that is represented in the list
rampGenerateBarPlot <- function(df){
  path_meta_list <- list()
  for (i in 1:nrow(df)){
    if (length(path_meta_list)==0){
      path_meta_list[[df[i,]$pathwaysourceId]] <- data.frame(metabolite = df[i,]$metabolite,stringsAsFactors = F)
    } else {
      path_meta_list[[df[i,]$pathwaysourceId]] <- 
		rbind(path_meta_list[[df[i,]$pathwaysourceId]],df[i,]$metabolite)
      path_meta_list[[df[i,]$pathwaysourceId]] <- 
	unique(path_meta_list[[df[i,]$pathwaysourceId]])
    }
  }
  return(path_meta_list)
}

#' Fast search given a list of metabolites
#' @param synonym a vector of synonym that need to be searched
#' @param find_synonym find all synonyms or just return same synonym (T/F)
#' @param conpass password for database access (string)
#' @param synonymOrIdS whether to return "synonyms" or "ids" (default is "ids")
#' @param dbname name of the mysql database (default is "ramp")
#' @param username username for database access (default is "root")
#' @return a list contains all metabolits as name and pathway inside.
#' 
##' Apply famil function...
#' @export
rampFastPathFromMeta<- function(synonym,
	find_synonym = FALSE,
	conpass=NULL,
	dbname="ramp",
	username="root",
	synonymOrIdS = "ids"){

  if(is.null(conpass)) {
        stop("Please define the password for the mysql connection")
  }

  now <- proc.time()

  if(synonymOrIdS == "synonyms"){
    synonym <- RaMP:::rampFindSynonymFromSynonym(synonym=synonym,
	find_synonym=find_synonym,
	conpass=conpass)
    list_metabolite <- unique(synonym)
    list_metabolite <- sapply(list_metabolite,shQuote)
    list_metabolite <- paste(list_metabolite,collapse = ",")
  } else if (synonymOrIdS == "ids"){
    source <- RaMP:::rampFindSourceRampId(sourceId=synonym, conpass=conpass)
    if (nrow(source)==0) {
	stop("Make sure you are actually inputting ids and not names (you have synonymOrIdS set to 'ids'")
	}
    list_metabolite <- unique(source)
    sourceIDTable <- list_metabolite
    list_metabolite <- list_metabolite$rampId
    list_metabolite <- sapply(list_metabolite,shQuote)
    list_metabolite <- paste(list_metabolite,collapse = ",")
  } else {
	stop("Make sure synonymOrIdS is set to 'synonyms' or 'ids'")
  }
  # Parse data to fit mysql
  # Can be simplified here
  if(list_metabolite=="") {
	stop("Unable to retrieve metabolites")
  }
 
 
  if(synonymOrIdS == "synonyms"){
    # Retrieve IDs for current name and all associated synonyms
    query1 <- paste0("select distinct Synonym,rampId from analytesynonym where Synonym in (",
                      list_metabolite,");")
    con <- connectToRaMP(dbname=dbname,username=username,conpass=conpass)
    df1<- DBI::dbGetQuery(con,query1)
    DBI::dbDisconnect(con)
    query2 <- paste0("select pathwayRampId,rampId from analytehaspathway where 
                      rampId in (select rampId from analytesynonym where Synonym in (",
                      list_metabolite,"));")
    con <- connectToRaMP(dbname=dbname,username=username,conpass=conpass)
    df2 <- DBI::dbGetQuery(con,query2)
    DBI::dbDisconnect(con)
  }else if (synonymOrIdS == "ids"){
    query2 <- paste0("select pathwayRampId,rampId from analytehaspathway where 
                      rampId in (",
                     list_metabolite,");")
    con <- connectToRaMP(dbname=dbname,username=username,conpass=conpass)
    df2 <- DBI::dbGetQuery(con,query2)
    DBI::dbDisconnect(con)
  }
  pathid_list <- df2$pathwayRampId
  pathid_list <- sapply(pathid_list,shQuote)
  pathid_list <- paste(pathid_list,collapse = ",")
  query3 <- paste0("select pathwayName,sourceId as pathwaysourceId,type as pathwaysource,pathwayRampId from pathway where pathwayRampId in (",
                    pathid_list,");")
  con <- connectToRaMP(dbname=dbname,username=username,conpass=conpass)
  df3 <- DBI::dbGetQuery(con,query3)
  DBI::dbDisconnect(con)
  if(synonymOrIdS == "synonyms"){
    mdf <- merge(df1,df2,all.x=T)
    mdf <- mdf[!is.na(mdf[,"pathwayRampId"]),]
    mdf <- merge(mdf,df3,all.x = T)
    colnames(mdf)[3] <- "metabolite"
    print("timing ...")
    print(proc.time()- now)
  } else if (synonymOrIdS == "ids"){
    mdf <- merge(df3,df2,all.x = T)
    mdf <- merge(mdf,source,all.x = T,by.y = "rampId")
  }
  return(mdf)
}

#' Generate data.frame from given files
#' 
#' identifing the file type, then it returns table output to 
#' shiny renderTable function as preview of searching data
#' 
#' @param infile a file object given from files 
#' @param synonymOrIdS whether to return "synonyms" or "ids" (default is "ids")
#' @param conpass password for database access (string)
#' @param dbname name of the mysql database (default is "ramp")
#' @param username username for database access (default is "root")
#' 
#' @return a data.frame either from multiple csv file
#' or search through by a txt file.
rampFileOfPathways <- function(infile,synonymOrIdS="ids",
	conpass=NULL,
	dbname="ramp",
	username="root"){
  name <- infile[[1,'name']]
    summary <- data.frame(pathway  = character(0),id = character(0),
                          source = character(0),metabolite = character(0))
    rampOut <- list()
    for (i in 1:length(infile[,1])){
      if(infile[[i,'type']]!="text/plain"){
        rampOut[[i]] <- utils::read.table(infile[[i,'datapath']])
        name <- infile[[i,'name']]
        print(infile[[i,'type']])
        rampOut[[i]]$new.col <- substr(name,1,nchar(name) - 4)
        colnames(rampOut[[i]]) <- c("pathway","id","source","metabolite")
        summary <- rbind(summary,rampOut[[i]])
      } else {
        rampOut <- readLines(infile[[i,'datapath']])
        summary <- RaMP:::rampFastPathFromMeta(rampOut,synonymOrIdS=synonymOrIdS,
		conpass=.conpass,username=username,dbname=dbname
		)
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

#' Calculate fisher test p-values from pathways returned when querying a list of genes
#' or metabolites (output of rampFastPathFromMeta)
#' @param rampOut The data frame generated by rampFastPathFromMeta
#' @param analyte_type specify whether to do test on 'metabolites' or 'genes'
#' @param total_analytes number of total genes or metabolites analyzed in the experiment (e.g. background) (default is 500, with assumption that analyte_type is "metabolite")
#' @param conpass password for database access (string)
#' @param dbname name of the mysql database (default is "ramp")
#' @param username username for database access (default is "root")
#' @return a data frame with three columns pathwayID, Number in pathway,
#' number out of pathway
rampFisherTestData <- function(rampOut,analyte_type="metabolites",totalanalytes=500,
	conpass=NULL,dbname="ramp",username="root"){
  if(analyte_type == "metabolites"){
    rampOut2 <- rampOut[grepl("RAMP_C_",rampOut$rampId),]
    total_analytes <- length(unique(rampOut$rampId[grepl("RAMP_C",rampOut$rampId)]))
  } else if (analyte_type == "genes") {
    rampOut2 <- rampOut[grepl("RAMP_G_",rampOut$rampId),]
    total_analytes <- length(unique(rampOut$rampId[grepl("RAMP_G",rampOut$rampId)]))
  } else {
    warning("Please input either 'metabolites' or 'genes' for analyte_type")
    return(NULL)
  }

  if(nrow(rampOut2)==0) {
	stop(paste0("You have selected ",analyte_type," yet no ",
		analyte_type," were found in the results table"))
	}
  fisher.pval <- c()
  for (i in unique(rampOut2$pathwaysourceId)) {
	tempOut <- rampOut2[which(rampOut2$pathwaysourceId==i),]
	fisher.pval <- c(fisher.pval,RaMP:::rampOneFisherTest(pathwaydf = tempOut,
		total_analytes=500,analyte_type=analyte_type,
		conpass=conpass,dbname=dbname,username=username))
	}
  fisher.adj.pval <- p.adjust(fisher.pval,method='fdr')
  # format output (retrieve pathway name for each unique source id first
  pathnames <- as.character(lapply(unique(rampOut2$pathwaysourceId), function(x)
	rampOut2$pathwayName[which(rampOut2$pathwaysourceId==x)[1]]))
  out=cbind(unique(rampOut2$pathwaysourceId),pathnames,fisher.pval,fisher.adj.pval)
  colnames(out)=c("SourceID","PathwayName","Fisher P-val","Fisher FDR-adj P-val")
  return(out)
}

