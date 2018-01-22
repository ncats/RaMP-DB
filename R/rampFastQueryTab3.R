#' Do fisher test for only one pathway from search result
#' clicked on highchart
#' @param pathwaydf a data frame resulting from rampFastPathFromMeta
#' @param total_metabolites number of metabolites analyzed in the experiment (e.g. background) (default is 1000; set to 'NULL' to retrieve total number of metabolites that map to any pathway in RaMP). Assumption that analyte_type is "metabolite")
#' @param total_genes number of genes analyzed in the experiment (e.g. background) (default is 20000, with assumption that analyte_type is "genes")
#' @param analyte_type "metabolites" or "genes" (default is "metabolites")
#' @param min_analyte if the number of analytes (gene or metabolite) in a pathway is
#' < min_analyte, do not report
#' @param conpass password for database access (string)
#' @param dbname name of the mysql database (default is "ramp")
#' @param username username for database access (default is "root")
#' @param host host name for database access (default is "localhost")
#' @return a dataframe with pathway enrichment (based on Fisher's test) results
#' @export
runFisherTest <- function(pathwaydf,total_metabolites=NULL,total_genes=20000,
                              analyte_type="metabolites",min_analyte=2,conpass=NULL,
                              dbname="ramp",username="root",
                          host = "localhost"){
    now <- proc.time()
  print("Fisher Testing ......")
  if(is.null(conpass)) {
    stop("Please define the password for the mysql connection")
  }

  if(analyte_type=="metabolites") {total_analytes=total_metabolites
	} else if (analyte_type=="genes") {
	total_analytes=total_genes
	} else {
	stop("Please define the analyte_type variable as 'metabolites' or 'genes'")
  }

  contingencyTb <- matrix(0,nrow = 2,ncol = 2)
  colnames(contingencyTb) <- c("In Pathway","Not In Pathway")
  rownames(contingencyTb) <- c("All Metabolites","User's Metabolites")

  # Get pathway ids that contain the user analytes
  pid <- unique(pathwaydf$pathwayRampId);
  list_pid <- sapply(pid,shQuote)
  list_pid <- paste(list_pid,collapse = ",")

  # Get the total number of metabolites that are mapped to pathways in RaMP (that's the default background)
   query <- "select * from analytehaspathway"
   con <- DBI::dbConnect(RMySQL::MySQL(), user = username,
         password = conpass,
         dbname = dbname,
         host = host)
   allids <- DBI::dbGetQuery(con,query)
   DBI::dbDisconnect(con)
   allids <- allids[!duplicated(allids),]

  if((analyte_type == "metabolites") && (is.null(total_metabolites))) {
	wiki_totanalytes <- length(unique(allids$rampId[grep("RAMP_C",allids[which(allids$type=="wiki"),"rampId"])]))
	react_totanalytes <- length(unique(allids$rampId[grep("RAMP_C",allids[which(allids$type=="reactome"),"rampId"])]))
	kegg_totanalytes <- length(unique(allids$rampId[grep("RAMP_C",allids[which(allids$type=="kegg"),"rampId"])]))
  }
  if(analyte_type=="genes") {
	wiki_totanalytes <- react_totanalytes <- kegg_totanalytes <- total_genes
  }

  print("Calculating p-values for pathways in input")
  # Retrieve the Ramp compound ids associated with the ramp pathway id and count them:
   query1 <- paste0("select rampId,pathwayRampId from analytehaspathway where pathwayRampId in (",
   list_pid,")")

   con <- DBI::dbConnect(RMySQL::MySQL(), user = username,
         password = conpass,
         dbname = dbname,
         host = host)
   cids <- DBI::dbGetQuery(con,query1)#[[1]]
   DBI::dbDisconnect(con)

   # Loop through each pathway, build the contingency table, and calculate Fisher's Exact
   # test p-value
   pval=totinpath=userinpath=pidused=c()
   for (i in pid) {
        curpathcids <- unique(cids[which(cids[,"pathwayRampId"]==i),"rampId"])
        if(analyte_type=="metabolites") {
                tot_in_pathway <- length(grep("RAMP_C",curpathcids))
        }else {
                tot_in_pathway <- length(grep("RAMP_G",curpathcids))
        }
	if(allids$type[which(allids$pathwayRampId==i)[1]] == "wiki") {
		total_analytes <- wiki_totanalytes
	} else if (allids$type[which(allids$pathwayRampId==i)[1]] == "reactome") {
                total_analytes <- react_totanalytes
        } else if (allids$type[which(allids$pathwayRampId==i)[1]] == "kegg") {
                total_analytes <- kegg_totanalytes
        } else {stop("Couldn't find pathway type for current pathway!")}

 	tot_out_pathway <- total_analytes - tot_in_pathway
	  # fill the rest of the table out
	  user_in_pathway <- nrow(pathwaydf[which(pathwaydf$pathwayRampId==i),])
	  user_out_pathway <- length(unique(pathwaydf$rampId)) - user_in_pathway
	  contingencyTb[1,1] <- tot_in_pathway
	  contingencyTb[1,2] <- tot_out_pathway
	  contingencyTb[2,1] <- user_in_pathway
	  contingencyTb[2,2] <- user_out_pathway
	  result <- stats::fisher.test(contingencyTb)
	  pval <- c(pval,result$p.value )
	  userinpath<-c(userinpath,user_in_pathway)
	  totinpath<-c(totinpath,tot_in_pathway)
	  pidused <- c(pidused,i)
  } # end for loop

  # Now run fisher's tests for all other pids
   query <- "select distinct(pathwayRampId) from analytehaspathway where type != 'hmdb';"
   con <- DBI::dbConnect(RMySQL::MySQL(), user = username,
         password = conpass,
         dbname = dbname,
         host = host)
   allpids <- DBI::dbGetQuery(con,query)
   DBI::dbDisconnect(con)
   pidstorun <- setdiff(allpids[,1],pid)
   pidstorunlist <- sapply(pidstorun,shQuote)
   pidstorunlist <- paste(pidstorunlist,collapse = ",")

   query2 <- paste0("select rampId,pathwayRampId from analytehaspathway where pathwayRampId in (",
   pidstorunlist,")")

   con <- DBI::dbConnect(RMySQL::MySQL(), user = username,
         password = conpass,
         dbname = dbname,
         host = host)
   restcids <- DBI::dbGetQuery(con,query2)#[[1]]
   DBI::dbDisconnect(con)

   query1 <- paste0("select rampId,pathwayRampId from analytehaspathway;")

   con <- DBI::dbConnect(RMySQL::MySQL(), user = username,
         password = conpass,
         dbname = dbname,
         host = host)
   allcids <- DBI::dbGetQuery(con,query1)#[[1]]
   DBI::dbDisconnect(con)

   print("Calculating p-values for all other pathways")
   print(paste0(length(pidstorun),"pathways"))
   # calculating p-values for all other pathways
   count=1;
   for (i in pidstorun) {
	if(( count %% 100) ==0) {print(paste0("Processed ",count))}
	count=count+1
	user_in_pathway=0
        curpathcids <- unique(allcids[which(allcids[,"pathwayRampId"]==i),"rampId"])
        if(analyte_type=="metabolites") {
                tot_in_pathway <- length(grep("RAMP_C",curpathcids))
        }else {
                tot_in_pathway <- length(grep("RAMP_G",curpathcids))
        }
        if(allids$type[which(allids$pathwayRampId==i)[1]] == "wiki") {
                total_analytes <- wiki_totanalytes
        } else if (allids$type[which(allids$pathwayRampId==i)[1]] == "reactome") {
                total_analytes <- react_totanalytes
        } else if (allids$type[which(allids$pathwayRampId==i)[1]] == "kegg") {
                total_analytes <- kegg_totanalytes
        } else if (allids$type[which(allids$pathwayRampId==i)[1]] == "hmdb") {
		total_analytes <- NULL
	} else {stop("Couldn't find pathway type for current pathway!")}

	if(is.null(total_analytes)) {next;}
        tot_out_pathway <- total_analytes - tot_in_pathway
          # fill the rest of the table out
          user_out_pathway <- length(unique(pathwaydf$rampId))
          #user_out_pathway <- total_analytes - user_in_pathway
          contingencyTb[1,1] <- tot_in_pathway
          contingencyTb[1,2] <- tot_out_pathway
          contingencyTb[2,1] <- user_in_pathway
          contingencyTb[2,2] <- user_out_pathway
          result <- stats::fisher.test(contingencyTb)
          pval <- c(pval,result$p.value )
          userinpath<-c(userinpath,user_in_pathway)
          totinpath<-c(totinpath,tot_in_pathway)
         # pidused <- c(pidused,i)
  } # end for loop

  fdr <- stats::p.adjust(pval,method="fdr")
  holm <- stats::p.adjust(pval,method="holm")
  print(paste0("Calculated p-values for ",length(pval)," pathways"))

  # format output (retrieve pathway name for each unique source id first
  out <- data.frame(pathwayRampId=c(pidused,pidstorun), Pval=pval,FDR.Adjusted.Pval=fdr,
	Holm.Adjusted.Pval=holm,
	Num_In_Path=userinpath,Total_In_Path=totinpath)
  out2 <- merge(out,pathwaydf[,c("pathwayName","pathwayRampId","pathwaysourceId",
	"pathwaysource","pathwayRampId")],
	by="pathwayRampId",all.x=TRUE)
  finout <- out2[,c("pathwayName", "Pval", "FDR.Adjusted.Pval",
	"Holm.Adjusted.Pval","pathwaysourceId",
	"pathwaysource","Num_In_Path","Total_In_Path","pathwayRampId")]
  finout=finout[!duplicated(finout),]

  return(finout[which(finout$Num_In_Path>=min_analyte),])
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
      path_meta_list[[df[i,]$pathwaysourceId]] <- data.frame(metabolite = df[i,]$rampId,stringsAsFactors = F)
    } else {
      path_meta_list[[df[i,]$pathwaysourceId]] <-
		rbind(path_meta_list[[df[i,]$pathwaysourceId]],df[i,]$rampId)
      path_meta_list[[df[i,]$pathwaysourceId]] <-
	unique(path_meta_list[[df[i,]$pathwaysourceId]])
    }
  }
  return(path_meta_list)
}

#' Function that search analytes (gene or compounds)  or a list of analytes and
#' returns associated pathways
#'
#' @param analytes a vector of analytes (genes or metabolites) that need to be searched
#' @param find_synonym find all synonyms or just return same synonym (T/F)
#' @param conpass password for database access (string)
#' @param NameOrIds whether input is "names" or "ids" (default is "ids")
#' @param host host name for database access (default is "localhost")
#' @param dbname name of the mysql database (default is "ramp")
#' @param username username for database access (default is "root")
#' @return a list contains all metabolits as name and pathway inside.
#' @export
rampFastPathFromMeta<- function(analytes,
	find_synonym = FALSE,
	conpass=NULL,
	host = "localhost",
	dbname="ramp",
	username="root",
	NameOrIds = "ids"){

  if(is.null(conpass)) {
        stop("Please define the password for the mysql connection")
  }

  now <- proc.time()

  if(NameOrIds == "names"){
    synonym <- RaMP:::rampFindSynonymFromSynonym(synonym=analytes,
	find_synonym=find_synonym,
	conpass=conpass)
    colnames(synonym)[1]="commonName"
    synonym$commonName <- tolower(synonym$commonName)
    if(nrow(synonym)==0) {
	stop("Could not find any matches to the analytes entered.  If pasting, please make sure the names are delimited by end of line (not analyte per line)\nand that you are selecting 'names', not 'ids'");
    }
    # Get all unique RaMP ids and call it list_metabolite
    list_metabolite <- unique(synonym$rampId)
    list_metabolite <- sapply(list_metabolite,shQuote)
    list_metabolite <- paste(list_metabolite,collapse = ",")
  } else if (NameOrIds == "ids"){
    sourceramp <- RaMP:::rampFindSourceRampId(sourceId=analytes, conpass=conpass)
    if (nrow(sourceramp)==0) {
	stop("Make sure you are actually inputting ids and not names (you have NameOrIds set to 'ids'. If you are, then no ids were matched in the RaMP database.")
	}
    # get all unique RaMP ids and call it list_metabolite
    list_metabolite <- unique(sourceramp$rampId)
    #sourceIDTable <- list_metabolite
    #list_metabolite <- list_metabolite$rampId
    list_metabolite <- sapply(list_metabolite,shQuote)
    list_metabolite <- paste(list_metabolite,collapse = ",")
  } else {
	stop("Make sure NameOrIds is set to 'names' or 'ids'")
  }
  # Parse data to fit mysql
  # Can be simplified here
  if(list_metabolite=="") {
	stop("Unable to retrieve metabolites")
  }

  # Now using the RaMP compound id, retrieve associated pathway ids
    query2 <- paste0("select pathwayRampId,rampId from analytehaspathway where
                      rampId in (",
                     list_metabolite,");")
    con <- connectToRaMP(dbname=dbname,username=username,conpass=conpass,host = host)
    df2 <- DBI::dbGetQuery(con,query2)
    DBI::dbDisconnect(con)
  pathid_list <- df2$pathwayRampId
  pathid_list <- sapply(pathid_list,shQuote)
  pathid_list <- paste(pathid_list,collapse = ",")
  # With pathway ids, retrieve pathway information
  query3 <- paste0("select pathwayName,sourceId as pathwaysourceId,type as pathwaysource,pathwayRampId from pathway where pathwayRampId in (",
                    pathid_list,");")
  con <- connectToRaMP(dbname=dbname,username=username,conpass=conpass)
  df3 <- DBI::dbGetQuery(con,query3)
  DBI::dbDisconnect(con)

  #Format output
  mdf <- merge(df3,df2,all.x = T)

  # And with rampIds (list_metabolite), get common names when Ids are input
  if(NameOrIds == "ids"){
     list_analytes <- sapply(analytes,shQuote)
     list_analytes <- paste(list_analytes,collapse = ",")
  query4 <-paste0("select sourceId,commonName,rampId from source where sourceId in (",list_analytes,");")

  con <- connectToRaMP(dbname=dbname,username=username,conpass=conpass,host = host)
  df4 <- DBI::dbGetQuery(con,query4)
  DBI::dbDisconnect(con)
  mdf <- merge(mdf,df4,,all.x = T,by.y = "rampId")
  mdf$commonName=tolower(mdf$commonName)
 } else{ # Just take on the name
  mdf <- merge(mdf,synonym,all.x = T,by.y = "rampId")
 }
  out<-mdf[!duplicated(mdf),]

  return(out[which(out$pathwaysource!="hmdb"),])
}

#' Generate data.frame from given files
#'
#' identifing the file type, then it returns table output to
#' shiny renderTable function as preview of searching data
#'
#' @param infile a file object given from files
#' @param NameOrIds whether to return "synonyms" or "ids" (default is "ids")
#' @param conpass password for database access (string)
#' @param dbname name of the mysql database (default is "ramp")
#' @param username username for database access (default is "root")
#' @param host host name for database access (default is "localhost")
#' @return a data.frame either from multiple csv file
#' or search through by a txt file.
rampFileOfPathways <- function(infile,NameOrIds="ids",
	conpass=NULL,
	dbname="ramp",
	username="root",
	host = "localhost"){
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
        summary <- RaMP::rampFastPathFromMeta(rampOut,NameOrIds=NameOrIds,
		conpass=conpass,username=username,dbname=dbname,
		host = host
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

#' Perform fuzzy multiple linkage partitioning clustering on pathways identified by
#' Fisher's test
#'
#' @param fishers_df The data frame generated by runFisherTest
#' @param analyte_type specify whether to use metabolite or gene pathway overlap matrix
#' @param perc_analyte_overlap Minimum overlap for pathways to be considered similar
#' (Default = 0.5)
#' @param min_pathway_tocluster Minimum number of 'similar' pathways required to start
#' a cluster (medoid) (Default = 3)
#' @param perc_pathway_overlap Minimum overlap for clusters to merge (Default = 0.5)
#' @param p_cutoff Filter pathways to cluster by FDR p value threshold (Default = 0.05)
#'
#' @return a list of clusters identified by the algorithm. Each entry of the list is
#' a cluster, containing a vector of pathways in the cluster
#' @export
find_clusters <- function(fishers_df,analyte_type,perc_analyte_overlap = 0.5,
	min_pathway_tocluster = 3,perc_pathway_overlap = 0.5,p_cutoff=0.05){
  if(perc_analyte_overlap <= 0 || perc_analyte_overlap >= 1 || perc_pathway_overlap <= 0 || perc_pathway_overlap >= 1){
    return(NULL)
  }
  similarity_matrix_list<-load_overlap_matrices()
  if(analyte_type=="metabolites"){
    similarity_matrix = similarity_matrix_list[[2]]
  } else if(analyte_type=="genes"){
    similarity_matrix = similarity_matrix_list[[1]]
  }
  pathway_list<-fishers_df[,9]
  #pathway_list<-pathway_list[which(fishers_df[,4] < p_cutoff)]

  pathway_indices<-match(pathway_list,rownames(similarity_matrix))

  if(length(which(is.na(pathway_indices)))>0){
    pathway_indices<-pathway_indices[-which(is.na(pathway_indices))]
  }

  pathway_matrix<-similarity_matrix[pathway_indices,pathway_indices]
  unmerged_clusters<-apply(pathway_matrix, 1, function(x){
    if(length(which(x>=perc_analyte_overlap))>(min_pathway_tocluster+1)){
      return(colnames(pathway_matrix)[which(x>=perc_analyte_overlap)])
    } else {
      return(NA)
    }
  })
  # Remove the unmerged clusters
  if(length(which(is.na(unmerged_clusters)))>0){
    unmerged_clusters<-unmerged_clusters[-which(is.na(unmerged_clusters))]
  }

  if(length(unmerged_clusters)==0){
    #stop("No medoids found, make perc_analyte_overlap or min_pathway_tocluster smaller")
    return(rep("Did not cluster",times = nrow(fishers_df)))
  }

  # Evaluate similarity between clusters
  cluster_similarity<-matrix(0,ncol = length(unmerged_clusters),nrow = length(unmerged_clusters))
  for(i in 1:length(unmerged_clusters)){
    for(j in 1:length(unmerged_clusters)){
      cluster_similarity[i,j]<-length(intersect(unmerged_clusters[[i]],unmerged_clusters[[j]]))/
        length(unique(c(unmerged_clusters[[i]],unmerged_clusters[[j]])))
    }
  }
  colnames(cluster_similarity)<-rownames(cluster_similarity)<-names(unmerged_clusters)
  unmerged_cluster_similarity<-cluster_similarity

  cluster_list<-unmerged_clusters

  # Merge Clusters
  count = 1
  while(length(which(cluster_similarity >= perc_pathway_overlap)) > nrow(cluster_similarity)){
    cluster_similarity_mod<-cluster_similarity
    for(i in 1:nrow(cluster_similarity_mod)){
      cluster_similarity_mod[i,i]<-0
    }

    clusters_to_merge<-which(cluster_similarity_mod == max(cluster_similarity_mod), arr.ind = TRUE)
    clusters_to_merge<-unique(t(apply(clusters_to_merge, 1, sort)))

    for(i in 1:nrow(clusters_to_merge)){
      if(!is.na(cluster_list[[clusters_to_merge[i,1]]])&&!is.na(cluster_list[[clusters_to_merge[i,2]]])){
        cluster_list[[clusters_to_merge[i,1]]]<-unique(unlist(cluster_list[c(clusters_to_merge[i,1],clusters_to_merge[i,2])]))
        cluster_list[[clusters_to_merge[i,2]]]<-NA
      }
    }

    if(length(which(is.na(cluster_list)))>0){
      cluster_list<-cluster_list[-which(is.na(cluster_list))]
    }

    cluster_similarity<-matrix(0,ncol = length(cluster_list),nrow = length(cluster_list))
    for(i in 1:length(cluster_list)){
      for(j in 1:length(cluster_list)){
        cluster_similarity[i,j]<-length(intersect(cluster_list[[i]],cluster_list[[j]]))/
          length(unique(c(cluster_list[[i]],cluster_list[[j]])))
      }
    }

    if(nrow(cluster_similarity)==1){
      #stop("Clusters converged, use larger perc_pathway_overlap")
      return(rep(1,times = nrow(fishers_df)))
    }
    count = count + 1
    if(count == length(unmerged_clusters)+1){
      #stop("ERROR: while loop failed to terminate")
      return(rep(1,times = nrow(fishers_df)))
    }
  }
  colnames(cluster_similarity) = rownames(cluster_similarity) = paste0("cluster_",c(1:length(cluster_list)))
  return(cluster_list)
}

#' Filter pathways by p-value cutoff for display and clustering
#' @param fishers_df The data frame generated by runFisherTest
#' @param p_holmadj_cutoff return pathways where Holm adjusted pvalues are < p_holmadj_cutoff
#' @param p_fdradj_cutoff return pathways where FDR adjusted pvalues are < p_fdradj_cutoff
#' @return dataframe of fishers results with only significant pathways
#' @export
FilterFishersResults<-function(fishers_df,p_holmadj_cutoff=NULL,p_fdradj_cutoff=NULL){
	if(!is.null(p_holmadj_cutoff)) {
  		return(fishers_df[which(fishers_df[,"Holm.Adjusted.Pval"] < p_holmadj_cutoff),])
	} else if (!is.null(p_fdradj_cutoff)) {
		return(fishers_df[which(fishers_df[,"FDR.Adjusted.Pval"] < p_fdradj_cutoff),])
	} else {
		stop("Please set a cutoff for Holm Adjusted pvalues (p_holmadj_cutoff paramter) or FDR Adjusted pvalues
			(p_fdradj_cutoff)")
	}
}
