# First tab panel
dataInput_name <- eventReactive(input$submit_compName,{
  progress <- shiny::Progress$new()

  on.exit(progress$close())

  progress$set(message = "Querying databases to find pathways ...", value = 0)
  progress$inc(0.3,detail = paste("Send Query ..."))

  rampOut <- RaMP::rampFastPathFromMeta(analytes=input$KW_synonym,
                                        NameOrIds=input$NameOrId,
                                        conpass=.conpass,
                                        host = .host)
  progress$inc(0.7,detail = paste("Done!"))
  return (rampOut)
})

summary_path_out<- eventReactive(input$submit_compName,{
  if (!is.null(nrow(dataInput_name()))){
    return (paste0("There are ",nrow(dataInput_name())," pathways returned for ",
                   input$KW_synonym))
  } else{
    return ("Given metabolites have no search result.")
  }
})

output$summary_path <- renderText({
  summary_path_out()
})

observe({
  if(input$NameOrId == "names"){
    choices <- kw_analyte[grepl(input$compName,kw_analyte,ignore.case=TRUE)]
    choices <- choices[order(nchar(choices),choices)]
    if(is.null(choices))
      return(NULL)
    if(length(choices) >10 ){
      choices <- choices[1:10]
    }
    isolate({
      updateSelectInput(session, "KW_synonym",
                        label = "Select from the list",
                        choices = choices, selected = head(choices,1)
      )
    })
  } else if (input$NameOrId == "ids"){
    choices <- kw_source[grepl(input$compName,kw_source,fixed = T)]
    choices <- choices[order(nchar(choices),choices)]
    if(is.null(choices))
      return(NULL)
    if(length(choices) >10 ){
      choices <- choices[1:10]
    }
    isolate({
      updateSelectInput(session, "KW_synonym",
                        label = "Select from the list",
                        choices = choices, selected = head(choices,1)
      )
    })
  }
})



output$result3 <- DT::renderDataTable({
  out_stc <- dataInput_name()
  out_stc[,c("pathwayName","pathwaysourceId","pathwaysource")]
})

output$comp_report <- downloadHandler(filename = function() {
  paste0(input$KW_synonym, ".csv")
},
content = function(file) {
  rampOut <- dataInput_name()
  rampOut <- data.frame(rampOut)
  write.csv(rampOut, file, row.names = FALSE, sep = ",")
})


# Second Tab
#
#
# rea_detector <- reactiveValues(num = NULL)
# 
# observe({
#   input$sub_mul_tab3
#   #isGeneMetabolites$content <- 'metabolites'
# })

data_mul_name <- eventReactive(input$sub_mul_tab3,{
  print(input$input_mul_tab3)
  parsedinput <- paste(strsplit(input$input_mul_tab3,"\n")[[1]])
  print(parsedinput)
  if(length(parsedinput)==0) {metabsearch=NULL} else{
  	metabsearch <- RaMP::rampFastPathFromMeta(analytes=parsedinput,
                             NameOrIds=input$NameOrSourcemult,
                             conpass=.conpass,
                             host = .host)
	  print(input$input_mul_tab3_genes)
  }

  parsedinputg <- paste(strsplit(input$input_mul_tab3_genes,"\n")[[1]])
  print(parsedinputg)
  if(length(parsedinputg)==0) {genesearch=NULL} else{
	  genesearch <- RaMP::rampFastPathFromMeta(analytes=parsedinputg,
                             NameOrIds=input$NameOrSourcemult_genes,
                             conpass=.conpass,
                             host = .host)
  }
  print(paste("metabsearch: ",ncol(metabsearch)))
  print(paste("genesearch: ",ncol(genesearch)))
  print(paste0("DIM of data_mul_name",nrow(rbind(metabsearch,genesearch))))
  rbind(metabsearch,genesearch)
})

# data_mul_file <- eventReactive(input$sub_file_tab3,{
#   infile <- input$inp_file_tab3
#   if (is.null(infile))
#     return(NULL)
# 
#   RaMP:::rampFileOfPathways(infile,conpass=.conpass,host = .host,NameOrIds=input$NameOrSourcemult)
# })
# 
# observe({
#   input$sub_file_tab3
#   rea_detector$num <- 2
# })

# Download table in a csv file.
output$tab3_mul_report <- downloadHandler(filename = function(){
  # if (rea_detector$num == 1){
      paste0("pathwayFromMetabolitesOutput.csv")
  # } else if (rea_detector$num == 2){
  #   infile <- input$inp_file_tab3
  #   paste0(infile[[1,'name']],"Output",".csv")
  # }

},
content = function(file) {
  # if (rea_detector$num == 1){
      rampOut <- data_mul_name()[,c("pathwayName","pathwaysourceId",
                                    "pathwaysource","commonName")]
  # } else if (rea_detector$num == 2){
  # }
  write.csv(rampOut,file,row.names = FALSE)
}
)

output$summary_mulpath_out<- DT::renderDataTable({
  if(is.null(data_mul_name())) {
    out <- data.frame(Query=NA,Freq=NA)
  } else {
      temp <- data_mul_name()
    }
    out <- as.data.frame(table(temp$commonName))
    colnames(out)[1] <- "Query"
  out
},rownames=FALSE)

output$preview_multi_names <- DT::renderDataTable({
  #if(is.null(isGeneMetabolites$content))
  #  return("Waiting for input")
  if(is.null(data_mul_name())) {
    return("No input found")
  } else {

  # if(rea_detector$num == 1){
    tb <- data_mul_name()[,c("pathwayName","pathwaysourceId",
                             "pathwaysource","commonName")]
#   } else if (rea_detector$num == 2) {
# 	    tb <- data_mul_file()[,c("pathwayName","pathwaysourceId",
#                              "pathwaysource","commonName")]
#   }
  return(tb)
 }
}
,rownames = FALSE)


fisherTestResult <- eventReactive(input$runFisher,{

    out <- RaMP::runCombinedFisherTest(req(data_mul_name()),
                               conpass=.conpass)
    print("Results generated")
    print(paste0("Fisher results size:",nrow(out[[1]])))
  out
})

output$summary_fisher <- DT::renderDataTable({
  if(!is.null(fisherTestResult())) {
    data <- fisherTestResult()
    out=as.data.frame(table(data$pathwaysource))
    colnames(out)[1]="Pathway_Source"
  } else {
    out <- data.frame(Pathway_Source=NA, Freq=NA)
  }
  out
},rownames=FALSE,filter="top")


output$num_mapped_namesids <- renderText({
	  data <- data_mul_name()
	  inputlist <- input$input_mul_tab3
	if(!is.null(data)) {
		parsedinput <- paste(strsplit(inputlist,"\n")[[1]])
		print(paste0("Found ",length(unique(data$commonName))," out of ",
		length(parsedinput)))
	}
})

fisherTestResultSignificant<-eventReactive(input$runFisher,{
  if(!is.null(fisherTestResult())){
    result<-RaMP::FilterFishersResults(fisherTestResult(),p_holmadj_cutoff=as.numeric(input$p_holmadj_cutoff))
    return(result)
  }else{
    return(NULL)
  }
})

cluster_output<-eventReactive(c(input$runFisher,input$runClustering),{
  data <- fisherTestResultSignificant()

  out<-RaMP::find_clusters(fishers_df=data,perc_analyte_overlap=as.numeric(input$perc_analyte_overlap),
	min_pathway_tocluster=as.numeric(input$min_pathway_tocluster),
                      perc_pathway_overlap=as.numeric(input$perc_pathway_overlap))
  if(length(unique(out))>1){
    print(paste0(length(out)," clusters found"))
  }else{
    print("Clustering failed")
  }
  return(out)
})

output$cluster_summary_text<-renderText(
  if(as.numeric(input$perc_analyte_overlap) <= 0 || as.numeric(input$perc_analyte_overlap) >= 1 || as.numeric(input$perc_pathway_overlap) <= 0 || as.numeric(input$perc_pathway_overlap) >= 1){
   print("Clustering warning: overlap thresholds must be a percentage greater than 0 and less than 1!")
  }else if(!is.null(cluster_output())){
    if(length(unique(cluster_output()))>1){
    paste0("Fuzzy clustering identified ",length(cluster_output()), " distinct cluster(s) of pathways")
    }else{
      print("Fuzzy clustering algorithm did not identify any clusters. Less stringent thresholds may help in identification, or there may not be enough pathways to cluster.")
    }
  }
)

output$cluster_summary_plot<-renderPlot(
  if(!is.null(cluster_output())&&length(unique(cluster_output()))>1){
    data<-as.numeric(lapply(cluster_output(),length))
    ylim<-c(0, 1.1*max(data))
    xx<- barplot(data, xaxt = 'n', xlab = '', width = 0.85, ylim = ylim, yaxt = 'n',
                 col = "steelblue")
    text(x = xx, y = data, label = data, pos = 3, cex = 0.8)
    axis(1, at=xx, labels=c(1:length(cluster_output())))
    #axis(4, at=c(1,round(max(data)/2),max(data)))
    title("Pathways per cluster", xlab="Cluster #")
    #mtext("# Pathways", side=4, line=-1.5)
  }
)

observe({
  updateSelectInput(session,"show_cluster","Display pathways in cluster:",
                    #choices = as.vector(na.exclude(c("All",ifelse(length(unique(cluster_output()))>1,1:length(cluster_output()),NA),"Did not cluster"),selected = "All")))
                    #choices = c("All",1:length(cluster_output()),"Did not cluster"),selected = "All")
                    choices = as.vector(na.exclude(c("All",ifelse(unique(cluster_output())!="Did not cluster",
			1:length(cluster_output()),NA),"Did not cluster"),selected = "All")))
})

results_fisher_clust <- reactive({
  # Turn this into a function in GeneralFunctions
  cluster_list<-cluster_output()
  if(is.null(fisherTestResult())) {
    data <- data.frame(Query=NA,Freq=NA)
  }else{
    data <- fisherTestResultSignificant()
  }
  if(is.null(cluster_list)){
    return(data)
  }else{
    # Need to remove RaMPID column
    fisher_df<-data[[1]]
    rampids<-fisher_df$pathwayRampId
    fisher_df$pathwayRampId<-NULL

    if(length(cluster_list)>1){
      cluster_assignment<-sapply(rampids,function(x){
        pathway<-x
        clusters<-""
        for(i in 1:length(cluster_list)){
          if(pathway %in% cluster_list[[i]]){
            clusters<-paste0(clusters,i,sep = ", ",collapse = ", ")
          }
        }
        if(clusters!=""){
          clusters=substr(clusters,1,nchar(clusters)-2)
        }else{
          clusters = "Did not cluster"
        }
        return(clusters)
      })
      data[[1]]<-cbind(fisher_df,cluster_assignment)
    }else{
      data[[1]]<-cbind(fisher_df,rep("Did not cluster",times=nrow(fisher_df)))
    }
    #data$Pval <- round(data$Pval,8)
    #data$Adjusted.Pval <- round(data$Adjusted.Pval,8)
    #colnames(data_2)<-c("Pathway Name", "Raw Fisher's P Value","FDR Adjusted P Value","Holm Adjusted P Value",
    #"Source ID","Source DB", "User Analytes in Pathway", "Total Analytes in Pathway", "In Cluster")
    data[[1]]$rampids<-rampids
    return(data)
  }
})

output$results_fisher <- DT::renderDataTable({
  results_fisher<-results_fisher_clust()
  results_fisher<-results_fisher$fishresults
  #results_fisher<-results_fisher[,c("pathwayName","Pval","Pval_FDR","Pval_Holm","pathwaysourceId","pathwaysource",
  #                                  "Num_In_Path","Total_In_Path","cluster_assignment","rampid")]
  results_fisher<-results_fisher[,c(6,1,4,5,7,8,2,3,9,10)]
  colnames(results_fisher)<-c("Pathway Name", "Raw Fisher's P Value","FDR Adjusted P Value","Holm Adjusted P Value",
  "Source ID","Source DB", "User Analytes in Pathway", "Total Analytes in Pathway", "In Cluster","rampids")
  rampids <- results_fisher$rampids
  rampids <- rampids[order(results_fisher[,"In Cluster"])]
  results_fisher$rampids <- NULL
  results_fisher <- results_fisher[order(results_fisher[,"In Cluster"]),]
  cluster_output<-cluster_output()
  if(input$show_cluster=="All"){
    results_fisher
  }else if(input$show_cluster=="Did not cluster"){
    results_fisher[which(results_fisher[,"In Cluster"]=="Did not cluster"),]
  }else{
    results_fisher[which(rampids %in% cluster_output[[as.numeric(input$show_cluster)]]),]
  }
  # data.frame(rep(NA, times = 9))
},rownames = FALSE,filter = "top")

output$fisher_stats_report <- downloadHandler(filename = function(){
  return("fisherText.csv")
},content = function(file){
  #print("Fisher Stats Output has some problems ...")
  rampOut <- results_fisher_clust()
  cluster_output <- cluster_output()
  if(!is.null(rampOut)&&length(unique(cluster_output))>1) {
    cluster_assignment<-apply(rampOut,1,function(x){
      pathway<-x[10]
      clusters<-c()
      for(i in 1:length(cluster_output)){
        if(pathway %in% cluster_output[[i]]){
          clusters<-c(clusters,i)
        }
      }
      return(clusters)
    })
    rampOut<-rampOut[,-10]
    rampOut[,9] <- rep(NA,times = nrow(rampOut))
    colnames(rampOut)[9]<-"In Cluster"
    duplicate_rows<-c()
    for(i in 1:nrow(rampOut)){
      if(is.null(cluster_assignment[[i]])){
        rampOut[i,9]="Did not cluster"
      } else if(length(cluster_assignment[[i]])>1){
        duplicate_rows<-c(duplicate_rows,i)
        for(j in cluster_assignment[[i]]){
          new_row<-c(rampOut[i,1:8],j)
          names(new_row)<-colnames(rampOut)
          rampOut<-rbind(rampOut,new_row)
        }
      }else{
        rampOut[i,9]=cluster_assignment[[i]]
      }
    }
    rampOut<-rampOut[-duplicate_rows,]
    rampOut <- do.call(cbind,rampOut)
  }else{
    rampOut<-rampOut[,-10]
  }
  if(!is.null(rampOut)){
    print(colnames(rampOut))
    rampOut<-rampOut[order(rampOut[,"Holm Adjusted P Value"]),]
    rampOut<-rampOut[!duplicated(rampOut),]
  write.csv(rampOut,file,row.names = FALSE)
  }else{
    write.csv(c("No significant results"),file,row.names = FALSE)
  }
})

