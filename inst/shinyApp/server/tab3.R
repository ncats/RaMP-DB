# First tab panel
dataInput_name <- eventReactive(input$submit_compName,{
  progress <- shiny::Progress$new()

  on.exit(progress$close())

  progress$set(message = "Querying databases to find pathways ...", value = 0)
  progress$inc(0.3,detail = paste("Send Query ..."))

  rampOut <- RaMP::getPathwayFromAnalyte(analytes=input$KW_synonym,
                                        NameOrIds=input$NameOrId,
                                        conpass=.conpass,
                                        host = .host, dbname = .dbname, username = .username)
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
    #choices <- kw_analyte[grepl(input$compName,kw_analyte,ignore.case=TRUE)]
    if(input$compName=="") {
	choices <- ""
    } else {
    	choices <- agrep(input$compName,kw_analyte,value=TRUE,ignore.case=TRUE)
	choices <- choices[order(nchar(choices),choices)]
    }
    if(is.null(choices))
      return(NULL)
    #if(length(choices) >10 ){
    #  choices <- choices[1:10]
    #}
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
},rownames = FALSE)

output$comp_report <- downloadHandler(filename = function() {
  paste0(input$KW_synonym, ".csv")
},
content = function(file) {
  rampOut <- dataInput_name()
  rampOut <- data.frame(rampOut)
  write.csv(rampOut, file, row.names = FALSE)
})

############
# Second Tab
#############
data_mul_name <- eventReactive(input$sub_mul_tab3,{
  #print(input$input_mul_tab3)
  parsedinput <- paste(strsplit(input$input_mul_tab3,"\n")[[1]])
  print(parsedinput)
  if(length(parsedinput)==0) {metabsearch=NULL} else{
  	metabsearch <- RaMP::getPathwayFromAnalyte(analytes=parsedinput,
                             NameOrIds=input$NameOrSourcemult,
                             conpass=.conpass,
                             host = .host,
                             dbname = .dbname, username = .username)
	  print(input$input_mul_tab3_genes)
  }

  parsedinputg <- paste(strsplit(input$input_mul_tab3_genes,"\n")[[1]])
  print(parsedinputg)
  if(length(parsedinputg)==0) {genesearch=NULL} else{
	  genesearch <- RaMP::getPathwayFromAnalyte(analytes=parsedinputg,
                             NameOrIds=input$NameOrSourcemult_genes,
                             conpass=.conpass,
                             host = .host,
                             dbname = .dbname, username = .username)
  }
  print(paste("metabsearch: ",ncol(metabsearch)))
  print(paste("genesearch: ",ncol(genesearch)))
  print(paste0("DIM of data_mul_name",nrow(rbind(metabsearch,genesearch))))
  rbind(metabsearch,genesearch)
})

# Download table in a csv file.
output$tab3_mul_report <- downloadHandler(filename = function(){
      paste0("pathwayFromMetabolitesOutput.csv")
}, content = function(file) {
  rampOut <- data_mul_name()[,c("pathwayName","pathwaysourceId",
                                    "pathwaysource","commonName")]
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
    colnames(out) <- c("Query","Num_Pathways")
  out
},rownames=FALSE)

output$preview_multi_names <- DT::renderDataTable({
  if(is.null(data_mul_name())) {
    return("No input found")
  } else {

    tb <- data_mul_name()[,c("pathwayName","pathwaysourceId",
                             "pathwaysource","commonName")]
  return(tb)
 }
}
,rownames = FALSE)


fisherTestResult <- eventReactive(input$runFisher,{

    out <- RaMP::runCombinedFisherTest(req(data_mul_name()),
                               conpass=.conpass,
                               dbname = .dbname, username = .username, host = .host)
    print("Results generated")
    print(paste0("Fisher results size:",nrow(out[[1]])))
  out
})

output$summary_fisher <- DT::renderDataTable({
  if(!is.null(fisherTestResult())) {
    data <- fisherTestResult()
    out=as.data.frame(table(data$fishresults$pathwaysource))
    colnames(out)[1]="Pathway_Source"
  } else {
    out <- data.frame(Pathway_Source=NA, Freq=NA)
  }
  out
},rownames=FALSE,filter="top")


output$num_mapped_namesids <- renderText({
	  data <- data_mul_name()
	  inputlist <- input$input_mul_tab3
	  inputlist2 <-input$input_mul_tab3_genes
	  inputsize=0
	  if(!is.null(inputlist) && !is.null(inputlist2)) {
	  	#inputsize <- length(inputlist)+length(inputlist2)
		inputsize <- length(strsplit(inputlist,"\n")[[1]])+length(strsplit(inputlist2,"\n")[[1]])
	  } else if (!is.null(inputlist) && is.null(inputlist2)) {
		#inputsize <- length(inputlist)
		inputsize <- length(strsplit(inputlist,"\n")[[1]])
          } else if (is.null(inputlist) && !is.null(inputlist2)) {
		#inputsize <- length(inputlist2)
		inputsize <- length(strsplit(inputlist2,"\n")[[1]])
	  }
	if(!is.null(data) && inputsize>0) {
		print(paste0("Found ",length(unique(data$commonName))," out of ",
		inputsize))
	}
})

output$fishersProgress<-renderText(
  if(!is.null(fisherTestResult())){
    print("Hit 'Filter and Cluster Results' to view below and download")
  }
)

fisherTestResultSignificant<-eventReactive(input$runClustering,{
  if(!is.null(fisherTestResult())){
    result<-RaMP::FilterFishersResults(fisherTestResult(),p_holmadj_cutoff=as.numeric(input$p_holmadj_cutoff))
    return(result)
  }else{
    return(NULL)
  }
})

cluster_output<-eventReactive(c(input$runFisher,input$runClustering),{
  data <- fisherTestResultSignificant()
  out<-RaMP::findCluster(fishers_df=data,perc_analyte_overlap=as.numeric(input$perc_analyte_overlap),
	min_pathway_tocluster=as.numeric(input$min_pathway_tocluster),
                      perc_pathway_overlap=as.numeric(input$perc_pathway_overlap))
  cluster_list<-out$cluster_list
  if(length(unique(cluster_list))>1){
    print(paste0(length(cluster_list)," clusters found"))
  }else{
    print("Clustering failed")
  }
  return(out)
})

cluster_list<-reactive({
  out<-cluster_output()
  if(!is.null(out)){
    cluster_list<-out$cluster_list
  } else {
    return('Nothing found based on given filters.')
  }
})

output$cluster_summary_text<-renderText(
  #out<-cluster_output()
  if(as.numeric(input$perc_analyte_overlap) <= 0 || as.numeric(input$perc_analyte_overlap) >= 1 || as.numeric(input$perc_pathway_overlap) <= 0 || as.numeric(input$perc_pathway_overlap) >= 1){
   print("Clustering warning: overlap thresholds must be a percentage greater than 0 and less than 1!")
  }else if(!is.null(cluster_output())){
    #cluster_list<-out$cluster_list
    if(length(unique(cluster_list()))>1){
    paste0("Fuzzy clustering identified ",length(cluster_list()), " distinct cluster(s) of pathways")
    }else{
      print("Fuzzy clustering algorithm did not identify any clusters. Less stringent thresholds may help in identification, or there may not be enough pathways to cluster.")
    }
  }
)

output$cluster_summary_plot<-renderPlot({
  out<-cluster_output()
  if(!is.null(out)&&length(unique(cluster_list()))>1){
    data<-as.numeric(lapply(cluster_list(),length))
    ylim<-c(0, 1.1*max(data))
    xx<- barplot(data, xaxt = 'n', xlab = '', width = 0.85, ylim = ylim, yaxt = 'n',
                 col = "steelblue")
    text(x = xx, y = data, label = data, pos = 3, cex = 0.8)
    axis(1, at=xx, labels=c(1:length(cluster_list())))
    #axis(4, at=c(1,round(max(data)/2),max(data)))
    title("Pathways per cluster", xlab="Cluster #")
    #mtext("# Pathways", side=4, line=-1.5)
  }
})

observe({
  out<-cluster_output()
  cluster_list<-out$cluster_list
  updateSelectInput(session,"show_cluster","Display pathways in cluster:",
                    choices = as.vector(na.exclude(c("All",ifelse(unique(cluster_list())!="Did not cluster",
			1:length(cluster_list()),NA),"Did not cluster"),selected = "All")))
})

results_fisher_clust <- reactive({
  cluster_out<-cluster_output()
  if(is.null(fisherTestResult())) {
    data <- data.frame(Query=NA,Freq=NA)
  }else{
    data <- fisherTestResultSignificant()
  }
  if(is.null(cluster_out)){
    return(data)
  }else{
    return(cluster_out)
    # # Need to remove RaMPID column
    # fisher_df<-data[[1]]
    # rampids<-fisher_df$pathwayRampId
    # fisher_df$pathwayRampId<-NULL
    #
    # if(length(cluster_list)>1){
    #   cluster_assignment<-sapply(rampids,function(x){
    #     pathway<-x
    #     clusters<-""
    #     for(i in 1:length(cluster_list)){
    #       if(pathway %in% cluster_list[[i]]){
    #         clusters<-paste0(clusters,i,sep = ", ",collapse = ", ")
    #       }
    #     }
    #     if(clusters!=""){
    #       clusters=substr(clusters,1,nchar(clusters)-2)
    #     }else{
    #       clusters = "Did not cluster"
    #     }
    #     return(clusters)
    #   })
    #   data[[1]]<-cbind(fisher_df,cluster_assignment)
    # }else{
    #   data[[1]]<-cbind(fisher_df,rep("Did not cluster",times=nrow(fisher_df)))
    # }
    # #data$Pval <- round(data$Pval,8)
    # #data$Adjusted.Pval <- round(data$Adjusted.Pval,8)
    # #colnames(data_2)<-c("Pathway Name", "Raw Fisher's P Value","FDR Adjusted P Value","Holm Adjusted P Value",
    # #"Source ID","Source DB", "User Analytes in Pathway", "Total Analytes in Pathway", "In Cluster")
    # data[[1]]$rampids<-rampids
    # return(data)
  }
})

output$results_fisher <- DT::renderDataTable({
  results_fisher_total<-results_fisher_clust()
  results_fisher<-results_fisher_total$fishresults
  if(nrow(results_fisher)==0){
    #return(data.frame(rep(NA, times = 9)))
    stop("No significant pathways identified. Your input analytes may be too small (e.g. number of analytes in each pathway is < 2)")
  }else{
    if(results_fisher_total$analyte_type=="both"){
      results_fisher<-results_fisher[,c("pathwayName","Num_In_Path.Metab","Total_In_Path.Metab",
                                        "Num_In_Path.Gene","Total_In_Path.Gene", "Pval_combined",
                                        "Pval_combined_FDR","Pval_combined_Holm","pathwaysourceId","pathwaysource",
                                        "cluster_assignment","rampids")]
      # Filtered:
      #"Pval.Metab","Pval.Gene",

      colnames(results_fisher)<-c("Pathway Name", "User Metabolites in Pathway",
                                  "Total Metabolites in Pathway","User Genes in Pathway",
                                  "Total Genes in Pathway","Raw Fisher's P Value (Combined)","FDR Adjusted P Value (Combined)",
                                  "Holm Adjusted P Value (Combined)","Source ID","Source DB","In Cluster","rampids")
      # Filtered:
      # "Raw Fisher's P Value (Metabolites)","Raw Fisher's P Value (Genes)",
    }else{
      #print(colnames(results_fisher))
      results_fisher<-results_fisher[,c("pathwayName","Pval","Pval_FDR","Pval_Holm","pathwaysourceId","pathwaysource",
                                        "Num_In_Path","Total_In_Path","cluster_assignment","rampids")]
      colnames(results_fisher)<-c("Pathway Name", "Raw Fisher's P Value","FDR Adjusted P Value","Holm Adjusted P Value",
                                  "Source ID","Source DB", "User Analytes in Pathway", "Total Analytes in Pathway",
                                  "In Cluster","rampids")
    }
    rampids <- results_fisher$rampids
    rampids <- rampids[order(results_fisher[,"In Cluster"])]
    results_fisher$rampids <- NULL
    results_fisher <- results_fisher[order(results_fisher[,"In Cluster"]),]
    cluster_output<-cluster_list()
    if(input$show_cluster=="All"){
      results_fisher
    }else if(input$show_cluster=="Did not cluster"){
      results_fisher[which(results_fisher[,"In Cluster"]=="Did not cluster"),]
    }else{
      results_fisher[which(rampids %in% cluster_output[[as.numeric(input$show_cluster)]]),]
    }
  }
},rownames = FALSE,filter = "top")

output$fisher_stats_report <- downloadHandler(filename = function(){
  return("fisherText.csv")
},content = function(file){
  out <- results_fisher_clust()
  rampOut <- out$fishresults
  cluster_output <- cluster_list()
  if(!is.null(rampOut)) {
    if(out$analyte_type=="both"){
      rampOut<-rampOut[,c("pathwayName","Pval.Metab","Num_In_Path.Metab","Total_In_Path.Metab",
                                        "Pval.Gene", "Num_In_Path.Gene","Total_In_Path.Gene", "Pval_combined",
                                        "Pval_combined_FDR","Pval_combined_Holm","pathwaysourceId","pathwaysource",
                                        "cluster_assignment","rampids")]

      colnames(rampOut)<-c("Pathway Name", "Raw Fisher's P Value (Metabolites)","User Metabolites in Pathway",
                                  "Total Metabolites in Pathway","Raw Fisher's P Value (Genes)","User Genes in Pathway",
                                  "Total Genes in Pathway","Raw Fisher's P Value (Combined)","FDR Adjusted P Value (Combined)",
                                  "Holm Adjusted P Value (Combined)","Source ID","Source DB","In Cluster","rampids")
      rampOut<-rampOut[order(rampOut[,"Holm Adjusted P Value (Combined)"]),]
    }else{
      results_fisher<-rampOut[,c("pathwayName","Pval","Pval_FDR","Pval_Holm","pathwaysourceId","pathwaysource",
                                        "Num_In_Path","Total_In_Path","cluster_assignment","rampids")]
      colnames(rampOut)<-c("Pathway Name", "Raw Fisher's P Value","FDR Adjusted P Value","Holm Adjusted P Value",
                                  "Source ID","Source DB", "User Analytes in Pathway", "Total Analytes in Pathway",
                                  "In Cluster","rampids")
      rampOut<-rampOut[order(rampOut[,"Holm Adjusted P Value"]),]
    }
  write.csv(rampOut,file,row.names = FALSE)
  }else{
    write.csv(c("No significant results"),file,row.names = FALSE)
  }
})

