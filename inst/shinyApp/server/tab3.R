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

#output$preview_tab3 <- renderUI({
#  input$submit_compName

#  isolate({
#    if(input$synonymOrSource == "synonyms"){
#      tables <- RaMP:::rampTablize(dataInput_name())
#      return(div(HTML(unlist(tables)),class = "shiny-html-output"))
#    } else {
#      return(NULL)
#    }
#  })
#})


# Second Tab
#
#
rea_detector <- reactiveValues(num = NULL)

observe({
  input$sub_mul_tab3

  rea_detector$num <- 1
})


# Batch Query
data_mul_name <- eventReactive(input$sub_mul_tab3,{
  print(input$input_mul_tab3)
  parsedinput <- paste(strsplit(input$input_mul_tab3,"\n")[[1]])
  print(parsedinput)
  RaMP::rampFastPathFromMeta(analytes=parsedinput,
                             NameOrIds=input$NameOrSourcemult,
                             conpass=.conpass,
                             host = .host)
})


data_mul_file <- eventReactive(input$sub_file_tab3,{
  infile <- input$inp_file_tab3
  if (is.null(infile))
    return(NULL)

  RaMP:::rampFileOfPathways(infile,conpass=.conpass,host = .host,NameOrIds=input$NameOrSourcemult)
})

observe({
  input$sub_file_tab3

  rea_detector$num <- 2
})

# Download table in a csv file.
output$tab3_mul_report <- downloadHandler(filename = function(){
  if (rea_detector$num == 1){
    paste0("pathwayOutput.csv")
  } else if (rea_detector$num == 2){
    infile <- input$inp_file_tab3
    paste0(infile[[1,'name']],"Output",".csv")
  }
},
content = function(file) {
  if (rea_detector$num == 1){
    rampOut <- data_mul_name()[,c("pathwayName","pathwaysourceId",
                                  "pathwaysource","commonName")]
    #colnames(rampOut)[4] <- "Analyte"
  } else if (rea_detector$num == 2){
    rampOut <- data_mul_file()[,c("pathwayName","pathwaysourceId",
                                  "pathwaysource","commonName")]
    #colnames(rampOut)[4] <- "Analyte"
  }
  write.csv(rampOut,file,row.names = FALSE)
}
)

output$summary_mulpath_out<- DT::renderDataTable({
  if(is.null(data_mul_name())) {
    out <- data.frame(Query=NA,Freq=NA)
  }
  else {
    temp <- data_mul_name()
    out <- as.data.frame(table(temp$commonName))
    colnames(out)[1] <- "Query"
  }
  out
},rownames=FALSE)

output$preview_multi_names <- DT::renderDataTable({
  if(is.null(rea_detector$num))
    return("Waiting for input")

  if(rea_detector$num == 1){
    tb <- data_mul_name()[,c("pathwayName","pathwaysourceId",
                             "pathwaysource","commonName")]
  } else if (rea_detector$num == 2) {
    tb <- data_mul_file()[,c("pathwayName","pathwaysourceId",
                             "pathwaysource","commonName")]
  }
  #colnames(tb)[4]="Analyte"
  tb
}
,rownames = FALSE)


#fisher_result_bar <- eventReactive(input$runFisher,{
#  if(rea_detector$num == 1){
#    RaMP:::runFisherTest(meta_path_list(),input$analyte_type,
#	total_metabolites=input$total_metabolites, total_genes=input$total_genes,
#	conpass=.conpass)
#  } else if (rea_detector$num == 2){
#    RaMP:::runFisherTest(meta_path_list(),input$analyte_type,
#	total_metabolites=input$total_metabolites, total_genes=input$total_genes,
#	conpass=.conpass)
#  } else {
#    return("No Input")
#  }
#})

#fisher_result_tab3 <- reactive({
#  if (is.null(rea_detector$num)){
#    return()
#  }
#  if (rea_detector$num == 1){
#    RaMP:::rampFisherTest(meta_path_list(),analyte_type=
#length(unique(data_mul_name()$metabolite)),
#	conpass=.conpass)
#  } else if (rea_detector$num == 2){
#    RaMP:::rampFisherTest(meta_path_list(),length(unique(data_mul_file()$metabolite)),
#	conpass=.conpass)
#  }
#})
#output$stats_fisher_tab3 <- renderTable({
#  if (is.null(input$hcClicked))
#    return("Click plots for fisher test...")

# stats <- fisher_result_tab3()
#
# display <- data.frame("Stats result" = c(stats[[input$hcClicked$name]]$p.value,
#                                                stats[[input$hcClicked$name]]$method))
#  stats <- fisher_result_bar()
#  display <- data.frame("Stats Result" = c(paste0("p-value:",stats$p.value),
#                                           paste0("Method: Fisher Exact Test")))
#  return(display)
#},
#rownames =F,
#striped = T)


fisherTestResult <- eventReactive(input$runFisher,{
  print("Generating fisher test result...")
  out <- RaMP::runFisherTest(req(data_mul_name()),analyte_type=input$analyte_type,
                             conpass=.conpass)
  print("Results generated")
  print(paste0("Fisher results size:",nrow(out)))
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
	if(!is.null(data)) {
		parsedinput <- paste(strsplit(input$input_mul_tab3,"\n")[[1]])
		print(paste0("Found ",length(unique(data$commonName))," out of ",
		length(parsedinput)))
	}
})	

fisherTestResultSignificant<-eventReactive(input$runFisher,{
  result<-RaMP::FilterFishersResults(fisherTestResult(),p_fdradj_cutoff=as.numeric(input$p_fdradj_cutoff))
  print(paste0(nrow(result)," significant pathways identified"))
  result
})

cluster_output<-eventReactive(input$runFisher,{
  data <- fisherTestResultSignificant()
  out<-RaMP::find_clusters(data,input$analyte_type, as.numeric(input$perc_analyte_overlap), as.numeric(input$min_pathway_tocluster),
                      as.numeric(input$perc_pathway_overlap),p_cutoff = as.numeric(input$p_fdradj_cutoff))
  if(length(unique(out))>1){
    print(paste0(length(out)," clusters found"))
  }else{
    print("Clustering failed")
  }
  out
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
                    choices = as.vector(na.exclude(c("All",ifelse(unique(cluster_output())!="Did not cluster",1:length(cluster_output()),NA),"Did not cluster"),selected = "All")))
})

total_results_fisher <- eventReactive(input$runFisher,{
  if(is.null(fisherTestResult())) {
    data <- data.frame(Query=NA,Freq=NA)
  }
  data <- fisherTestResultSignificant()
  # Need to remove RaMPID column
  rampids<-data[,9]
  data<-data[,-9]

  cluster_list<-cluster_output()
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
        #&&length(unique(clusters))>1
        clusters=substr(clusters,1,nchar(clusters)-2)
      }else{
        clusters = "Did not cluster"
      }
      return(clusters)
    })
    data_2<-cbind(data,cluster_assignment)
  } else{
    data_2<-cbind(data,rep("Did not cluster",times=nrow(data)))
  }
  #data$Pval <- round(data$Pval,8)
  #data$Adjusted.Pval <- round(data$Adjusted.Pval,8)
  colnames(data_2)<-c("Pathway Name", "Raw Fisher's P Value","FDR Adjusted P Value","Holm Adjusted P Value",
                      "Source ID","Source DB", "User Analytes in Pathway", "Total Analytes in Pathway", "In Cluster")
  data_2<-cbind(data_2,rampids)
})

output$results_fisher <- DT::renderDataTable({
  if(!is.null(total_results_fisher)){
  results_fisher<-total_results_fisher()
  cluster_output<-cluster_output()
  if(input$show_cluster=="All"){
    results_fisher[,-10]
  }else if(input$show_cluster=="Did not cluster"){
    results_fisher[which(results_fisher[,9]=="Did not cluster"),-10]
  }else{
    results_fisher[which(results_fisher[,10] %in% cluster_output[[as.numeric(input$show_cluster)]]),-10]
  }
  }else{
    data.frame(rep(NA, times = 9))
  }
},rownames = FALSE,filter = "top")

output$fisher_stats_report <- downloadHandler(filename = function(){
  return("fisherText.csv")
},content = function(file){
  #print("Fisher Stats Output has some problems ...")
  rampOut <- total_results_fisher()
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
    rampOut<-rampOut[order(rampOut[,"Source DB"]),]
  write.csv(rampOut,file,row.names = FALSE)
  }else{
    write.csv(c("No significant results"),file,row.names = FALSE)
  }
})



#fisherHeatMap <- reactive({
#  print("Generate data for heatmap ...")
#  fisher <- fisherTestResult()
#  if (is.null(fisher))
#    return()
#  hc_data <- data.frame(y = NULL, pathway = NULL)
#
#  for (pathway in names(fisher)){
#    if(fisher[[pathway]]$p.value < as.numeric(input$pvalue_fisher)){
#      print("pathway:")
#      print(pathway)
#      hc_data <- rbind(hc_data,
#                       data.frame(y = fisher[[pathway]]$p.value,
#                                          pathway = pathway))
#    }
#  }
#  return(hc_data)
#})



#output$heatmap_pvalue <- highcharter::renderHighchart({
#  data <- fisherHeatMap()
#  data <- data[order(data$y),]
#  pathway <- as.vector(data$pathway)
#  pvalue <- data$y
#  heatmap_data <- data.frame(v1 = rep(1,length(pathway)),v2 = 1:length(pathway))
#  heatmap_data$value <- pvalue
#  heatmap_data <- list_parse2(heatmap_data)
#
#  fntltp <- highcharter::JS(
#    "function(){
#    return this.series.yAxis.categories[this.point.y] + ' p='+this.point.value;
#  }")
#  hc <- highcharter::highchart() %>%
#    highcharter::hc_chart(type = "heatmap",
#             borderColor = '#ceddff',
#             borderRadius = 10,
#             borderWidth = 2,
#             zoomType = "y",
#             backgroundColor = list(
#               linearGradient = c(0, 0, 500, 500),
#               stops = list(
#                 list(0, 'rgb(255, 255, 255)'),
#                 list(1, 'rgb(219, 228, 252)')
#               ))) %>%
#    hc_title(text = "P-value for Fisher Exact Test") %>%
#    hc_xAxis(categories = c(".",
#                            enable = FALSE),
#             visible = FALSE) %>%
#    hc_yAxis(categories = c("",pathway),
#             visible = TRUE) %>%
#    hc_add_series(name = "pvalue",data = heatmap_data) %>%
#    hc_tooltip(formatter = fntltp, valueDecimals = 2) %>%
#    hc_exporting(enabled = TRUE)
#
#  hc_colorAxis(hc,minColor ="#FFFFFF", maxColor = "#F44242")
#
#})

# Format data from querying database and provide appropriate layout to generate
# bar plot for highcharter.
#meta_path_list <- reactive({
#  if(rea_detector$num == 1){
#      bar_plot_info <-
#        RaMP:::rampGenerateBarPlot(data_mul_name()[,c("pathwayName",
#                "pathwaysourceId","pathwaysource","rampId")])
#  } else if (rea_detector$num == 2){
#      bar_plot_info <-
#        RaMP:::rampGenerateBarPlot(data_mul_file()[,c("pathwayName",
#                "pathwaysourceId","pathwaysource","rampId")])
#  }
#  bar_plot_info <- bar_plot_info[order(sapply(bar_plot_info,nrow),decreasing =TRUE)]
#})

# highchart
# order data in decreasing...
# 12/12 change it to display the log(p) value of each pathways.
#output$tab3_hc_output <- highcharter::renderHighchart({
#  if (is.null(rea_detector$num) && is.null(input$inp_file_tab3))
#    return()
#
#  hc_data <- meta_path_list()
#
#  myClickFunc <- highcharter::JS("function(event) {Shiny.onInputChange('hcClicked',event.point.category);}")
#  freq <- lapply(hc_data,nrow)
#  x_data <- names(freq)
#  detail <- sapply(hc_data,as.vector)
#  detail <- lapply(detail,paste,collapse = " ")
#  names(detail) <- NULL
#
#  names(freq) <- NULL
#  y_data <- data.frame(y = unlist(freq),detail = unlist(detail))
#  hc <- RaMP:::rampHcOutput(x_data,y_data,"column",myClickFunc)
#  return(hc)
#})

## interactive plot displays information of a bar.
#detail_of_bar <- reactive({
#  if (is.null(input$runFisher))
#    return(NULL)
#  output <- meta_path_list()
#  string <- paste(output[[input$hcClicked$name]][[1]],collapse = ' ')
#  output <- paste0("The pathway ",input$hcClicked$name," has ",
#                   length(output[[input$hcClicked$name]][[1]])," metabolites:",string)
#  output <- paste0("Total ",length(unique(names(meta_path_list())))," pathways.",output)
#  return(output)
#})

#output$summary_Fisher <- renderText({
#  detail_of_bar()
#})

#fisher_result_bar <- eventReactive(input$runFisher,{
#  if(rea_detector$num == 1){
#    RaMP:::runFisherTest(meta_path_list(),input$analyte_type,
#       total_metabolites=input$total_metabolites, total_genes=input$total_genes,
#       conpass=.conpass)
#  } else if (rea_detector$num == 2){
#    RaMP:::runFisherTest(meta_path_list(),input$analyte_type,
#       total_metabolites=input$total_metabolites, total_genes=input$total_genes,
#       conpass=.conpass)
#  } else {
#    return("No Input")
#  }
#})
