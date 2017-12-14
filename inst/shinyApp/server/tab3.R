# First tab panel
dataInput_name <- eventReactive(input$submit_compName,{
  progress <- shiny::Progress$new()
  
  on.exit(progress$close())
  
  progress$set(message = "Querying databases to find pathways ...", value = 0)
  progress$inc(0.3,detail = paste("Send Query ..."))

      rampOut <- RaMP::rampFastPathFromMeta(synonym=input$KW_synonym,
                synonymOrIdS=input$synonymOrSource,
                conpass=.conpass)
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
   if(input$synonymOrSource == "synonyms"){
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
   } else if (input$synonymOrSource == "ids"){
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
  RaMP:::rampFastPathFromMeta(synonym=input$input_mul_tab3,
	synonymOrIdS=input$synonymOrSourcemult,
	conpass=.conpass)
})


data_mul_file <- eventReactive(input$sub_file_tab3,{
  infile <- input$inp_file_tab3
  if (is.null(infile))
    return(NULL)
  
  RaMP:::rampFileOfPathways(infile,conpass=.conpass,synonymOrIdS=input$synonymOrSource)
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
		"pathwaysource","metabolite")]
    colnames(rampOut)[4] <- "Analyte"
  } else if (rea_detector$num == 2){
    rampOut <- data_mul_file()[,c("pathwayName","pathwaysourceId",
                "pathwaysource","metabolite")]
    colnames(rampOut)[4] <- "Analyte"
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
    out <- as.data.frame(table(temp$metabolite))
    print(dim(out))
    colnames(out)[1] <- "Query"
  }
  out
},rownames=FALSE)

#output$mulsummary_path <- renderText({
#  summary_mulpath_out()
#})



output$preview_multi_names <- DT::renderDataTable({
  if(is.null(rea_detector$num))
    return("Waiting for input")
  
  if(rea_detector$num == 1){
      tb <- data_mul_name()[,c("pathwayName","pathwaysourceId",
                "pathwaysource","metabolite")]
  } else if (rea_detector$num == 2) {
      tb <- data_mul_file()[,c("pathwayName","pathwaysourceId",
                "pathwaysource","metabolite")]
  }
  colnames(tb)[4]="Analyte"
  tb
}
,rownames = FALSE)

# Format data from querying database and provide appropriate layout to generate
# bar plot for highcharter.
meta_path_list <- reactive({
  if(rea_detector$num == 1){
      bar_plot_info <- RaMP:::rampGenerateBarPlot(data_mul_name()[,c("pathwayName",
		"pathwaysourceId","pathwaysource","metabolite")])
  } else if (rea_detector$num == 2){
      bar_plot_info <- RaMP:::rampGenerateBarPlot(data_mul_file()[,c("pathwayName",
                "pathwaysourceId","pathwaysource","metabolite")])
  }
  bar_plot_info <- bar_plot_info[order(sapply(bar_plot_info,nrow),decreasing =TRUE)]
})

# bar plot
# highchart
# order data in decreasing...
# 12/12 change it to display the log(p) value of each pathways.
output$tab3_hc_output <- highcharter::renderHighchart({
  if (is.null(rea_detector$num) && is.null(input$inp_file_tab3))
    return()
  
  hc_data <- meta_path_list()
  
  myClickFunc <- highcharter::JS("function(event) {Shiny.onInputChange('hcClicked',event.point.category);}")
  freq <- lapply(hc_data,nrow)
  x_data <- names(freq)
  detail <- sapply(hc_data,as.vector)
  detail <- lapply(detail,paste,collapse = " ")
  names(detail) <- NULL
  
  names(freq) <- NULL
  y_data <- data.frame(y = unlist(freq),detail = unlist(detail))
  hc <- RaMP:::rampHcOutput(x_data,y_data,"column",myClickFunc)          
  return(hc)
})


# interactive plot displays information of a bar.
detail_of_bar <- reactive({
  if (is.null(input$hcClicked))
    return(NULL)
  output <- meta_path_list()
  string <- paste(output[[input$hcClicked$name]][[1]],collapse = ' ')
  output <- paste0("The pathway ",input$hcClicked$name," has ",
                   length(output[[input$hcClicked$name]][[1]])," metabolites:",string)
  output <- paste0("Total ",length(unique(names(meta_path_list())))," pathways.",output)
  return(output)
})

output$hc_click_output <- renderText({
  detail_of_bar()
})

fisher_result_bar <- eventReactive(input$hcClicked,{
  if(rea_detector$num == 1){
    RaMP:::rampOneFisherTest(meta_path_list(),input$hcClicked$name,
	length(unique(data_mul_name()$metabolite)),conpass=.conpass)
  } else if (rea_detector$num == 2){
    RaMP:::rampOneFisherTest(meta_path_list(),input$hcClicked$name,
	length(unique(data_mul_file()$metabolite)),conpass=.conpass)
  } else {
    return("No Input")
  }
})
fisher_result_tab3 <- reactive({
  if (is.null(rea_detector$num)){
    return()
  }
  if (rea_detector$num == 1){
    RaMP:::rampFisherTest(meta_path_list(),length(unique(data_mul_name()$metabolite)),
	conpass=.conpass)
  } else if (rea_detector$num == 2){
    RaMP:::rampFisherTest(meta_path_list(),length(unique(data_mul_file()$metabolite)),
	conpass=.conpass)
  }
})
output$stats_fisher_tab3 <- renderTable({
  if (is.null(input$hcClicked))
    return("Click plots for fisher test...")

  # stats <- fisher_result_tab3()
  # 
  # display <- data.frame("Stats result" = c(stats[[input$hcClicked$name]]$p.value,
  #                                                stats[[input$hcClicked$name]]$method))
  stats <- fisher_result_bar()
  display <- data.frame("Stats Result" = c(paste0("p-value:",stats$p.value),
                                           paste0("Method: Fisher Exact Test")))
  return(display)
},
rownames =F,
striped = T)


fisherTestResult <- eventReactive(input$generateFisherTest,{
  print("Generating fisher test result...")
  hc_data <- meta_path_list()
 
  if (rea_detector$num == 1){
    RaMP:::rampFisherTest(hc_data,length(unique(data_mul_name()$metabolite)),
                   FisherPathwayTable = FisherPathwayTable,conpass=.conpass)
  } else if (rea_detector$num == 2){
    RaMP:::rampFisherTest(hc_data,length(unique(data_mul_file()$metabolite)),
                   FisherPathwayTable = FisherPathwayTable,conpass=.conpass)
  }
})
fisherHeatMap <- reactive({
  print("Generate data for heatmap ...")
  fisher <- fisherTestResult()
  if (is.null(fisher))
    return()
  hc_data <- data.frame(y = NULL, pathway = NULL)

  for (pathway in names(fisher)){
    if(fisher[[pathway]]$p.value < as.numeric(input$pvalue_fisher)){
      print("pathway:")
      print(pathway)
      hc_data <- rbind(hc_data,
                       data.frame(y = fisher[[pathway]]$p.value,
                                          pathway = pathway))
    }
  }
  return(hc_data)
})


output$summary_fisher <- DT::renderDataTable({
  
  data <- fisherHeatMap()
  data <- data[,c('pathway','y')]
  colnames(data)[2] <- "pvalue" 
  data$pvalue <- round(data$pvalue,8)
  data
},rownames = FALSE,filter = "top")

output$heatmap_pvalue <- highcharter::renderHighchart({
  data <- fisherHeatMap()
  data <- data[order(data$y),]
  pathway <- as.vector(data$pathway)
  pvalue <- data$y
  heatmap_data <- data.frame(v1 = rep(1,length(pathway)),v2 = 1:length(pathway))
  heatmap_data$value <- pvalue
  heatmap_data <- list_parse2(heatmap_data)
  
  fntltp <- highcharter::JS(
    "function(){
    return this.series.yAxis.categories[this.point.y] + ' p='+this.point.value;
  }")
  hc <- highcharter::highchart() %>%
    highcharter::hc_chart(type = "heatmap",
             borderColor = '#ceddff',
             borderRadius = 10,
             borderWidth = 2,
             zoomType = "y",
             backgroundColor = list(
               linearGradient = c(0, 0, 500, 500),
               stops = list(
                 list(0, 'rgb(255, 255, 255)'),
                 list(1, 'rgb(219, 228, 252)')
               ))) %>%
    hc_title(text = "P-value for Fisher Exact Test") %>%
    hc_xAxis(categories = c(".",
                            enable = FALSE),
             visible = FALSE) %>%
    hc_yAxis(categories = c("",pathway),
             visible = TRUE) %>%
    hc_add_series(name = "pvalue",data = heatmap_data) %>%
    hc_tooltip(formatter = fntltp, valueDecimals = 2) %>%
    hc_exporting(enabled = TRUE)
  
  hc_colorAxis(hc,minColor ="#FFFFFF", maxColor = "#F44242")
  
})

output$stats_report <- downloadHandler(filename = function(){
  return("fisherText.csv")
},content = function(file){
  print("Fisher Stats Output has some problems ...")
  rampOut <- fisherTestResult()
  rampOut <- do.call(cbind,rampOut)
  write.csv(rampOut,file,row.names = FALSE)
})
