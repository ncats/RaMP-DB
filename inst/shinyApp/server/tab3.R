# First tab panel
dataInput_name <- eventReactive(input$submit_compName,{
  progress <- shiny::Progress$new()
  
  on.exit(progress$close())
  
  progress$set(message = "Querying databases to find pathways ...", value = 0)
  progress$inc(0.3,detail = paste("Send Query ..."))
  
  rampOut <- rampPathFromMeta(input$KW_synonym, 99999)
  progress$inc(0.7,detail = paste("Done!"))
  return (rampOut)
})

summary_path_out<- eventReactive(input$submit_compName,{
  if (!is.null(nrow(dataInput_name()))){
    return (paste0("There are(is) ",nrow(dataInput_name())," relevent items in databases."))
  } else{
    return ("Given metabolites have no search result.")
  }
})

output$summary_path <- renderText({
  summary_path_out()
})
# output$KWsearchComp <- renderUI({
#   selectInput("KW_synonym", "", choices = NULL)
# })


observe({

  choices <- kw_analyte[grepl(input$compName,kw_analyte,fixed = T)]
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
})



output$result3 <- DT::renderDataTable({
  out_stc <- dataInput_name()
})

output$comp_report <- downloadHandler(filename = function() {
  paste0(input$KW_synonym, ".csv")
},
content = function(file) {
  rampOut <- dataInput_name()
  rampOut <- data.frame(rampOut)
  write.csv(rampOut, file, row.names = FALSE, sep = ",")
})

output$preview_tab3 <- renderUI({
  input$submit_compName
  
  isolate({
    tables <- rampTablize(dataInput_name())
    return(div(HTML(unlist(tables)),class = "shiny-html-output"))
  })
})


# Second Tab
#
#
rea_detector <- reactiveValues(num = NULL)

observe({
  input$sub_mul_tab3
  
  rea_detector$num <- 1
})
data_mul_name <- eventReactive(input$sub_mul_tab3,{
  rampFastPathFromMeta(input$input_mul_tab3)
})
data_mul_file <- eventReactive(input$sub_file_tab3,{
  infile <- input$inp_file_tab3
  if (is.null(infile))
    return(NULL)
  
  rampFileOfPathways(infile)
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
    rampOut <- data_mul_name()
  } else if (rea_detector$num == 2){
    rampOut <- data_mul_file()
  }
  write.csv(rampOut,file,row.names = FALSE)
}
)
output$preview_multi_names <- DT::renderDataTable({
  if(is.null(rea_detector$num))
    return("Waiting for input")
  
  if(rea_detector$num == 1){
      tb <- data_mul_name()
  } else if (rea_detector$num == 2) {
      tb <- data_mul_file()
  }
  tb
}
,rownames = FALSE)
meta_path_list <- reactive({
  if(rea_detector$num == 1){
      bar_plot_info <- rampGenerateBarPlot(data_mul_name())
  } else if (rea_detector$num == 2){
      bar_plot_info <- rampGenerateBarPlot(data_mul_file())
  }
  bar_plot_info <- bar_plot_info[order(sapply(bar_plot_info,nrow),decreasing =TRUE)]
})
# bar plot
# highchart
# order data in decreasing...
output$tab3_hc_output <- renderHighchart({
  if (is.null(rea_detector$num) && is.null(input$inp_file_tab3))
    return()
  
  hc_data <- meta_path_list()
  
  myClickFunc <- JS("function(event) {Shiny.onInputChange('hcClicked',event.point.category);}")
  freq <- lapply(hc_data,nrow)
  x_data <- names(freq)
  detail <- sapply(hc_data,as.vector)
  detail <- lapply(detail,paste,collapse = " ")
  names(detail) <- NULL
  
  names(freq) <- NULL
  y_data <- data.frame(y = unlist(freq),detail = unlist(detail))
  hc <- rampHcOutput(x_data,y_data,"column",myClickFunc)          
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
    rampOneFisherTest(meta_path_list(),input$hcClicked$name,length(unique(data_mul_name()$metabolite)))
  } else if (rea_detector$num == 2){
    rampOneFisherTest(meta_path_list(),input$hcClicked$name,length(unique(data_mul_file()$metabolite)))
  } else {
    return("No Input")
  }
})
fisher_result_tab3 <- reactive({
  if (is.null(rea_detector$num)){
    return()
  }
  if (rea_detector$num == 1){
    rampFisherTest(meta_path_list(),length(unique(data_mul_name()$metabolite)))
  } else if (rea_detector$num == 2){
    rampFisherTest(meta_path_list(),length(unique(data_mul_file()$metabolite)))
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
    rampFisherTest(hc_data,length(unique(data_mul_name()$metabolite)),
                   FisherPathwayTable = FisherPathwayTable)
  } else if (rea_detector$num == 2){
    rampFisherTest(hc_data,length(unique(data_mul_file()$metabolite)),
                   FisherPathwayTable = FisherPathwayTable)
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

output$heatmap_pvalue <- renderHighchart({
  data <- fisherHeatMap()
  data <- data[order(data$y),]
  pathway <- as.vector(data$pathway)
  pvalue <- data$y
  heatmap_data <- data.frame(v1 = rep(1,length(pathway)),v2 = 1:length(pathway))
  heatmap_data$value <- pvalue
  heatmap_data <- list_parse2(heatmap_data)
  
  fntltp <- JS(
    "function(){
    return this.series.yAxis.categories[this.point.y] + ' p='+this.point.value;
  }")
  hc <- highchart() %>%
    hc_chart(type = "heatmap",
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
