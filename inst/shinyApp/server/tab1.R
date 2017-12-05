
# Panel 1 
observe({
  if (input$buttonstop > 0) {
    stopApp()
  }
})
# output$KWsearch <- renderUI({
#   selectInput("kw_search", "", choices = NULL)
# })

data_stored <- eventReactive(input$subText,{
  progress <- shiny::Progress$new()
  
  on.exit(progress$close())
  
  progress$set(message = "Querying databases to find metabolites...", value = 0)
  progress$inc(0.3,detail = paste("Send Query ..."))
  # rampOut<- rampGenesFromComp(input$kw_search, 90000, input$geneOrComp1)
  rampOut <- RaMP:::rampFastMetaInPath2(input$kw_search)
  
  progress$inc(0.7,detail = paste("Down!"))
  return(rampOut)
})


textOut <- eventReactive(input$subText,{
  if(is.character(data_stored())){
    return("There is no relevant item in database.")
  }
  paste0("Search result: ",nrow(data_stored())," items in databases")
})
output$numOfitems <- renderText({
  validate(
    need(data_stored(),"There is no relevant results in database.")
  )
  textOut()
})

# observe({
#   x <- rampKWsearch(input$singleInput, "analytesynonym")
#   
#   isolate({
#     if (is.null(x)) {
#       x <- character(0)
#     } else {
#       x <- as.matrix(x)
#       x <- matrix(x, ncol = ncol(x), dimnames = NULL)
#     }
#     updateSelectInput(session, "kw_search",
#                       label = "Select from the list",
#                       choices = x, selected = head(x,1)
#     )
#   })
#   
# })

observe({
  choices <- kw_analyte[grepl(input$singleInput,kw_analyte,fixed = T)]
  choices <- choices[order(nchar(choices),choices)]
  if(is.null(choices))
    return(NULL)
  if(length(choices) >10 ){
    choices <- choices[1:10]
  }
  isolate({
    updateSelectInput(session, "kw_search",
                      label = "Select from the list",
                      choices = choices, selected = head(choices,1)
    )
  })
  
})

output$result <- renderDataTable({
  out_src <- data_stored()
  if (is.character(out_src)) {
    return("Not Found")
  }
  
  return(out_src)
  
})

output$report <- downloadHandler(filename = function() {
  paste0(input$kw_search, ".csv")
}, content = function(file) {
  ramp <- data_stored()
  write.csv(ramp, file, append = TRUE, row.names = FALSE)
})

output$preview_tab1 <- renderUI({
  input$subText
  # rampout <- data_stored()
  # 
  # tb_title <- unique(rampout$source_type)
  # tables <- list()
  # for (item in tb_title){
  # }
  tb_output <- rampTablize(data_stored())
  return(div(HTML(unlist(tb_output)),class = "shiny-html-output"))

})


# Panel 2

detector_tab1 <- reactiveValues(num = NULL)

observe({
  input$sub_mul_tab1
  
  detector_tab1$num <- 1
})
data_mul_name_tab1 <- eventReactive(input$sub_mul_tab1,{
  if(is.null(input$sub_mul_tab1))
    return(NULL)
  rampFastPathFromMeta(input$input_mul_tab1)
})
data_mul_file_tab1 <- eventReactive(input$sub_file_tab1,{
  infile <- input$inp_file_tab1
  if (is.null(infile))
    return(NULL)
  
  rampFileOfPathways(infile)
})
observe({
  input$sub_file_tab1
  
  detector_tab1$num <- 2
})


# Download table in a csv file.
output$tab1_mul_report <- downloadHandler(filename = function(){
  if (detector_tab1$num == 1){
    paste0("pathwayOutput.csv")
  } else if (detector_tab1$num == 2){
    infile <- input$inp_file_tab1
    paste0(infile[[1,'name']],"Output",".csv")
  }
},
content = function(file) {
  if (detector_tab1$num == 1){
    rampOut <- data_mul_name_tab1()
  } else if (detector_tab1$num == 2){
    rampOut <- data_mul_file_tab1()
  }
  write.csv(rampOut,file,row.names = FALSE)
}
)

tb_data_tab1 <- reactive({
  if(detector_tab1$num == 1){
    tb <- data_mul_name_tab1()
  } else if (detector_tab1$num == 2) {
    tb <- data_mul_file_tab1()
  }
  tb
})
output$preview_multi_names_tab1 <- renderDataTable({
  if(is.null(detector_tab1$num))
    return("Waiting for input")
  
  tb_data_tab1()
}
)
meta_path_list_tab1 <- reactive({
  if(detector_tab1$num == 1){
    rampGenerateBarPlot(data_mul_name_tab1())
  } else if (detector_tab1$num == 2){
    rampGenerateBarPlot(data_mul_file_tab1())
  }
})


output$tab1_hc_output <- renderHighchart({
  if (is.null(detector_tab1$num) && is.null(input$inp_file_tab1))
    return()
  
  hc_data <- meta_path_list_tab1()
  hc_data <- hc_data[order(sapply(hc_data,nrow),decreasing = TRUE)]
  myClickFunc <- JS("function(event) {Shiny.onInputChange('hcClicked_tab1',event.point.category);}")
  freq <- lapply(hc_data,nrow)
  x_data <- names(freq)
  detail <- sapply(hc_data,as.vector)
  detail <- lapply(detail,paste,collapse = " ")
  names(detail) <- NULL
  
  names(freq) <- NULL
  y_data <- data.frame(y =unlist(freq),detail = unlist(detail))
  hc <- rampHcOutput(x_data,y_data,"column",myClickFunc)          
  return(hc)
})


# interactive plot displays information of a bar.
observe({
  tb <- tb_data_tab1()
  
  isolate({
    pathway_list <- unlist(unique(tb[,1]))
    pathway_list <- sort(pathway_list)
    updateSelectizeInput(session,"pathway_from_hc",
                         choices = pathway_list,
                         server = T,
                         selected = input$hcClicked_tab1$name)
  })
})

output$trying <- renderText({
  input$hcClicked_tab1$name
})

meta_in_pathways <- eventReactive(input$meta_from_path_fire,{
  rampout <- rampFastMetaFromPath(input$pathway_from_hc)
})

output$meta_from_path_preview <- renderDataTable({
  meta_in_pathways()
})
output$meta_from_path_download <- downloadHandler(filename = function(){
    paste0("metabolitesInPathway.csv")
},
content = function(file) {
  rampOut <-  meta_in_pathways()
  write.csv(rampOut,file,row.names = FALSE)
}
)
