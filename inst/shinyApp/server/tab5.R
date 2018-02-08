dataInput_onto <- eventReactive(input$subText_onto,{
  progress <- shiny::Progress$new()
  
  on.exit(progress$close())
  
  progress$set(message = "Querying databases to find pathways ...", value = 0)
  progress$inc(0.3,detail = paste("Send Query ..."))
  
  # rampOut <- rampOntoOut(input$KW_onto, 99999)
  rampOut <- rampFastBiofluid(input$KW_onto,input$metaOrOnto)
  progress$inc(0.7,detail = paste("Done!"))
  return (rampOut[,1:3])
})
summary_onto_out<- eventReactive(input$subText_onto,{
  if (!is.null(nrow(dataInput_onto()))){
    return (paste0("There are(is) ",nrow(dataInput_onto())," relevent items in databases."))
  } else{
    return ("Given metabolites have no search result.")
  }
})


observe({
  if (input$metaOrOnto == "ontology") {
    # x <- rampKWsearch(input$ontoInput, "ontology")
    choices <- kw_biofluid[grepl(input$ontoInput,kw_biofluid,fixed = T)]
    choices <- choices[order(nchar(choices),choices)]
  } else {
    # x <- rampKWsearch(input$ontoInput, "analytesynonym")
    choices <- kw_analyte[grepl(input$ontoInput,kw_analyte)]
    choices <- choices[order(nchar(choices),choices)]
  }
  if(length(choices) == length(kw_pathway))
    return(NULL)
  isolate({
    if(is.null(choices))
      return(NULL)
    if(length(choices) >10 ){
      choices <- choices[1:10]
    }
    
    updateSelectInput(session,
                      "KW_onto", 
                      label = "Select from list", 
                      choices = choices, 
                      selected = head(choices,1))
  })
})

output$summary_onto <- renderText(
  {
    summary_onto_out()
  }
)
output$result_onto <- DT::renderDataTable({
  dataInput_onto()
}, rownames = FALSE)

output$report_onto <- downloadHandler(filename = function() {
  paste0(input$KW_onto, ".csv")
}, content = function(file) {
  rampOut <- dataInput_onto()
  rampOut <- data.frame(rampOut)
  write.csv(rampOut, file, row.names = FALSE, sep = ",")
})



output$preview_tab5 <- renderUI({
  input$subText_onto
  
  isolate({
    rampout <- dataInput_onto()
    if(length(rampout) == 0){
      return(HTML("<strong> No Summary due to small dataset</strong>"))
    }
    tables <- rampTablize(rampout)
    return(div(HTML(unlist(tables)),class = "shiny-html-output"))
  })
})

# second sub Tab
detector_tab5 <- reactiveValues(num = NULL)

observe({
  input$sub_mul_tab5
  
  detector_tab5$num <- 1
})
data_mul_name_tab5 <- eventReactive(input$sub_mul_tab5,{
  rampFastBiofluid(input$input_mul_tab5)
})
data_mul_file_tab5 <- eventReactive(input$sub_file_tab5,{
  infile <- input$inp_file_tab5
  if (is.null(infile))
    return(NULL)
  
  rampOut <- rampFileOfBiofluid(infile)
})
observe({
  input$sub_file_tab5
  
  detector_tab5$num <- 2
})

output$preview_multi_names_tab5_A <- DT::renderDataTable({
  if(detector_tab5$num == 1){
    rampout <- data_mul_name_tab5()
    rampout[['analyte']]
  } else if (detector_tab5$num == 2){
    rampout <- data_mul_file_tab5()
    rampout[['analyte']]
  } else {
    return(NULL)
  }
},rownames = F)

output$preview_multi_names_tab5_B <- DT::renderDataTable({
  if(detector_tab5$num == 1){
    rampout <- data_mul_name_tab5()
    rampout[['biofluid']]
  } else if (detector_tab5$num == 2) {
    rampout <- data_mul_file_tab5()
    rampout[['biofluid']]
  } else {
    return(NULL)
  }
},rownames = F)

output$tab5_mul_report_A <- downloadHandler(filename = function() {
  if(detector_tab5$num == 1){
    paste0("textAreaInput", ".csv")
  } else if (detector_tab5$num == 2) {
    paste0("fileInput",".csv")
  }
}, content = function(file) {
  if(detector_tab5$num == 1){
    rampout <- data_mul_name_tab5()
    rampout <- rampout[['analyte']]
  } else if (detector_tab5$num == 2){
    rampout <- data_mul_file_tab5()
    rampout <- rampout[['analyte']]
  } else {
    return(NULL)
  }
  
  write.csv(rampout, file, row.names = FALSE)
})

output$tab5_mul_report_B <- downloadHandler(filename = function() {
  if(detector_tab5$num == 1){
    paste0("textAreaInput", ".csv")
  } else if (detector_tab5$num == 2) {
    paste0("fileInput",".csv")
  }
}, content = function(file) {
  if(detector_tab5$num == 1){
    rampout <- data_mul_name_tab5()
    rampout <- rampout[['biofluid']]
  } else if (detector_tab5$num == 2){
    rampout <- data_mul_file_tab5()
    rampout <- rampout[['biofluid']]
  } else {
    return(NULL)
  }
  
  write.csv(rampout, file, row.names = FALSE)
})
