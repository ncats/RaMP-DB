dataInput_onto <- eventReactive(input$subText_onto,{
  progress <- shiny::Progress$new()
  
  on.exit(progress$close())
  
  progress$set(message = "Querying databases to find pathways ...", value = 0)
  progress$inc(0.3,detail = paste("Send Query ..."))
  
  # rampOut <- rampOntoOut(input$KW_onto, 99999)
  if(input$metaOrOnto %in% c('source','name')){
    rampOut <- RaMP:::rampFastOntoFromMeta(input$KW_onto,
                                          conpass =.conpass,
                                          host = .host,
                                          sourceOrName = input$metaOrOnto)
  } else if(input$metaOrOnto == 'ontology'){
    rampOut <- RaMP:::rampFastMetaFromOnto(input$KW_onto,
                                          conpass = .conpass,
                                          host = .host)
  }
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
  if(input$metaOrOnto == 'ontology'){
    choices <- kw_biofluid[grepl(input$ontoInput,kw_biofluid,fixed = T)]
    choices <- choices[order(nchar(choices),choices)]
    choices <- c(choices,kw_biofluid[!(kw_biofluid %in% choices)])
    updateSelectInput(session,
                      "KW_onto", 
                      label = "Select from list", 
                      choices = choices, 
                      selected = head(choices,1))
  }
})
observe({
  if(input$metaOrOnto == 'ontology'){
    choices <- kw_biofluid[grepl(input$ontoInput,kw_biofluid,fixed = T)]
    choices <- choices[order(nchar(choices),choices)]
    choices <- c(choices,kw_biofluid[!(kw_biofluid %in% choices)])
    
  }
  else if (input$metaOrOnto == 'name') {
    # x <- rampKWsearch(input$ontoInput, "analytesynonym")
    choices <- kw_analyte[grepl(input$ontoInput,kw_analyte,
                                fixed = T)]
    choices <- choices[order(nchar(choices),choices)]
  } else if (input$metaOrOnto == 'source'){
    choices <- kw_source[grepl(input$ontoInput,kw_source,
                               fixed = T)]
    choices <- choices[order(nchar(choices),choices)]
  } 

  isolate({
    if(is.null(choices))
      return(NULL)
    if(length(choices) >10 & input$metaOrOnto != 'ontology'){
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


# second sub Tab
clicked <- reactiveValues(content = NULL)
# Identify which button user clicked
observe({
  input$sub_mul_tab5
  
  clicked$content <- 'name'
  
})
observe({
  input$sub_mul_tab5_sourceid
  
  clicked$content <- 'source'
  
})
observe({
  input$sub_mul_tab5_biofluid
  
  clicked$content <- 'ontology'
  
})

data_mul_name_tab5 <- eventReactive(input$sub_mul_tab5,{
  RaMP:::rampFastOntoFromMeta(input$input_mul_tab5,
                              conpass = .conpass,
                              host =.host,
                              username = .username,
                              dbname = .dbname,
                              sourceOrName ='name')
})
data_mul_source_tab5 <- eventReactive(input$sub_mul_tab5_sourceid,{
  RaMP:::rampFastOntoFromMeta(input$input_mul_tab5_sourceid,
                              conpass = .conpass,
                              host =.host,
                              username = .username,
                              dbname = .dbname,
                              sourceOrName ='source')
})
data_mul_biofluid_tab5 <- eventReactive(input$sub_mul_tab5_biofluid,{
  RaMP:::rampFastMetaFromOnto(input$input_mul_tab5_biofluid,
                              conpass = .conpass,
                              host =.host,
                              username = .username,
                              dbname = .dbname
                              )
})
# data_mul_file_tab5 <- eventReactive(input$sub_file_tab5,{
#   infile <- input$inp_file_tab5
#   if (is.null(infile))
#     return(NULL)
#   
#   rampOut <- rampFileOfBiofluid(infile)
# })
# observe({
#   input$sub_file_tab5
#   
#   detector_tab5$num <- 2
# })

output$preview_multi_names_tab5 <- DT::renderDataTable({
  if(is.null(clicked$content))
    return(NULL)
  if(clicked$content == 'source'){
    data_mul_source_tab5()
  } else if (clicked$content == 'name'){
    data_mul_name_tab5()
  } else if (clicked$content == 'ontology' ){
    data_mul_biofluid_tab5()
  }
},rownames = F)

output$tab5_mul_report <- downloadHandler(filename = function() {
  if(detector_tab5$num == 1){
    paste0("textAreaInput", ".csv")
  } else if (detector_tab5$num == 2) {
    paste0("fileInput",".csv")
  }
}, content = function(file) {
  if(clicked$content == 'name'){
    rampout <- data_mul_name_tab5()
  } else if (clicked$content == 'source'){
    rampout <- data_mul_source_tab5() 
    
  } else if (clicked$conttent == 'ontology') {
    rampout <- data_mul_biofluid_tab5()
  } else {
    return(NULL)
  }
  
  write.csv(rampout, file, row.names = FALSE)
})