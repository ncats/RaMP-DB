dataInput_onto <- eventReactive(input$subText_onto,{
  progress <- shiny::Progress$new()

  on.exit(progress$close())

  progress$set(message = "Querying databases to find pathways ...", value = 0)
  progress$inc(0.3,detail = paste("Send Query ..."))

  # rampOut <- rampOntoOut(input$KW_onto, 99999)
  if(input$metaOrOnto %in% c('ids','name')){
    rampOut <- RaMP:::getOntoFromMeta(input$KW_onto,
                                          conpass =.conpass,
                                          host = .host,
                                          NameOrIds = input$metaOrOnto)
  } else if(input$metaOrOnto == 'ontology'){
    rampOut <- RaMP:::getMetaFromOnto(input$KW_onto,
                                          conpass = .conpass,
                                          host = .host)
  }
  progress$inc(0.7,detail = paste("Done!"))
  return (rampOut)
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
  } else if (input$metaOrOnto == 'ids'){
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
  df <- dataInput_onto()
  if(input$metaOrOnto %in% c('ids','name')){
    df <- transform(df,sourceId = paste(IDtype,sourceId,sep = ':'))
    print(colnames(df))
    df <- df[c('Metabolites','Ontology','biofluidORcellular','sourceId')]
    df2 <- aggregate(.~Metabolites+biofluidORcellular,df,FUN = function(x){
      paste(unique(x),collapse = ',')
    })
    colnames(df2) <- c('Metabolites','Ontology_type','Ontology','metabolites_source')
    return(df2)
  } else if(input$metaOrOnto == 'ontology'){
    df <- transform(df,sourceId = paste(IDtype,sourceId,sep = ':'))
    df <- df[c('Metabolites','Ontology','biofluidORcellular','sourceId')]
    df2 <- aggregate(.~Metabolites+biofluidORcellular,df,FUN = function(x){
      paste(unique(x),collapse = ',')
    })
    colnames(df2) <- c('Metabolites','Ontology_type','Ontology','metabolites_source')
    return(df2)
  }

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

data_mul_tab5 <- eventReactive(input$sub_mul_tab5,{
  if(input$input_categories_tab5 %in% c('ids','name')){
    RaMP:::getOntoFromMeta(input$input_mul_tab5,
                                conpass = .conpass,
                                host =.host,
                                username = .username,
                                dbname = .dbname,
                                NameOrIds =input$input_categories_tab5)
  } else if(input$input_categories_tab5 == 'ontology'){
    RaMP:::getMetaFromOnto(input$input_mul_tab5,
                                conpass = .conpass,
                                host = .host,
                                username = .username,
                                dbname = .dbname)
  }
})

output$preview_multi_names_tab5 <- DT::renderDataTable({
  if(is.null(clicked$content))
    return(NULL)
  df <- data_mul_tab5()
  df <- transform(df,sourceId = paste(IDtype,sourceId,sep = ':'))
  df <- df[c('Metabolites','Ontology','biofluidORcellular','sourceId')]
  df2 <- aggregate(.~Metabolites+biofluidORcellular,df,FUN = function(x){
    paste(unique(x),collapse = ',')
  })
  colnames(df2) <- c('Metabolites','Ontology_type','Ontology','metabolites_source')
  return(df2)
},rownames = F)

output$tab5_mul_report <- downloadHandler(filename = function() {
  if(is.null(clicked$content))
    return(NULL)
  return(paste0(input$input_categories_tab5,'report.csv'))
}, content = function(file) {
  if(is.null(clicked$content))
    return(NULL)
  rampout <- data_mul_tab5()

  write.csv(rampout, file, row.names = FALSE)
})
