observe({
  options <- input$givenDataType
  
  if(options == "metabolites"||options == "genes"){
    choices <- kw_analyte[grepl(input$kwInput,kw_analyte)]
    choices <- choices[order(nchar(choices),choices)]
    if(length(choices) == length(kw_analyte)|| length(choices) == 0)
      return(NULL)
  } else if(options == "pathways"){
    choices <- kw_pathway[grepl(input$kwInput,kw_pathway)]
    choices <- choices[order(nchar(choices),choices)]
    if(length(choices) == length(kw_pathway)|| length(choices) == 0)
      return(NULL)
  }else if(options == "biofluid"){
    choices <- kw_biofluid[grepl(input$kwInput,kw_biofluid)]
    choices <- choices[order(nchar(choices),choices)]
    if(length(choices) == length(kw_biofluid)|| length(choices) == 0)
      return(NULL)
  }else{
    message("Wrong Input")
    return(NULL)
  }
  
  
  isolate({
    if(is.null(choices))
      return(NULL)
    if(length(choices) >10 ){
      choices <- choices[1:10]
    }
    
    updateSelectInput(session,
                      "kwSearch", 
                      label = "Select from list", 
                      choices = choices, 
                      selected = head(choices,1)
                      )
  })
})

kw_search_Result <- eventReactive(input$submitKW,{
  if(is.null(input$kwSearch)){
    return(NULL)
  }
  keywords <- input$kwSearch
})

observeEvent(input$submitKW,{
  givenData <- input$givenData
  givenData <- paste0(givenData,collapse = ",")
  print(givenData)
  appendData <- kw_search_Result()
  appendData <- paste0(appendData,collapse = ",")
  print(appendData)
  if(givenData != ""){
    givenData <- paste(givenData,appendData,sep = ",")
  } else{
    givenData <- appendData
  }
  updateTextAreaInput(session,"givenData",value = givenData)
})

observeEvent(input$clearGivenData,{
  updateTextAreaInput(session,"givenData",value = "")
})

searchingResult <- eventReactive(input$searchGivenData,{
  if(input$givenDataType == "metabolites"|| input$givenDataType == "genes"){
    if(input$givenMetabolites == "Synonym"){
      synonym <- rampFindSynonymFromSynonym(input$givenData, T)
      source <- rampFindSourceFromId(synonym,T)
      result <- merge(synonym,source)
      result <- unique(result[,2:ncol(result)])
    }else if(input$givenMetabolites == "Biofluid"){
      synonym <- rampFindSynonymFromSynonym(input$givenData,F)
      biof <- rampFastBiofluid(synonym,"analyteSynonym")
      biof
    }else if (input$givenMetabolites == "Catalyzation"){
      cata <- rampFastMulCata(input$givenData)
    }else if (input$givenMetabolites == "Pathway"){
      path <- rampFastPathFromMeta(input$givenData)
    }
  } else if (input$givenDataType == "pathways"){
    meta <- rampFastMetaFromPath(input$givenData)
  } else if (input$givenDataType == "biofluid"){
    meta <- rampFastBiofluid(input$givenData,'ontology')
  }
})

output$searchingResult <- DT::renderDataTable({
  searchingResult()
},
rownames = FALSE)
