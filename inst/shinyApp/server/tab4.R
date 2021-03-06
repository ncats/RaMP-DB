dataInput_cata <- eventReactive(input$subText_cata,{
  tryCatch({
    progress <- shiny::Progress$new()
    
    on.exit(progress$close())
    
    progress$set(message = "Querying databases based on catalyzation ...", value = 0)
    progress$inc(0.3,detail = paste("Send Query ..."))
    
    # rampOut <- rampCataOut(input$KW_cata, 99999)
    rampOut <- RaMP::rampFastCata(input$KW_cata,conpass=.conpass,host = .host,
                                  dbname = .dbname, username = .username,
                                  NameOrIds = input$CataInput_choices)
    progress$inc(0.7,detail = paste("Done!"))
    return (rampOut)
  }, error = function(e) return())
  
})

# summary_cata_out<- eventReactive(input$subText_cata,{
#   if (!is.null(nrow(dataInput_cata()))){
#     return (paste0("There are ",nrow(dataInput_cata())," analytes in the same reaction as input"))
#   } else{
#     return ("Input analyte is not involved in any reactions.")
#   }
# })
# 
# output$summary_cata <- renderText({
#   summary_cata_out()
# })

#Summary
summary_cata_out<- eventReactive(input$subText_cata,{
  if (!is.null(nrow(dataInput_cata()))){
    return (paste0("There are ",nrow(dataInput_cata())," analytes in the same reaction as input"))
  } 
})

summary_cata_out_empty<- eventReactive(input$subText_cata,{
  if (is.null(nrow(dataInput_cata()))){
    return ("Input analyte is not involved in any reactions (or) Analyte name or Source-Id not found.")
  }
})

output$summary_cata <- renderText({
  if (!is.null(summary_cata_out())) {
    return(summary_cata_out())
  } else {
    return(summary_cata_out_empty())
  }
})

#status box
output$statusbox_tab4_subtab1 <- shinydashboard::renderInfoBox({
  if (!is.null(summary_cata_out())) {
    shinydashboard::infoBox(
      "Status",
      HTML(paste("Successful")),
      icon = icon("thumbs-up", lib = "glyphicon"),
      color = "green", fill = TRUE)
  } else if(!is.null(summary_cata_out_empty())) {
    shinydashboard::infoBox(
      "Status",
      HTML(paste("Not-Found")),
      icon = icon("thumbs-down", lib = "glyphicon"),
      color = "yellow", fill = TRUE)
  }
})


observe({
  if(input$CataInput_choices == 'names'){
    if(input$CataInput=="") {
	choices <- ""
    } else {
	#choices <- kw_analyte[grepl(input$CataInput,kw_analyte,fixed = T)]
	choices <- agrep(input$CataInput,kw_analyte,ignore.case=T,value=T)
	choices <- choices[order(nchar(choices),choices)]
    }
  } else if(input$CataInput_choices == 'ids'){
    choices <- kw_source[grepl(input$CataInput,kw_source,fixed = T)]
    choices <- choices[order(nchar(choices),choices)]
  }
  if(is.null(choices))
    return(NULL)
  if(length(choices) >10 ){
    choices <- choices[1:10]
  }
  isolate({
    
    updateSelectInput(session, "KW_cata",
                      label = "Select from the list",
                      choices = choices, selected = head(choices,1)
    )
    })
  })

output$result_cata <- DT::renderDataTable({
  out_src <- dataInput_cata()
},rownames = F)

output$report_cata <- downloadHandler(filename = function() {
  paste0(input$KW_cata, ".csv")
}, content = function(file) {
  rampOut <- dataInput_cata()
  rampOut <- data.frame(rampOut)
  write.csv(rampOut, file, row.names = FALSE)
})


# second sub tab

detector_tab4 <- reactiveValues(num = NULL)

observe({
  input$sub_mul_tab4
  detector_tab4$num <- 1
})

data_mul_name_tab4 <- eventReactive(input$sub_mul_tab4,{
  tryCatch({
    RaMP::rampFastCata(input$input_mul_tab4,conpass=.conpass,
                       dbname = .dbname, username = .username,
                       host = .host,
                       NameOrIds = input$input_mul_tab4_choices)
  }, error = function(e) return())
  
})

#Summary
summary_cata_out_tab2<- eventReactive(input$sub_mul_tab4,{
  if (!is.null(nrow(data_mul_name_tab4()))){
    return ("Result Found")
  } 
})

summary_cata_out_tab2_empty<- eventReactive(input$sub_mul_tab4,{
  if (is.null(nrow(data_mul_name_tab4()))){
    return ("Input analyte is not involved in any reactions (or) Analyte name or Source-Id not found.")
  }
})

output$summary_cata_tab2 <- renderText({
  if (!is.null(summary_cata_out_tab2())) {
    return(summary_cata_out_tab2())
  } else {
    return(summary_cata_out_tab2_empty())
  }
})

#status box
output$statusbox_tab4_subtab2 <- shinydashboard::renderInfoBox({
  if (!is.null(summary_cata_out_tab2())) {
    shinydashboard::infoBox(
      "Status",
      HTML(paste("Successful")),
      icon = icon("thumbs-up", lib = "glyphicon"),
      color = "green", fill = TRUE)
  } else if(!is.null(summary_cata_out_tab2_empty())) {
    shinydashboard::infoBox(
      "Status",
      HTML(paste("Not-Found")),
      icon = icon("thumbs-down", lib = "glyphicon"),
      color = "yellow", fill = TRUE)
  }
})

data_mul_file_tab4 <- eventReactive(input$sub_file_tab4,{
  infile <- input$inp_file_tab4
  print(paste0("Input file is ",infile))
  if (is.null(infile))
    return(NULL)
  
  rampOut <- rampFileOfAnalytes_tab4(infile,conpass=.conpass,host = .host)
  
})

observe({
  input$sub_file_tab4
  detector_tab4$num <- 2
})

# Download table in a csv file.
output$tab4_mul_report <- downloadHandler(filename = function(){
  if (detector_tab4$num == 1){
    paste0("pathwayOutput.csv")
  } else if (detector_tab4$num == 2){
    infile <- input$inp_file_tab4
    paste0(infile[[1,'name']],"Output",".csv")
  }
},
content = function(file) {
  if (detector_tab4$num == 1){
    rampOut <- data_mul_name_tab4()
  } else if (detector_tab4$num == 2){
    rampOut <- data_mul_file_tab4()
  } else{
    message("no data to write ")
  }
  write.csv(rampOut,file,row.names = FALSE)
}
)
output$preview_multi_names_tab4 <- DT::renderDataTable({
  if(is.null(detector_tab4$num))
    return("Waiting for input")
  
  if(detector_tab4$num == 1){
    tb <- data_mul_name_tab4()
  } else if (detector_tab4$num == 2) {
    tb <- data_mul_file_tab4()
  }
  tb
}
)

  output$network <- visNetwork::renderVisNetwork({
	if(!is.null(dataInput_cata())) {
		RaMP::plotCataNetwork(dataInput_cata())
	}

})

   output$networkmulti <- visNetwork::renderVisNetwork({
        if(!is.null(data_mul_name_tab4())) {
                RaMP::plotCataNetwork(data_mul_name_tab4())
        }

})


