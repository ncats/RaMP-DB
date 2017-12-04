dataInput_cata <- eventReactive(input$subText_cata,{
  progress <- shiny::Progress$new()
  
  on.exit(progress$close())
  
  progress$set(message = "Querying databases based on catalyzation ...", value = 0)
  progress$inc(0.3,detail = paste("Send Query ..."))
  
  # rampOut <- rampCataOut(input$KW_cata, 99999)
  rampOut <- rampFastOneCata(input$KW_cata)
  progress$inc(0.7,detail = paste("Done!"))
  return (rampOut[,1:4])
})
summary_cata_out<- eventReactive(input$subText_cata,{
  if (!is.null(nrow(dataInput_cata()))){
    return (paste0("There are(is) ",nrow(dataInput_cata())," relevent items in databases."))
  } else{
    return ("Given metabolites have no search result.")
  }
})

output$summary_cata <- renderText({
  summary_cata_out()
})


observe({
  # x <- rampKWsearch(input$CataInput, "analytesynonym")
  # if (is.null(x)) {
  #   x <- character(0)
  # } else {
  #   x <- as.matrix(x)
  #   x <- matrix(x, ncol = ncol(x), dimnames = NULL)
  # }
  # 
  # updateSelectInput(session, "KW_cata", label = "Select from the list", choices = x, selected = head(x,1))
  choices <- kw_analyte[grepl(input$CataInput,kw_analyte,fixed = T)]
  choices <- choices[order(nchar(choices),choices)]
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
  write.csv(rampOut, file, row.names = FALSE, sep = ",")
})

output$preview_tab4 <- renderUI({
  input$subText_cata
  
  isolate({
    rampout <- dataInput_cata()
    if(length(rampout) == 0){
      return(HTML("<strong> No Summary due to small dataset</strong>"))
    }
    tables <- rampTablize(rampout)
    return(div(HTML(unlist(tables)),class = "shiny-html-output"))
  })
})


# second sub tab

detector_tab4 <- reactiveValues(num = NULL)

observe({
  input$sub_mul_tab4
  
  detector_tab4$num <- 1
})
data_mul_name_tab4 <- eventReactive(input$sub_mul_tab4,{
  rampFastMulCata(input$input_mul_tab4)
})
data_mul_file_tab4 <- eventReactive(input$sub_file_tab4,{
  infile <- input$inp_file_tab4
  if (is.null(infile))
    return(NULL)
  
  rampOut <- rampFileOfAnalytes_tab4(infile)
  
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