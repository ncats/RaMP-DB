dataInput_path <- eventReactive(input$subText2,{
  progress <- shiny::Progress$new()
  
  on.exit(progress$close())
  
  progress$set(message = "Querying databases to find metabolites ...", value = 0)
  progress$inc(0.3,detail = paste("Send Query ..."))
  
  rampOut <- rampFastMetaFromPath(input$KW_path)
  rampOut <- rampOut[,1:(ncol(rampOut) - 1)]
  if(input$geneOrComp2 != "both"){
    rampOut <- rampOut[rampOut[,2] == input$geneOrComp2,]
  }
  progress$inc(0.7,detail = paste("Down!"))
  
  return(rampOut)
})

path_text_out<- eventReactive(input$subText2,{
  if(!is.null(dataInput_path())){
    return (paste0("There are(is) ",nrow(dataInput_path())," relevent items in the databases."))
  } else{
    return ("Given pathway not found.")
  }
  
})

output$summary_synonym <- renderText({
  path_text_out()
})

output$result2 <- DT::renderDataTable({
  out_src <- dataInput_path()
  if (is.null(out_src))
    return(NULL)
  
  
  return(out_src)
},rownames = FALSE)


observe({
  
  choices <- kw_pathway[grepl(input$singleInput2,kw_pathway,fixed = T)]
  choices <- choices[order(nchar(choices),choices)]
  if(is.null(choices))
    return(NULL)
  if(length(choices) >10 ){
    choices <- choices[1:10]
  }
  isolate({
    updateSelectInput(session, "KW_path",
                      label = "Select from list",
                      choices = choices, selected = head(choices,1))
    })
})


output$path_report <- downloadHandler(filename = function() {
  paste0(input$KW_path, ".csv")
}, content = function(file) {
  ramp <- dataInput_path()
  write.csv(ramp, file, append = TRUE, row.names = FALSE)
})


output$preview_tab2 <- renderUI({
  input$subText2
  isolate({
    tables <- rampTablize(dataInput_path())
    return(div(HTML(unlist(tables)),class = "shiny-html-output"))
  })
})


# Second Panel 

detector_tab2 <- reactiveValues(num = NULL)

observe({
  input$sub_mul_tab2
  
  isolate({
    detector_tab2$num <- 1
  })
})
data_mul_name_tab2 <- eventReactive(input$sub_mul_tab2,{
  if(is.null(input$sub_mul_tab2))
    return(NULL)
  rampFastMetaFromPath(input$input_mul_tab2)
})
data_mul_file_tab2 <- eventReactive(input$sub_file_tab2,{
  infile <- input$inp_file_tab2
  if (is.null(infile))
    return(NULL)
  
  rampFileOfPathways_tab2(infile)
})


observe({
  input$sub_file_tab2
  
  detector_tab2$num <- 2
})

# Download table in a csv file.
output$tab2_mul_report <- downloadHandler(filename = function(){
  if (detector_tab2$num == 1){
    paste0("pathwayOutput.csv")
  } else if (detector_tab2$num == 2){
    infile <- input$inp_file_tab2
    paste0(infile[[1,'name']],"Output",".csv")
  }
},
content = function(file) {
  if (detector_tab2$num == 1){
    rampOut <- data_mul_name_tab2()
  } else if (detector_tab2$num == 2){
    rampOut <- data_mul_file_tab2()
  }
  write.csv(rampOut,file,row.names = FALSE)
}
)

tb_data_tab2 <- reactive({
  if(detector_tab2$num == 1){
    tb <- data_mul_name_tab2()
  } else if (detector_tab2$num == 2) {
    tb <- data_mul_file_tab2()
  }
  tb
})
output$preview_multi_names_tab2 <- DT::renderDataTable({
  if(is.null(detector_tab2$num))
    return("Waiting for input")
  
  tb_data_tab2()
},
rownames = FALSE
)
output$tryingTab2 <- renderText({
  input$input_mul_tab2
})
