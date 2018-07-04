# Function that runs getAnalyteFromPathway to return analytes from input pathway(s)
dataInput_path <- eventReactive(input$subText2,{
  tryCatch({
    progress <- shiny::Progress$new()
    
    on.exit(progress$close())
    
    progress$set(message = "Querying databases to find metabolites ...", value = 0)
    progress$inc(0.3,detail = paste("Send Query ..."))
    
    rampOut <- RaMP::getAnalyteFromPathway(input$KW_path,
                                           conpass=.conpass,host = .host, 	
                                           dbname = .dbname, username = .username)
    if(input$geneOrComp2 != "both"){
      rampOut <- rampOut[rampOut$geneOrCompound == input$geneOrComp2,]
    }
    progress$inc(0.7,detail = paste("Down!"))
    
    print(dim(rampOut))
    return(rampOut)
  }, error = function(e) return())
})

# Counts the number of lines in resulting query
# path_text_out<- eventReactive(input$subText2,{
#   if(!is.null(dataInput_path())){
#     return (paste0("There are(is) ",nrow(dataInput_path())," ",
# 	input$geneOrComp2," returned."))
#   } else{
#     return ("Given pathway not found.")
#   }
# 
# })
# 
# output$summary_search <- renderText({
#   path_text_out()
# })

# Counts the number of lines in resulting query
path_text_out<- eventReactive(input$subText2,{
  if(!is.null(dataInput_path())){
    return (paste0("There are(is) ",nrow(dataInput_path())," ",
                   input$geneOrComp2," returned."))
  } 
})

path_text_out_empty<- eventReactive(input$subText2,{
  if(is.null(dataInput_path())){
    return ("Given pathway not found.")
  } 
})

output$summary_search <- renderText({
  if(!is.null(path_text_out())) {
    return(path_text_out())
  } else {
    path_text_out_empty()
  }
})

#Status Box
output$statusbox_tab2_subtab1 <- shinydashboard::renderInfoBox({
  if (!is.null(path_text_out())) {
    shinydashboard::infoBox(
      "Status",
      HTML(paste("Successful")),
      icon = icon("thumbs-up", lib = "glyphicon"),
      color = "green", fill = TRUE)
  } else if (!is.null(path_text_out_empty())) {
    shinydashboard::infoBox(
      "Status",
      HTML(paste("Not-Found")),
      icon = icon("thumbs-down", lib = "glyphicon"),
      color = "yellow", fill = TRUE)
  } 
})

# outputs results of query to a table
output$result <- DT::renderDataTable({
  out_src <- dataInput_path()
  if (is.null(out_src))
    return(NULL)
  return(out_src[colnames(out_src) != 'pathwayCategory'])
},rownames = FALSE)


observe({
  #choices <- kw_pathway[grepl(input$singleInput2,kw_pathway,ignore.case=TRUE)]
    if(input$singleInput2=="") {
        choices <- ""
    } else {
	choices <- agrep(input$singleInput2,kw_pathway,ignore.case=TRUE,value=TRUE)
	choices <- choices[order(nchar(choices),choices)]
    }
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

# Outputs query results to a file
output$result_file <- downloadHandler(filename = function() {
  paste0(input$KW_path, ".csv")
}, content = function(file) {
  ramp <- dataInput_path()
  write.csv(ramp, file, append = TRUE, row.names = FALSE)
})

############################################
# Second Panel
###########################################

detector_tab2 <- reactiveValues(num = NULL)

observe({
  input$sub_mul_tab2

  isolate({
    detector_tab2$num <- 1
  })
})

# Runs query given a list of pathways
data_mul_name_tab2 <- eventReactive(input$sub_mul_tab2,{
  tryCatch({
    if(is.null(input$sub_mul_tab2))
      return(NULL)
    RaMP::getAnalyteFromPathway(input$input_mul_tab2,
                                conpass=.conpass,host = .host,dbname=.dbname,username=.username)
  }, error = function(e) return())
})
# Download table in a csv file.
output$tab2_mul_report <- downloadHandler(filename = function(){
  if (detector_tab2$num == 1){
    paste0("pathwayOutput.csv")
  }
},
content = function(file) {
  if (detector_tab2$num == 1){
    rampOut <- data_mul_name_tab2()
  } 
  write.csv(rampOut,file,row.names = FALSE)
}
)

tb_data_tab2 <- reactive({
  if(detector_tab2$num == 1){
    tb <- data_mul_name_tab2()
  } 
  tb[colnames(tb) != 'pathwayCategory']
})

#Summary
path_text_out_tab2<- eventReactive(input$sub_mul_tab2,{
  if(!is.null(tb_data_tab2())){
    return("Result Found")
  } 
})

path_text_out_tab2_empty<- eventReactive(input$sub_mul_tab2,{
  if(is.null(tb_data_tab2())){
    return("Given pathway not found.")
  } 
})


output$summary_search_tab2 <- renderText({
  if(!is.null(path_text_out_tab2())) {
    return(path_text_out_tab2())
  } else {
    path_text_out_tab2_empty()
  }
})

#Status box
output$statusbox_tab2_subtab2 <- shinydashboard::renderInfoBox({
  if (!is.null(path_text_out_tab2())) {
    shinydashboard::infoBox(
      "Status",
      HTML(paste("Successful")),
      icon = icon("thumbs-up", lib = "glyphicon"),
      color = "green", fill = TRUE)
  } else if(!is.null(path_text_out_tab2_empty())) {
    shinydashboard::infoBox(
      "Status",
      HTML(paste("Not-Found")),
      icon = icon("thumbs-down", lib = "glyphicon"),
      color = "yellow", fill = TRUE)
  } 
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
