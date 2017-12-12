server <- function(input, output, session) {
  
  session$onSessionEnded(function(){
    RaMP:::killDbConnections()
    shiny::stopApp()
  }) # close shiny app when close browser.
  
  # Tab 1 for convert synonym to synonyms that have connection to it.
  #source(file.path("server","tab1.R"), local = TRUE)$value
  
  # Tab 2 given a pathway name, output all synonym
  source(file.path("server","tab2.R"),local = TRUE)$value
  
  # Tab 3 given the metabolites, return all of pathway's name that get involved in that compound
  source(file.path("server","tab3.R"),local = TRUE)$value
 
  # Tab 4 return metabolites based on catalyzation server functions are defined below.
  source(file.path("server","tab4.R"),local = TRUE)$value
  
  # Tab 5 Find the ontology relationship between compounds/gene and ontology location
  # source(file.path("server","tab5.R"),local = TRUE)$value
  
  # Tab 6 Customize query to database 
  # source(file.path("server","tab6.R"),local = TRUE)$value
  
}
