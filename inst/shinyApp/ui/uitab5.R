# Tab 5 Find the ontology relationship between compounds/gene and ontology location
# 
# 
# 
tabItem5<- shinydashboard::tabItem(
  tabName = "geneCompOnto",
  shinydashboard::tabBox(width = 12,
         id = "tabset5",
         height = "100%",
         shiny::tabPanel(strong("Input one by one."),"",
              shinydashboard::box(width=6,
                  title = strong("Input metabolites or ontology location"),
                  solidHeader = T,
                  height = "100%",status = "primary",
                  fluidRow(
                    column(12,
                           helpText("If given metabolites' synonym, it returns all of ontology location.
                                    If given ontology, it returns all of metabolites"),
                           textInput("ontoInput","", placeholder = "Input metabolites or ontology location"),
                           selectInput("KW_onto", "Select from list", choices = NULL),
                           radioButtons("metaOrOnto",
                                        "Define the type of given data",
                                        choices = c(
                                          'Analytes Source ID' = 'source',
                                          "Analytes Synonyms" = "name",
                                          "Biofluid Location" = "ontology")),
                           actionButton("subText_onto","Submit")
                           )
                  )
              ),
              shinydashboard::box(width = 6,title = strong("Search Result:"),solidHeader = T,
                  status = "primary",
                  fluidRow(
                    div(style = "margin:25px;",
                        downloadButton("report_onto","Generate Report"),
                        hr(),
                        textOutput("summary_onto"),
                        div(style = "height:300px;overflow-x:auto;overflow-y:auto;",
                            DT::dataTableOutput("result_onto")
                        )
                    )
                  )
              )
         ),
         shiny::tabPanel(
            title = strong("Input a list of analytes or biofluid location"),
            fluidRow(
              shinydashboard::box(width = 8,
                  solidHeader = T,
                  status = "primary",
                  title = strong("Input a list of analyte or biofluid location"),
                  textAreaInput("input_mul_tab5",label = "",
                                placeholder = "Input list of analytes"),
                  radioButtons('input_categories_tab5',label = 'Select from list',
                               choices = c(
                                 'Analytes Source ID' = 'source',
                                 "Analytes Synonyms" = "name",
                                 "Biofluid Location" = "ontology"
                               )),
                  actionButton("sub_mul_tab5",label = "Submit")
                  )
            ),
            hr(),
            fluidRow(
              HTML("<div id='database-group-output'>"),
              shinydashboard::box(width = 12,
                  solidHeader = T,
                  status = "primary",
                  title = strong("Data preview"),
                  downloadButton("tab5_mul_report",
                                 label = "Download"),
                  hr(),
                  DT::dataTableOutput("preview_multi_names_tab5")
              ),
              HTML("</div>")
            )
          )
        )
)




