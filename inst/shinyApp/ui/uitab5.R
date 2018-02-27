# Tab 5 Find the ontology relationship between compounds/gene and ontology location
#
#
#
tabItem5<- shinydashboard::tabItem(
  tabName = "geneCompOnto",
  shinydashboard::tabBox(width = 12,
         id = "tabset5",
         height = "100%",
         shiny::tabPanel(strong("Input one metabolite or one ontology"),"",
              shinydashboard::box(width=6,
                  title = strong("Input a metabolite or ontology"),
                  solidHeader = T,
                  height = "100%",status = "primary",
                  fluidRow(
                    column(12,
                           helpText("When metabolites are input, it returns associated ontologies (e.g. biofluid location, cellular location, origins (e.g. drug, food, microbial), and/or tissue location"),
                           helpText("When an ontology name is input, it returns all associated metabolites"),
                           textInput("ontoInput","", placeholder = "Input metabolites or ontology location"),
                           selectInput("KW_onto", "Select from list", choices = NULL),
                           radioButtons("metaOrOnto",
                                        "Define the type of given data",
                                        choices = c(
                                          'Analyte Source ID' = 'ids',
                                          "Analyte Name" = "name",
                                          "Ontology Name" = "ontology")),
                           actionButton("subText_onto","Submit")
                           )
                  )
              ),
              shinydashboard::box(width = 6,title = strong("Preview Result"),solidHeader = T,
                  status = "primary",
                  fluidRow(
                    div(style = "margin:25px;",
                        downloadButton("report_onto","Download Results"),
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
                  title = strong("Input a list of metabolites or ontologies"),
                  textAreaInput("input_mul_tab5",label = "",
                                placeholder = "Input list of metabolites"),
                  radioButtons('input_categories_tab5',label = 'Select from list',
                               choices = c(
                                 'Analyte Source IDs' = 'ids',
                                 "Analyte Names" = "name",
                                 "Ontology Names" = "ontology"
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
                  title = strong("Preview Results"),
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




