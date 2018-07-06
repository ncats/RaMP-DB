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
                           h5(strong("IMPORTANT NOTE about inputting source IDs:")),
                           h5("When inputting source IDs, it is important to add a prefix to denote the id type.  This is important because it is possible for two different metabolites to have the same IDs, although each ID may be from a different database source."),
                           h5("Metabolites can be searched with the following ID types: CAS, chebi, chemspider, hmdb, kegg, LIPIDMAPS, and pubchem.  To search for a metabolite, the ID type must be added as a prefix.  For example, the compound 'HMDB0000562' must be searched by 'hmdb:HMDB0000562', the compound '16737' must be searched by 'chebi:16737'."),
                           textInput("ontoInput","", placeholder = "Input metabolites or ontology location"),
                           selectInput("KW_onto", "Select from list", choices = NULL),
                           radioButtons("metaOrOnto",
                                        "Define the type of given data",
                                        choices = c(
                                          'Analyte Source ID' = 'ids',
                                          "Analyte Name" = "name",
                                          "Ontology Name" = "ontology")),
                           actionButton("subText_onto","Submit"),
                           br(),
                           br(),
                           shinydashboard::infoBoxOutput("statusbox_tab5_subtab1", width = NULL)
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
                  actionButton("sub_mul_tab5",label = "Submit"),
                  br(),
                  br(),
                  textOutput("summary_onto_tab2"),
                  br(),
                  shinydashboard::infoBoxOutput("statusbox_tab5_subtab2", width = NULL)
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




