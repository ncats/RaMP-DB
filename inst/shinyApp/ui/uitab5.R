# Tab 5 Find the ontology relationship between compounds/gene and ontology location
#
#
#
tabItem5<- shinydashboard::tabItem(
  tabName = "geneCompOnto",
  shinydashboard::tabBox(width = 12,
         id = "tabset5",
         height = "100%",
         shiny::tabPanel(
           title = strong("Input one metabolite or one ontology"),
           "",
           #box for important note
           shinydashboard::box(width = 12,
                               solidHeader = T,
                               collapsible = TRUE,
                               collapsed = TRUE,
                               status = "primary",
                               title = strong("IMPORTANT NOTE/INSTRUCTIONS"),
                               shiny::mainPanel(
                                 width = 12,
                                 h5("When inputting source IDs, it is important to add a prefix to denote the id type.  This is important because it is possible for two different metabolites to have the same IDs, although each ID may be from a different database source."),
                                 h5("Metabolites can be searched with the following ID types: CAS, chebi, chemspider, hmdb, kegg, LIPIDMAPS, and pubchem.  To search for a metabolite, the ID type must be added as a prefix.  For example, the compound 'HMDB0000562' must be searched by 'hmdb:HMDB0000562', the compound '16737' must be searched by 'chebi:16737'.")
                               )),
           #box for Query Analytes,Summary and Results
           shinydashboard::box(width = 12,
                               solidHeader = T,
                               collapsible = TRUE,
                               collapsed = FALSE,
                               status = "primary",
                               title = strong("Query Analytes"),
                               shiny::mainPanel(
                                 width = 6,
                                 h4(strong('Input a metabolite or ontology')),
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
                                 actionButton("subText_onto","Submit"),
                                 br(),
                                 br(),
                                 shinydashboard::infoBoxOutput("statusbox_tab5_subtab1", width = NULL)
                               ),
                               shiny::mainPanel(
                                 width = 6,
                                 h4(strong('Summary')),
                                 textOutput("summary_onto"),
                                 br(),
                                 downloadButton("report_onto","Download Results")
                               ),
                               shiny::mainPanel(
                                 width = 12,
                                 h4(strong('Query Result')),
                                 helpText("Table output display"),
                                 fluidRow(
                                  DT::dataTableOutput("result_onto")
                                 )
                               ))
         ),#end tabpanel

         shiny::tabPanel(
           title = strong("Input a list of analytes or biofluid location"),
           "",
           #box for important note
           shinydashboard::box(width = 12,
                               solidHeader = T,
                               collapsible = TRUE,
                               collapsed = TRUE,
                               status = "primary",
                               title = strong("IMPORTANT NOTE/INSTRUCTIONS"),
                               shiny::mainPanel(
                                 width = 12,
                                 h5("When inputting source IDs, it is important to add a prefix to denote the id type.  This is important because it is possible for two different metabolites to have the same IDs, although each ID may be from a different database source.")
                               )),
           #box for Query Analytes,Summary and Results
           shinydashboard::box(width = 12,
                               solidHeader = T,
                               collapsible = TRUE,
                               collapsed = FALSE,
                               status = "primary",
                               title = strong("Query Analytes"),
                               shiny::mainPanel(
                                 width = 6,
                                 h4(strong('Input a list of metabolites or ontologies')),
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
                                 shinydashboard::infoBoxOutput("statusbox_tab5_subtab2", width = NULL)
                               ),
                               shiny::mainPanel(
                                 width = 6,
                                 h4(strong('Summary')),
                                 textOutput("summary_onto_tab2"),
                                 br(),
                                 downloadButton("tab5_mul_report",label = "Download Results")
                               ),
                               shiny::mainPanel(
                                 width = 12,
                                 h4(strong('Query Result')),
                                 helpText("Table output display"),
                                 fluidRow(
                                   HTML("<div id='database-group-output'>"),
                                   DT::dataTableOutput("preview_multi_names_tab5"),
                                   HTML("</div>")
                                 )
                               ))
         )#end tabpanel
         
        )
)




