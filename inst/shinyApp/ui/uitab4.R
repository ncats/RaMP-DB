# Tab 4
#
#
tabItem4<-  shinydashboard::tabItem(
  tabName = "geneCataComp",
  shinydashboard::tabBox(width = 12,id = "tabset4",height = "100%",
                         shiny::tabPanel(
                           title = strong("Input metabolite or gene one by one"),
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
                                                 h4(strong('Input metabolite or gene')),
                                                 helpText("If a compound is input, all genes that catalyze reactions involving the compound are returned."),
                                                 helpText("Conversely, if a gene is input, all compounds in reactions that are catalyzed by that gene are returned"),
                                                 textInput("CataInput","", placeholder = "Input compound name or id"),
                                                 radioButtons('CataInput_choices',label = 'Search by name or source id',
                                                              choices = c('Names' = 'names','Source ID' = 'ids')),
                                                 selectInput("KW_cata", "Select from list", choices = NULL),
                                                 actionButton("subText_cata","Submit"),
                                                 br(),
                                                 br(),
                                                 shinydashboard::infoBoxOutput("statusbox_tab4_subtab1", width = NULL)
                                               ),
                                               shiny::mainPanel(
                                                 width = 6,
                                                 h4(strong('Summary')),
                                                 textOutput("summary_cata"),
                                                 br(),
                                                 downloadButton("report_cata","Download Results")
                                               ),
                                               shiny::mainPanel(
                                                 width = 12,
                                                 h4(strong('Query Result')),
                                                 helpText("Table output display"),
                                                 fluidRow(
                                                  DT::dataTableOutput("result_cata")
                                                 )
                                               )),
                           #box for interaction network display
                           shinydashboard::box(width = 12,
                                               height= "1000px",
                                               solidHeader = T,
                                               status = "primary",
                                               collapsible = T,
                                               collapsed = F,
                                               title = strong("Visuazlize gene-metabolite interaction network"),
                                               visNetwork::visNetworkOutput("network",height = '1000px')) # end box
                         ),#end tabpanel
                         shiny::tabPanel(
                           title = strong("Input a list of genes or metabolites"),
                           "",
                           #box for Important note
                           shinydashboard::box(width = 12,
                                               solidHeader = T,
                                               collapsible = TRUE,
                                               collapsed = TRUE,
                                               status = "primary",
                                               title = strong("IMPORTANT NOTE/INSTRUCTIONS"),
                                               shiny::mainPanel(
                                                 width = 12,
                                                 h5("When inputting source IDs, it is important to add a prefix to denote the id type.  This is important because it is possible for two different metabolites to have the same IDs, although each ID may be from a different database source."),
                                                 h5("Metabolites can be searched with the following ID types: CAS, chebi, chemspider, hmdb, kegg, LIPIDMAPS, and pubchem.  To search for a metabolite, the ID type must be added as a prefix.  For example, the compound 'HMDB0000562' must be searched by 'hmdb:HMDB0000562', the compound '16737' must be searched by 'chebi:16737'."),
                                                 h5("Genes can be searched with the following ID types:  enzymeNomenclature, ensembl, entrez, hmdb, kegg, uniprot. Similar to metabolites, prefix ID types must be added to the ID for searching."),
                                                 h5(strong("Only one type of ID can be input: either source ID or name, not both"))
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
                                                 h4(strong('Input a list of genes or metabolites')),
                                                 textAreaInput("input_mul_tab4",label = "",
                                                               placeholder = "Input list of genes or metabolites, one per line"),
                                                 radioButtons('input_mul_tab4_choices',label = 'Search by source id or name?',
                                                              choices = c('Name' = 'names','Source ID' = 'ids')),
                                                 actionButton("sub_mul_tab4",label = "Submit"),
                                                 br(),
                                                 br(),
                                                 shinydashboard::infoBoxOutput("statusbox_tab4_subtab2", width = NULL)
                                               ),
                                               shiny::mainPanel(
                                                 width = 6,
                                                 h4(strong('Summary')),
                                                 textOutput("summary_cata_tab2"),
                                                 br(),
                                                 downloadButton("tab4_mul_report",label = "Download Results")
                                               ),
                                               shiny::mainPanel(
                                                 width = 12,
                                                 h4(strong('Query Result')),
                                                 helpText("Table output display"),
                                                 fluidRow(
                                                   DT::dataTableOutput("preview_multi_names_tab4")
                                                 )
                                               )),
                           #box for interaction network display
                           shinydashboard::box(width = 12,
                                               height="1000px",
                                               solidHeader = T,
                                               status = "primary",
                                               collapsible = T,
                                               collapsed = F,
                                               title = strong("Visuazlize gene-metabolite interaction network"),
                                               visNetwork::visNetworkOutput("networkmulti",height = '1000px')
                           )
                         )#end tabpanel
                         
     )# end tabBox
  )# end tabItem





