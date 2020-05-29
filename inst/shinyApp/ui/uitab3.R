tabItem3 <-  shinydashboard::tabItem(
  tabName = "pathFromMeta",
  shinydashboard::tabBox(width = 12, id = "tabset3",
                         shiny::tabPanel(
                           title = strong("Input analyte one by one"),
                           label ="sub_tab_1_TAB3",
                           #Box for Important note
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
                                                 h5("Genes can be searched with the following ID types:  enzymeNomenclature, ensembl, entrez, hmdb, kegg, uniprot. Similar to metabolites, prefix ID types must be added to the ID for searching.")
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
                                                 h4(strong('Input compound name or source id')),
                                                 helpText("Given a gene/metabolite name or ID,retrieve pathways in which the compound is involved in."),
                                                 textInput("compName","",placeholder = "Input compound name or source id"),
                                                 radioButtons("NameOrId","Search by common names or source IDs?",
                                                              choices = c("Names" = "names", "Source ID" = "ids"),selected = "ids"),
                                                 selectInput("KW_synonym", "Select from list", choices = NULL),
                                                 actionButton("submit_compName","Submit"),
                                                 br(),
                                                 br(),
                                                 shinydashboard::infoBoxOutput("statusbox_tab3_subtab1", width = NULL)
                                               ),
                                               shiny::mainPanel(
                                                 width = 6,
                                                 h4(strong('Summary')),
                                                 textOutput("summary_path"),
                                                 br(),
                                                 downloadButton("comp_report","Download Results")
                                               ),
                                               shiny::mainPanel(
                                                 width=12,
                                                 h4(strong('Query Result')),
                                                 fluidRow(
                                                   div(style = "width:100%;height:100%;align:center;",
                                                       helpText("Preview of output only display first 20 items."),
                                                       DT::dataTableOutput("result3")
                                                   )
                                                 )

                                               ))
                         ),
                        # end tabPanel
                         shiny::tabPanel(
                           title = strong("Input multiple analytes (batch query)"),
                           label = "sub_tab_2_TAB3",

                           #Box for Important note

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
                                               title = strong("1.Query Analytes"),
                                               shiny::mainPanel(
                                                 width=6,
                                                 h4(strong('Input a list of metabolites, one per line')),
                                                 textAreaInput("input_mul_tab3",label = "",
                                                               placeholder = "Input list of metabolites, one per line"),
                                                 radioButtons("NameOrSourcemult","Search by names or source IDs?",
                                                              choices = c("Names" = "names", "Source ID" = "ids"),selected = "ids"),
                                                 h4(strong('Input a list of genes, one per line')),
                                                 textAreaInput("input_mul_tab3_genes",label = "",
                                                               placeholder = "Input list of genes, one per line"),
                                                 radioButtons("NameOrSourcemult_genes","Search by names or source IDs?",
                                                              choices = c("Names" = "names", "Source ID" = "ids"),selected = "ids"),
                                                 actionButton("sub_mul_tab3",label = "Submit Query"),
                                                 br(),
                                                 br(),
                                                 shinydashboard::infoBoxOutput("statusbox_tab3_subtab2", width = NULL)
                                               ),
                                               shiny::mainPanel(
                                                 width=6,
                                                 h4(strong("Summary")),
                                                 textOutput("summary_path_tab2"),
                                                 p("The following number of names/ids were mapped:"),
                                                 textOutput("num_mapped_namesids"),
                                                 p("The query returned the following number of pathways per query gene/metabolite:"),
                                                 DT::dataTableOutput("summary_mulpath_out"),
                                                 p("Query results should be visible here"),
                                                 br(),
                                                 br(),
                                                 downloadButton("tab3_mul_report",label = "Download Results")
                                               ),
                                               shiny::mainPanel(
                                                 width = 12,
                                                 HTML("<div id='database-group-output'>"),
                                                 h4(strong("Query Result")),
                                                 helpText("Table output display"),
                                                 div(style = 'width:100%;height:100%;align:center;',
                                                     DT::dataTableOutput("preview_multi_names")),
                                                 HTML("</div>")
                                               )),

                           #Box for Pathway Enrichment Analysis

                           shinydashboard::box(width = 12,
                                               solidHeader = T,
                                               collapsible = TRUE,
                                               collapsed = FALSE,
                                               status = "primary",
                                               title = strong("2.Run Pathway Enrichment Analysis"),
                                               shiny::mainPanel(
                                                 width = 12,
                                                 h4(strong("Set Pathway Enrichment Parameters")),
                                                 tags$head(tags$style(type="text/css", "
             		  	                                  #loadmessage {
                                                                      position: fixed;
                                                                      top: 0px;
                                                                      left: 0px;
                                                                      width: 100%;
                                                                      padding: 5px 0px 5px 0px;
                                                                      text-align: center;
                                                                      font-weight: bold;
                                                                      font-size: 100%;
                                                                      color: #000000;
                                                                      background-color: #FFA500;
                                                                      z-index: 105;
                                                                      }
                                                                      ")),
                                                 conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                                  tags$div("Processing request...",id="loadmessage")),
                                                 actionButton("runFisher","Run Pathway Enrichment (please be patient!)"),
                                                 HTML("<br>"),
                                                 HTML("<br>"),
                                                 textOutput("FisherTestResultWithoutFilter_AnalyteType"),
                                                 HTML("<br>"),
                                                 DT::dataTableOutput("FisherTestResultWithoutFilter_fishresults",width = '100%'),
                                                 textOutput("fishersProgress")
                                               ),
                                               shiny::mainPanel(
                                                 width = 12
                                                 # h4(strong("Results of Pathway Enrichment Analysis")),
                                                 # column(width = 6,
                                                 #        selectInput("show_cluster","Display pathways in cluster:",choices = 1),
                                                 #        textOutput("cluster_summary_text")
                                                 # ), # end columnb
                                                 # column(width = 6,
                                                 #        plotOutput("cluster_summary_plot",height = "300px")
                                                 # ), # end column
                                                 # DT::dataTableOutput("results_fisher",width = '100%')
                                               )),

                           #Box for Filter and Cluster pathways

                           shinydashboard::box(width = 12,
                                               solidHeader = T,
                                               collapsible = TRUE,
                                               collapsed = FALSE,
                                               status = "primary",
                                               title = strong("3.Filter and cluster Pathways"),
                                               shiny::mainPanel(
                                                 width = 12,
                                                 h4(strong("Set Filtering/Clustering Parameters")),
                                                 numericInput("p_holmadj_cutoff",
                                                              "Select cutoff for Holm Adjusted p-values (will only return pathways with Holm adj p-values < cutoff)",
                                                              value = 0.01,
                                                              min=0,max=1,
                                                              width = "80%"),
                                                 h4("Clustering Parameters:"),
                                                 numericInput("perc_analyte_overlap", "Overlap threshold (0-1) for pathways to be considered similar (for medoid establishment):",min = 0.01, max = 1, value = 0.2, step = 0.01),
                                                 numericInput("min_pathway_tocluster", "Number of similar neighbors required (for medoid establishment):",min = 1, max = 100, value = 2, step = 1),
                                                 numericInput("perc_pathway_overlap", "Overlap threshold (0-1) for pathways to be clustered:",min = 0.01, max = 1, value = 0.75, step = 0.01),
                                                 actionButton("runClustering","Filter and Cluster Results"),
                                                 h4(strong("Results of Pathway Enrichment Analysis")),
                                                 column(width = 12,
                                                        selectInput("show_cluster","Display pathways in cluster:",choices = 1),
                                                        textOutput("cluster_summary_text")
                                                 ), # end columnb
                                                 br(),
                                                 br(),
                                                 column(width = 12,
                                                 shinydashboard::infoBoxOutput("results_status", width = NULL)
                                                 ),
                                                 # end column
                                                 textOutput("filterStatus"),
                                                 DT::dataTableOutput("results_fisher",width = '100%'),
                                                 column(width = 6,
                                                        plotOutput("cluster_summary_plot",height = "300px")
                                                 )
                                               )),

                            #Box for Download Results.

                           shinydashboard::box(width = 12,
                                               solidHeader = T,
                                               collapsible = TRUE,
                                               collapsed = FALSE,
                                               status = "primary",
                                               title = strong("4.Download Results"),
                                               shiny::mainPanel(
                                                 p("Significant pathways are returned under 'Results of Pathway Enrichment Analysis' and can be downloaded by clicking 'Download Results'"),
                                                 p("Note that only pathways that contain at least 3 analytes from the user input will be output"),
                                                 #textOutput("summary_Fisher"),
                                                 p("The following number of pathways per database were processed:"),
                                                 DT::dataTableOutput("summary_fisher",width = '80%'),
                                                 downloadButton("fisher_stats_report",label = "Download Results")
                                               ))


                         )#end tabpanel
  ) # end Tab Box
) # end Tabitem
