tabItem2<-  shinydashboard::tabItem(
  tabName = "metaFromPath",
  shinydashboard::tabBox(width = 12,id ="tabset2",
                         
                         shiny::tabPanel(
                           title = strong("Input pathway one by one"),
                           label = "sub_tab_1_TAB2",
                           #box for Query Analytes,Summary and Results
                           shinydashboard::box(width = 12,
                                               solidHeader = T,
                                               collapsible = TRUE,
                                               collapsed = FALSE,
                                               status = "primary",
                                               title = strong("Query Analytes"),
                                               shiny::mainPanel(
                                                 width = 6,
                                                 h4(strong('Input pathway you want to search for')),
                                                 helpText("Given a pathway name, it returns all compounds contained in that pathway"),
                                                 textInput("singleInput2","",placeholder ="Input one pathway name"),
                                                 selectInput("KW_path", "Select from list", choices = NULL),
                                                 actionButton("subText2","Submit"),
                                                 br(),
                                                 br(),
                                                 shinydashboard::infoBoxOutput("statusbox_tab2_subtab1", width = NULL),
                                                 hr()
                                               ),
                                               shiny::mainPanel(
                                                 width = 6,
                                                 h4(strong('Summary')),
                                                 br(),
                                                 textOutput("summary_search"),
                                                 br(),
                                                 downloadButton("result_file","Generate Report")
                                               ),
                                               shiny::mainPanel(
                                                 width = 12,
                                                 h4(strong("Select type of analytes to be returned")),
                                                 radioButtons("geneOrComp2","Return compound or gene", choices = c(
                                                   "Compound" = "compound", "Gene" = "gene","Both" = "both"
                                                 ),selected = "compound"),
                                                 br()
                                               ),
                                               shiny::mainPanel(
                                                 width = 12,
                                                 h4(strong('Query Result')),
                                                 fluidRow(
                                                   HTML("<div id='database-group-output'>"),
                                                   DT::dataTableOutput("result",width = "100%",height = "100%"),
                                                   HTML("</div>")
                                                 )
                                               ))
                         ),#end tabpanel
                         shiny::tabPanel(
                           title = strong("Input a list of pathways"),
                           "",
                           #box for Query Analytes,Summary and Results
                           shinydashboard::box(width = 12,
                                               solidHeader = T,
                                               collapsible = TRUE,
                                               collapsed = FALSE,
                                               status = "primary",
                                               title = strong("Query Analytes"),
                                               shiny::mainPanel(
                                                 width = 6,
                                                 h4(strong('Input a list of pathways, one per line')),
                                                 textAreaInput("input_mul_tab2",label = "",placeholder = "Input list of pathways, one per line"),
                                                 actionButton("sub_mul_tab2",label = "Submit"),
                                                 br(),
                                                 br(),
                                                 shinydashboard::infoBoxOutput("statusbox_tab2_subtab2", width = NULL)
                                                 #                    shinydashboard::box(
                                                 #                      width = 6,
                                                 #                      title = strong("Upload the file"),
                                                 #                      fileInput("inp_file_tab2",label = "",
                                                 #                                multiple = FALSE,
                                                 #                                accept = c("text/csv","text/comma-separated-values,/text/plain",".csv",".txt"),
                                                 #                                buttonLabel = "Browse..."),
                                                 #                      actionButton("sub_file_tab2",label = "Submit")
                                                 #                    )
                                               ),
                                               shiny::mainPanel(
                                                 width = 6,
                                                 h4(strong('Summary')),
                                                 br(),
                                                 textOutput("summary_search_tab2"),
                                                 br(),
                                                 downloadButton("tab2_mul_report",label = "Download Results")
                                               ),
                                               shiny::mainPanel(
                                                 width = 12,
                                                 h4(strong('Query Result')),
                                                 fluidRow(
                                                   DT::dataTableOutput("preview_multi_names_tab2",
                                                                       width = "100%",
                                                                       height = "100%"),
                                                   highcharter::highchartOutput("tab2_hc_output"),
                                                   textOutput("hc_click_output_tab2"),
                                                   tableOutput("stats_fisher_tab2")
                                                 )
                                               ))
                         )#end tabpanel
                         
  )
)



