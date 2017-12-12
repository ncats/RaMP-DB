tabItem2<-  shinydashboard::tabItem(
  tabName = "metaFromPath",
  shinydashboard::tabBox(width = 12,id ="tabset2",
         shiny::tabPanel(title = strong("Input pathway one by one"),
                  label = "sub_tab_1_TAB2",
                  fluidRow(
                    shinydashboard::box(
                      width = 6,
                      title = strong("Input pathway you want to search for"),
                      solidHeader = T,
                      status = "primary",
                      fluidRow(
                        column(12,
                               helpText("Given pathway's name, it returns 
                                        all of compounds' synonym contained in that pathway"),
                               textInput("singleInput2","",placeholder ="Input one pathway name"),
                               selectInput("KW_path", "Select from list", choices = NULL),
                               actionButton("subText2","Submit"),
                               hr())
                        ),
                      fluidRow(
                        shinydashboard::box(width = 12,
                            title = "Select type of metabolites returned by databases",
                            radioButtons("geneOrComp2","Return compound or gene", choices = c(
                              "Compound" = "compound", "Gene" = "gene","Both" = "both"
                            ),selected = "compound"),
                            status = "primary",
                            solidHeader = T,
                            collapsible = TRUE
                        )
                      )
                    ),
                    shinydashboard::box(width = 6,
                        title = strong("Search Result:"),
                        solidHeader = T,
                        status = "primary",
                        fluidRow(
                          div(
                            style = "margin:25px",
                            downloadButton("path_report","Generate Report"),
                            hr(),
                            textOutput("summary_synonym"),
                            div(style = "height:300px;overflow-y:auto;overflow-x:auto;",
                                helpText("Preview of output only display first 20 items."),
                                DT::dataTableOutput("result2",width = "80%",height = "100%")
                            )
                          )
                        )
                    )
                  ),
                  fluidRow(
                    HTML("<div id='database-group-output'>"),
                    shinydashboard::box(
                      width = 12,
                      status = "info",
                      solidHeader = T,
                      collapsible = T,
                      title = strong("Summary"),
                      uiOutput("preview_tab2")
                    ),
                    HTML("</div>")
                  )
         ),
         shiny::tabPanel(title = "Input a list of pathways",
                  fluidRow(
                    shinydashboard::box(width = 6,
                        title = strong("Input a list of pathways:"),
                        textAreaInput("input_mul_tab2",label = "",
                                      placeholder = "Input list of metabolites in lines or separated by \",\""),
                        actionButton("sub_mul_tab2",label = "Submit")
                    ),
                    shinydashboard::box(
                      width = 6,
                      title = strong("Upload the file"),
                      fileInput("inp_file_tab2",label = "",
                                multiple = FALSE,
                                accept = c("text/csv","text/comma-separated-values,/text/plain",".csv",".txt"),
                                buttonLabel = "Browse..."),
                      actionButton("sub_file_tab2",label = "Upload")
                    )
                  ),
                  hr(),
                  fluidRow(
                    HTML("<div id='database-group-output'>"),
                    shinydashboard::box(width = 12,
                        title = strong("Data preview"),
                        downloadButton("tab2_mul_report",label = "Download Table"),
                        hr(),
                        DT::dataTableOutput("preview_multi_names_tab2",
                                            width = "100%",
                                            height = "100%"),
                        highcharter::highchartOutput("tab2_hc_output"),
                        textOutput("hc_click_output_tab2"),
                        tableOutput("stats_fisher_tab2")
                    ),
                    HTML("</div>")
                  )
         )
  )
)
