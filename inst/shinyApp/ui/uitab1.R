tabItem1<-shinydashboard::tabItem(
  tabName = "genesFromCompounds",
  shinydashboard::tabBox(width = 12,id ="tabset1",
                         shinydashboard::tabPanel(title = strong("Input analyte one by one"),
                  label = "sub_tab_1_TAB1",
                  fluidRow(
                    box(width=6,
                        title = strong("Input synonym that you want to search for"),
                        solidHeader = T,
                        height = "100%",
                        status = "primary",
                        fluidRow(
                          column(12,
                                 helpText("Given synonym, it returns all of genes or compounds 
                                          which have same pathway"),
                                 textInput("singleInput","",placeholder = "Input compound's synonym"),
                                 selectInput("kw_search", "Select from list", choices = NULL),
                                 actionButton("subText","Submit")
                                 )
                        )
                    ),
                    box(width = 6,
                        title = strong("Search Result:"),
                        solidHeader = T,
                        height = "100%",
                        status = "primary",
                        fluidRow(
                          div(style = "margin:25px;",
                              downloadButton("report","Generate Report"),
                              hr(),
                              textOutput("numOfitems"),
                              div(style = "height:300px;overflow-x:auto;overflow-y:auto;",
                                  helpText("Preview of output only display first 20 items."),
                                  dataTableOutput("result")
                              )
                          )
                          
                        )
                    )
                  ),
                  fluidRow(
                    HTML("<div id ='databases-group-output'>"),
                    box(width = 12,
                        solidHeader = T,
                        status = "info",
                        collapsible = T,
                        title = strong("Summary"),
                        uiOutput("preview_tab1")
                    ),
                    HTML("</div>")
                  )
            ),
            shinydashboard::tabPanel(title = "Input a list of analytes",
                  fluidRow(
                    box(width = 6,
                        status = "primary",
                        solidHeader = T,
                        title = strong("Input a list of metabolites"),
                        textAreaInput("input_mul_tab1",label = "",
                                      placeholder = "Input list of metabolites in lines or separated by \",\""),
                        actionButton("sub_mul_tab1",label = "Submit")
                    ),
                    box(
                      width = 6,
                      solidHeader = T,
                      status = "primary",
                      title = strong("Upload the file"),
                      fileInput("inp_file_tab1",label = "",
                                multiple = TRUE,
                                accept = c("text/csv","text/comma-separated-values,/text/plain",".csv",".txt"),
                                buttonLabel = "Browse..."),
                      actionButton("sub_file_tab1",label = "Upload")
                    )
                  ),
                  hr(),
                  fluidRow(
                    HTML("<div id='database-group-output'>"),
                    box(width = 12,
                        solidHeader = T,
                        status = "primary",
                        title = strong("Data preview"),
                        downloadButton("tab1_mul_report",label = "Download Table"),
                        hr(),
                        dataTableOutput("preview_multi_names_tab1"),
                        highchartOutput("tab1_hc_output"),
                        conditionalPanel(
                          condition = "output.tab1_hc_output != null",
                          selectizeInput(inputId = "pathway_from_hc","",
                                         multiple = T,
                                         choice = NULL),
                          textOutput("trying"),
                          div(style="display:inline-block;",
                              actionButton("meta_from_path_fire",
                                           "Search metabolites in above pathways")
                          ),
                          div(style="display:inline-block;margin:10px;",
                              downloadButton("meta_from_path_download",
                                             label = "Download Reprot")),
                          hr(),
                          dataTableOutput("meta_from_path_preview")
                        )
                    ),
                    HTML("</div>")
                  )
         )
    )
)


