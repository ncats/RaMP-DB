

# Tab 4
# 
# 
tabItem4<-  shinydashboard::tabItem(
  tabName = "geneCataComp",
  shinydashboard::tabBox(width = 12,id = "tabset4",height = "100%",
         shiny::tabPanel(strong("Input metabolite or gene one by one."),
                  "",
                  shinydashboard::box(width=6,
                      title = strong("Input synonym of gene or compound"),
                      solidHeader = T,
                      height = "100%",
                      status = "primary",
                      fluidRow(
                        column(12,
                               helpText("If given compound's synonym, it returns all of genes that catalyzed it.
                                        If given gene's synonym, it returns all of compounds catalyzed by this gene"),
                               textInput("CataInput","", placeholder = "Input compound's synonym"),
                               selectInput("KW_cata", "Select from list", choices = NULL),
                               actionButton("subText_cata","Submit")
                               )
                      )
                  ),
                  shinydashboard::box(width = 6,title = strong("Search Result:"),solidHeader = T,
                      status = "primary",
                      fluidRow(
                        div(style = "margin:25px;",
                            downloadButton("report_cata","Generate Report"),
                            hr(),
                            textOutput("summary_cata"),
                            div(style = "height:300px;overflow-x:auto;overflow-y:auto;",
                                helpText("Preview of output only display first 20 items."),
                                DT::dataTableOutput("result_cata")
                            )
                        )
                        
                      )
                  ),
                  hr(),
                  fluidRow(
                    shinydashboard::box(
                      width = 12,
                      solidHeader = T,
                      status = "info",
                      collapsible = T,
                      collapsed = F,
                      title = strong("Summary"),
                      uiOutput("preview_tab4")
                    )
         )),
         shiny::tabPanel(
           title = strong("Input a list of gene or metabolites"),
           fluidRow(
             shinydashboard::box(width = 6,
                 solidHeader = T,
                 status = "primary",
                 title = strong("Input a list of gene or metabolites:"),
                 textAreaInput("input_mul_tab4",label = "",
                               placeholder = "Input list of metabolites or gene in lines or separated by \",\""),
                 actionButton("sub_mul_tab4",label = "Submit")
             ),
             shinydashboard::box(
               width = 6,
               title = strong("Upload the file"),
               solidHeader = T,
               status = "primary",
               fileInput("inp_file_tab4",label = "",
                         multiple = FALSE,
                         accept = c("text/csv","text/comma-separated-values,/text/plain",".csv",".txt"),
                         buttonLabel = "Browse..."),
               actionButton("sub_file_tab4",label = "Upload")
             )
           ),
           hr(),
           fluidRow(
             HTML("<div id='database-group-output'>"),
             shinydashboard::box(width = 12,
                 title = strong("Data preview"),
                 solidHeader = T,
                 status = "primary",
                 downloadButton("tab4_mul_report",label = "Download Table"),
                 hr(),
                 DT::dataTableOutput("preview_multi_names_tab4"),
                 highcharter::highchartOutput("tab4_hc_output")
             ),
             HTML("</div>")
           )
         )
  )
  )




  
