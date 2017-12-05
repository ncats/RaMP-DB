tabItem3 <-  shinydashboard::tabItem(
  tabName = "pathFromMeta",
  shinydashboard::tabBox(width = 12, id = "tabset3",
         tabPanel(
           title = strong("Input analyte one by one"),
           label ="sub_tab_1_TAB3",
           box(
             width = 6,
             title = strong("Input synonym"),
             solidHeader = T,
             status = "primary",
             fluidRow(
               column(12,
                      helpText("Given compound's name, 
                               it returns pathways in which the compound get involved in."),
                      textInput("compName","",placeholder = "Input compound synonym"),
                      selectInput("KW_synonym", "Select from list", choices = NULL),
                      actionButton("submit_compName","Submit")
                      )
             )
           ),
           box(width = 6,
               title = strong("Search Result:"),
               solidHeader = T,
               status = "primary",
               fluidRow(
                 div(
                   style = "margin:25px",
                   downloadButton("comp_report","Generate Report"),
                   hr(),
                   textOutput("summary_path"),
                   div(style = "height:300px;overflow-x:auto;overflow-y:scroll",
                       helpText("Preview of output only display first 20 items."),
                       DT::dataTableOutput("result3")
                   )
                 )
               )
           ),
           hr(),
           fluidRow(
             box(
               width = 12,
               solidHeader = T,
               status = "info",
               collapsible = T,
               collapsed = F,
               title = strong("Summary"),
               uiOutput("preview_tab3")
             )
           )
         ),
         tabPanel(
           title = strong("Input list of metabolites"),
           label = "sub_tab_2_TAB3",
           fluidRow(
             box(width = 6,
                 solidHeader = T,
                 status = "primary",
                 title = strong("Input a list of metabolites:"),
                 textAreaInput("input_mul_tab3",label = "",
                               placeholder = "Input list of metabolites in lines or separated by \",\""),
                 actionButton("sub_mul_tab3",label = "Submit")
             ),
             box(
               width = 6,
               solidHeader = T,
               status = "primary",
               title = strong("Upload the file"),
               fileInput("inp_file_tab3",label = "",
                         multiple = TRUE,
                         accept = c("text/csv","text/comma-separated-values,/text/plain",".csv",".txt"),
                         buttonLabel = "Browse..."),
               actionButton("sub_file_tab3",label = "Upload")
             )
           ),
           hr(),
           fluidRow(
             HTML("<div id='database-group-output'>"),
             box(width = 12,
                 title = strong("Data preview"),
                 solidHeader = T,
                 status = "primary",
                 downloadButton("tab3_mul_report",label = "Download Table"),
                 hr(),
                 dataTableOutput("preview_multi_names"),
                 highchartOutput("tab3_hc_output"),
                 textOutput("hc_click_output"),
                 tableOutput("stats_fisher_tab3"),
                 conditionalPanel(condition = "output.stats_fisher_tab3 != null",
                                  selectInput("pvalue_fisher",
                                              "Select significant level for Fisher Exact Test",
                                              choices = c(0.01,0.05,0.1)),
                                  actionButton("generateFisherTest",
                                               "Pathway Enrichment Analysis"),
                                  downloadButton("stats_report","Download Fiser Test Result"),
                                  hr(),
                                  box(title = strong("Enriched pathways identified by Fisher Test"),
                                      collapsible = T,
                                      collapsed = T,
                                      solidHeader = T,
                                      width = 12,
                                      status = "info",
                                      column(width = 6,
                                        # htmlOutput("summary_fisher")
                                        textOutput("text_fisher"),
                                        DT::dataTableOutput("summary_fisher")
                                      ),
                                      column(width = 6,
                                             highchartOutput("heatmap_pvalue")
                                             )
                                  )
                 )
             ),
             HTML("</div>")
           )
         )
         )
  )
