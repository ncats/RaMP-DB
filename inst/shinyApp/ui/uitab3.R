tabItem3 <-  shinydashboard::tabItem(
  tabName = "pathFromMeta",
  shinydashboard::tabBox(width = 12, id = "tabset3",
         shiny::tabPanel(
           title = strong("Input analyte one by one"),
           label ="sub_tab_1_TAB3",
           shinydashboard::box(
             width = 6,
             title = strong("Input synonym"),
             solidHeader = T,
             status = "primary",
             fluidRow(
               column(12,
                      helpText("Given compound's name, 
                               it returns pathways in which the compound get involved in."),
                      textInput("compName","",placeholder = "Input compound synonym or source"),
                      radioButtons("synonymOrSource","Search by synonyms or source IDs?", choices = c(
                        "Synonym" = "synonyms", "Source ID" = "ids"
                      ),selected = "ids"),
                      selectInput("KW_synonym", "Select from list", choices = NULL),
                      actionButton("submit_compName","Submit")
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
                   downloadButton("comp_report","Download Results"),
                   hr(),
                   textOutput("summary_path")
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
#               uiOutput("preview_tab3")
                   div(style = "height:900px;overflow-x:auto;overflow-y:scroll",
                       helpText("Preview of output only display first 20 items."),
                       DT::dataTableOutput("result3")
                   )
             )
           )
         ),
         shiny::tabPanel(
           title = strong("Input multiple analytes (batch query)"),
           label = "sub_tab_2_TAB3",
           fluidRow(
	     #HTML("<div id='database-group-output'>"),
             shinydashboard::box(width = 6,
                 solidHeader = T,
                 status = "primary",
                 title = strong("Input a list of genes or metabolites, or a file"),
		 h4("Input list of metabolites, one per line:"),
                 textAreaInput("input_mul_tab3",label = "",
                               placeholder = "Input list of analytes, one per line"),
		 radioButtons("synonymOrSourcemult","Search by synonyms or source IDs?", choices = c(
                        "Synonym" = "synonyms", "Source ID" = "ids"
                      ),selected = "ids"),
                 actionButton("sub_mul_tab3",label = "Submit"),
		 br(),
                 br(),
		 h4("Upload file, one analyte per line"),
                 fileInput("inp_file_tab3",label = "",
                         multiple = TRUE,
                         accept = c("text/csv","text/comma-separated-values,/text/plain",".csv",".txt"),
                         buttonLabel = "Browse..."),
                 actionButton("sub_file_tab3",label = "Upload")
             ), # end box
             shinydashboard::box(width = 6,
                 title = strong("Summary:"),
                 solidHeader = T,
                 status = "primary",
		 p("The query returned the following number of pathways per query gene/metabolite:"),
		 DT::dataTableOutput("summary_mulpath_out"),
		 p("Query results should be visible here"),
                 downloadButton("tab3_mul_report",label = "Download Results")
	      )
	   #HTML("</div>")
         ), # end fluid row
          fluidRow(
	     HTML("<div id='database-group-output'>"),
             shinydashboard::box(width = 12,
                 title = strong("Search Results:"),
                 solidHeader = T,
                 status = "primary",
                 DT::dataTableOutput("preview_multi_names")
              ),
             HTML("</div>")
           ), # end fluidRow
          fluidRow(
            #shinydashboard::box(
            #  width = 6,
            #  shiny::fileInput(
            #    "fisher_test_background_tab3_2",
            #     label = "Input Background for Fisher Test",
            #     placeholder = "Specific format required",
            #     buttonLabel = "Upload"
            #     )
            #  ),
            shinydashboard::box(
                 width = 6,
                  numericInput("pvalue_fisher",
                        "Select significant level for Fisher Exact Test",
                         value = 0.01,
			 min=0,max=1,
                         width = "80%"),
		   actionButton("generateFisherTest",
                     "Pathway Enrichment Analysis")
                )
              )
             #),
#             conditionalPanel(condition = "output.preview_multi_names != null",
#                              highcharter::highchartOutput("tab3_hc_output"),
#                              textOutput("hc_click_output"),
#                              tableOutput("stats_fisher_tab3"),
#                              downloadButton("stats_report","Download Fiser Test Result"),
#                              hr(),
#                              shinydashboard::box(title = strong("Enriched pathways identified by Fisher Test"),
#                                                  collapsible = T,
#                                                  collapsed = T,
#                                                  solidHeader = T,
#                                                  width = 12,
#                                                  status = "info",
#                                                  column(width = 6,
#                                                         # htmlOutput("summary_fisher")
#                                                         textOutput("text_fisher"),
#                                                         DT::dataTableOutput("summary_fisher")
#                                                  ),
#                                                  column(width = 6,
#                                                         highcharter::highchartOutput("heatmap_pvalue")
#                                                  )
#                              )
#             )
#           ),

        ) # end tabPanel
) # end Tab Box
) # end Tabitem
