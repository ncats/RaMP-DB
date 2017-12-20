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
                      radioButtons("NameOrId","Search by common names or source IDs?", choices = c(
                        "Names" = "names", "Source ID" = "ids"
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
		 radioButtons("NameOrSourcemult","Search by synonyms or source IDs?", choices = c(
                        "Names" = "names", "Source ID" = "ids"
                      ),selected = "ids"),
                 actionButton("sub_mul_tab3",label = "Submit"),
		 br()
#                 br(),
#		 h4("Upload file, one analyte per line"),
#                 fileInput("inp_file_tab3",label = "",
#                         multiple = TRUE,
#                         accept = c("text/csv","text/comma-separated-values,/text/plain",".csv",".txt"),
#                         buttonLabel = "Browse..."),
#                 actionButton("sub_file_tab3",label = "Upload")
             ), # end box
             shinydashboard::box(width = 6,
                 title = strong("Summary:"),
                 solidHeader = T,
                 status = "primary",
		 p("The query returned the following number of pathways per query gene/metabolite:"),
		 DT::dataTableOutput("summary_mulpath_out"),
		 p("Query results should be visible here"),
		 br(),
                 br(),
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
            shinydashboard::box(
              solidHeader = T,
              status = "primary",
           	  title = strong("Set Parameters to run Pathway Enrichment Analysis"),
	          width = 6,
                  numericInput("p_fdradj_cutoff",
                        "Select cutoff for FDR Adjusted p-values (will only return pathways with FDR adj p-values < cutoff)",
                         value = 0.01,
			 min=0,max=1,
                         width = "80%"),
		  selectInput("analyte_type", "Type of analyte",
			choices = c(
                            "Metabolites" = "metabolites",
			    "Genes" = "genes")),
			h4("Set parameters for functional pathway clustering:"),
			numericInput("perc_analyte_overlap", "Overlap threshold (0-1) for pathways to be considered similar (for medoid establishment):",min = 0.01, max = 1, value = 0.2, step = 0.01),
			numericInput("min_pathway_tocluster", "Number of similar neighbors required (for medoid establishment):",min = 1, max = 100, value = 2, step = 1),
			numericInput("perc_pathway_overlap", "Overlap threshold (0-1) for pathways to be clustered:",min = 0.01, max = 1, value = 0.2, step = 0.01),

		  #numericInput("total_genes",
		  #	"Input the total number of genes measured in experiment (to be used as background)",
		  #	value=20000,
		  #	min=1,max=100000,width="80%"),
		   #numericInput("total_metabolites",
                  #      "Input the total number of metabolites measured in experiment (to be used as background)",
                  #      value=1000,
                  #      min=1,max=100000,width="80%"),
		actionButton("runFisher","Run Pathway Enrichment (please be patient!)")
                ),#end of Box
            shinydashboard::box(
                 width = 6,
                 solidHeader = T,
                 status = "primary",
		 title = strong("Summary:"),
		 p("Significant pathways are returned below and can be downloaded by clicking 'Download Results'"),
		 p("Note that only pathways that contain at least 2 analytes from the user input will be output"),
                #textOutput("summary_Fisher"),
		DT::dataTableOutput("summary_fisher"),
                 downloadButton("fisher_stats_report",label = "Download Results")
              ) #end box
	  ), # end of fluidRow
	  fluidRow(
            shinydashboard::box(
                  title=strong("Results of Pathway Enrichment Analysis"),
                  solidHeader = T,
                  status = "primary",
                  width = 12,
                  column(width = 6,
          selectInput("show_cluster","Display pathways in cluster:",choices = 1),
          textOutput("cluster_summary_text")
                  ),
          column(width = 6,plotOutput("cluster_summary_plot",height = "300px")),
	    	  DT::dataTableOutput("results_fisher")
	    ) # end box
          ) # end of fluidRow
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
