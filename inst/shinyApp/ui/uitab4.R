# Tab 4
# 
# 
tabItem4<-  shinydashboard::tabItem(
  tabName = "geneCataComp",
  shinydashboard::tabBox(width = 12,id = "tabset4",height = "100%",
         shiny::tabPanel(strong("Input metabolite or gene one by one"),
                  "",
                  shinydashboard::box(width=6,
                      title = strong("Input metabolite or gene"),
                      solidHeader = T,
                      height = "100%",
                      status = "primary",
                      fluidRow(
                        column(12,
                               helpText("If a compound is input, 
		all genes that catalyze reactions involving the compound are returned."),
			       helpText("Conversely,if a gene is input, all compounds in reactions that are catalyzed by that gene are returned"),
                               textInput("CataInput","", placeholder = "Input compound name or id"),
                               radioButtons('CataInput_choices',label = 'Search by name or source id',
                                            choices = c('Names' = 'names','Source ID' = 'ids')),
                               selectInput("KW_cata", "Select from list", choices = NULL),
                               actionButton("subText_cata","Submit")
                        ) # end column
                      ) # end fluidRow
                    ),# end box
                  shinydashboard::box(width = 6,title = strong("Search Result:"),solidHeader = T,
                      status = "primary",
                      fluidRow(
                        div(style = "margin:25px;",
                            downloadButton("report_cata","Download Results"),
                            hr(),
                            textOutput("summary_cata"),
                            DT::dataTableOutput("result_cata")
                        )
                      ) #end fluidRow
                  ), # end box
                  hr(),
#                   fluidRow(
#                     shinydashboard::box(
#                       width = 12,
#                       solidHeader = T,
#                       status = "info",
#                       collapsible = T,
#                       collapsed = F,
#                       title = strong("Summary"),
# 		      DT::dataTableOutput("result_cata")
#                       #uiOutput("preview_tab4")
#                     ) # end box
# 		  ), # end fluid row
		  hr(),
		  fluidRow(
                    shinydashboard::box(
                      width = 12,
		      height= "1000px",
                      solidHeader = T,
                      status = "primary",
                      collapsible = T,
                      collapsed = F,
                      title = strong("Visuazlize gene-metabolite interaction network"),
                      visNetwork::visNetworkOutput("network",height = '1000px')
                    ) # end box
		   ) # end fluidRow
         ),
         shiny::tabPanel(
           title = strong("Input a list of genes or metabolites"),
           fluidRow(
             shinydashboard::box(width = 6,
                 solidHeader = T,
                 status = "primary",
                 title = strong("Input a list of genes or metabolites:"),
                 textAreaInput("input_mul_tab4",label = "",
                               placeholder = "Input list of genes or metabolites, one per line"),
                 radioButtons('input_mul_tab4_choices',label = 'Search by source id or name?',
                              choices = c('Name' = 'names','Source ID' = 'ids')),
                 actionButton("sub_mul_tab4",label = "Submit")
             ), # end box
             shinydashboard::box(width = 6,
                 solidHeader = T,
                 status = "primary",
                 title = strong("Search Result"),
		 DT::dataTableOutput("preview_multi_names_tab4"),
		fluidRow(
			div(style = "margin:25px;",
 	 			downloadButton("tab4_mul_report",label = "Download Results"),
				hr()
			) # end div
		) # end fluidRow
	     ) # end box

#             shinydashboard::box(
#               width = 6,
#               title = strong("Upload File"),
#               solidHeader = T,
#               status = "primary",
#               fileInput("inp_file_tab4",label = "",
#                         multiple = FALSE,
#                         accept = c("text/csv","text/comma-separated-values,/text/plain",".csv",".txt"),
#                         buttonLabel = "Browse..."),
#               actionButton("sub_file_tab4",label = "Upload File and Run Query")
#             )
           ), # end fluidRow (with two boxes)
           hr(),
           fluidRow(
             #HTML("<div id='database-group-output'>"),
             #shinydashboard::box(width = 12,
             #    title = strong("Summary"),
             #    solidHeader = T,
             #    status = "primary",
             #    #downloadButton("tab4_mul_report",label = "Download Results"),
             #    #hr(),
             #    DT::dataTableOutput("preview_multi_names_tab4")
                 #highcharter::highchartOutput("tab4_hc_output")
             #), # end box
             #hr(),
             shinydashboard::box(
             width = 12,height="1000px",
            	solidHeader = T,
            	status = "primary",
            	collapsible = T,
            	collapsed = F,
            	title = strong("Visuazlize gene-metabolite interaction network"),
            	visNetwork::visNetworkOutput("networkmulti",height = '1000px')
             )
           ) # end 2nd fluidRow
         ) # end tab
     )# end tabBox
  )# end tabItem




  
