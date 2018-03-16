tabItem3 <-  shinydashboard::tabItem(
  tabName = "pathFromMeta",
  shinydashboard::tabBox(width = 12, id = "tabset3",
         shiny::tabPanel(
           title = strong("Input analyte one by one"),
           label ="sub_tab_1_TAB3",
           shinydashboard::box(
             width = 6,
             title = strong("Input analyte"),
             solidHeader = T,
             status = "primary",
             fluidRow(
               column(12,
                      h5("Given a gene/metabolite name or ID,
                               retrieve pathways in which the compound is involved in."),
	              h5(strong("IMPORTANT NOTE about inputting source IDs:")),
		      h5("When inputting source IDs, it is important to add a prefix to denote the id type.  This is important because it is possible for two different metabolites to have the same IDs, although each ID may be from a different database source."),
                      h5("Metabolites can be searched with the following ID types: CAS, chebi, chemspider, hmdb, kegg, LIPIDMAPS, and pubchem.  To search for a metabolite, the ID type must be added as a prefix.  For example, the compound 'HMDB0000562' must be searched by 'hmdb:HMDB0000562', the compound '16737' must be searched by 'chebi:16737'."),
                      h5("Genes can be searched with the following ID types:  enzymeNomenclature, ensembl, entrez, hmdb, kegg, uniprot. Similar to metabolites, prefix ID types must be added to the ID for searching."),
		      br(),
                      column(width = 6),
                      textInput("compName","",placeholder = "Input compound name or source id"),
                      radioButtons("NameOrId","Search by common names or source IDs?",
			choices = c("Names" = "names", "Source ID" = "ids"),selected = "ids"),
                      selectInput("KW_synonym", "Select from list", choices = NULL),
                      actionButton("submit_compName","Submit")
               ) # end column
             )# end fluidRow
           ),# end box
           shinydashboard::box(width = 6,
               title = strong("Summary"),
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
           ), # end box
           hr(),
           shinydashboard::box(
               width = 12,
               solidHeader = T,
               status = "primary",
               collapsible = T,
               collapsed = F,
               title = strong("Query Result"),
	       fluidRow(
                   div(style = "width:80%;height:100%;align:center;margin-right:20px;margin-left:20px;",
                       helpText("Preview of output only display first 20 items."),
                       DT::dataTableOutput("result3")
                   )
		)
             )
         ), # end tabPanel
         shiny::tabPanel(
           title = strong("Input multiple analytes (batch query)"),
           label = "sub_tab_2_TAB3",
           fluidRow(
             shinydashboard::box(width = 6,
                 solidHeader = T,
                 status = "primary",
                 title = strong("1. Input a list of genes and/or metabolites"),
                 # h4(strong('Input a list of genes, one per line')),
                 shiny::mainPanel(
                   width = 12,
		   h5(strong("IMPORTANT NOTE about inputting source IDs:")), 
		   h5("When inputting source IDs, it is important to add a prefix to denote the id type.  This is important because it is possible for two different metabolites to have the same IDs, although each ID may be from a different database source."),
		   h5("Metabolites can be searched with the following ID types: CAS, chebi, chemspider, hmdb, kegg, LIPIDMAPS, and pubchem.  To search for a metabolite, the ID type must be added as a prefix.  For example, the compound 'HMDB0000562' must be searched by 'hmdb:HMDB0000562', the compound '16737' must be searched by 'chebi:16737'."),
		   h5("Genes can be searched with the following ID types:  enzymeNomenclature, ensembl, entrez, hmdb, kegg, uniprot. Similar to metabolites, prefix ID types must be added to the ID for searching."),
		   br(),
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
                   br()
                  ) # mainPanel
             ), # end box
             shinydashboard::box(width = 6,
                 title = strong("Summary"),
                 solidHeader = T,
                 status = "primary",
		 p("The following number of names/ids were mapped:"),
		 textOutput("num_mapped_namesids"),
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
                 title = strong("Query Result"),
                 solidHeader = T,
                 status = "primary",
                 div(style = 'width:80%;height:100%;align:center;',
                  DT::dataTableOutput("preview_multi_names")
                 )
              ), # end box
             HTML("</div>")
           ), # end fluidRow
	fluidRow(
	 	column(width = 6,
         		shinydashboard::box(
           			solidHeader = T,
		  	        status = "primary",
           			title = strong("2. Run pathway enrichment analysis"),
           			width = NULL,
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
           			textOutput("fishersProgress")

           		),#end of box

            shinydashboard::box(
           	solidHeader = T,
           	status = "primary",
           	title = strong("3. Set parameters for signifance filtering and functional pathway clustering:"),
           	width = NULL,
		br(),
		br(),
           	numericInput("p_holmadj_cutoff",
                        "Select cutoff for Holm Adjusted p-values (will only return pathways with Holm adj p-values < cutoff)",
                        value = 0.01,
                        min=0,max=1,
                        width = "80%"),
           	h4("Clustering Parameters:"),
           	numericInput("perc_analyte_overlap", "Overlap threshold (0-1) for pathways to be considered similar (for medoid establishment):",min = 0.01, max = 1, value = 0.2, step = 0.01),
           	numericInput("min_pathway_tocluster", "Number of similar neighbors required (for medoid establishment):",min = 1, max = 100, value = 2, step = 1),
           	numericInput("perc_pathway_overlap", "Overlap threshold (0-1) for pathways to be clustered:",min = 0.01, max = 1, value = 0.75, step = 0.01),
           	actionButton("runClustering","Filter and Cluster Results")
              ) # end box
	    ),#end of column

	  column(width = 6,
        	 shinydashboard::box(
        	   width = "100%",
        	   solidHeader = T,
        	   status = "primary",
        	   title = strong("4. Download results:"),
        	   p("Significant pathways are returned under 'Results of Pathway Enrichment Analysis' and can be downloaded by clicking 'Download Results'"),
        	   p("Note that only pathways that contain at least 3 analytes from the user input will be output"),
        	   #textOutput("summary_Fisher"),
               	   p("The following number of pathways per database were processed:"),
           	   DT::dataTableOutput("summary_fisher",width = '80%'),
           	   downloadButton("fisher_stats_report",label = "Download Results")
          	) # end box
	   )   #end column
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
                  ), # end columnb
          	  column(width = 6,
			plotOutput("cluster_summary_plot",height = "300px")
		   ), # end column
	    	   DT::dataTableOutput("results_fisher",width = '80%')
	    ) # end box
          ) # end of fluidRow
        ) # end tabPanel
  ) # end Tab Box
) # end Tabitem
