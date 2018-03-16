tabItem_about <- shinydashboard::tabItem(
  tabName = "About",
  shinydashboard::box(width = 12,
      title = tags$h1(strong("RaMP: Relational Database of Metabolic Pathways")),
      footer = p("For questions or issues, please submit an issue on our GitHub site at https://github.com/Mathelab/RaMP-DB"),
      solidHeader = T,
      div(
        style ="margin:25px;",
        p("To use this software on your local machine, you must have the RaMP SQL database installed locally.  For instructions on how to do so, please refer to the RaMP-DB GitHub repository ",
          tags$a(
            href = "https://github.com/Mathelab/RaMP-DB",
            rel="noopener noreferrer",
            target = "_blank",
            "RaMP-DB GitHub repository"
          ), "."),
        p(" RaMP currently integrates information from 4 different databases:"),
        tags$ul(
          id = "tab1link",
          tags$li(tags$a(href = "http://www.genome.jp/kegg/",
                         rel="noopener noreferrer",
                         target = "_blank",
                         "Kegg")),
          tags$li(tags$a(href = "http://www.reactome.org/",
                         rel="noopener noreferrer",
                         target = "_blank",
                         "Reactome")),
          tags$li(tags$a(href = "http://www.hmdb.ca/",
                         rel="noopener noreferrer",
                         target = "_blank",
                         "HMDB")),
          tags$li(tags$a(href = "http://www.wikipathways.org/index.php/WikiPathways",
                         rel="noopener noreferrer",
                         target = "_blank",
                         "wikipathways"))
        ),
        p("Of note, RaMP does not import all the information contained in these databases.  Instead, RaMP includes pathway names and the genes and metabolites they include.  RaMP is a work in progress and we will be expanding its content (reaction-level information, non-human compounds and genes, etc.).
          Currently, the following queries are supported through this interface: "),
        tags$ul(
          id = "ramp-function",
          tags$li("Tab1: Given one or multiple  pathway names, retrieve all genes and/or
		metabolites contained in the pathway(s)"),
          tags$li("Tab2: Given a list of metabolite(s) or gene(s), retrieve all pathways
		that they are involved in. This tab also support pathway 
		overrepresentation analysis"),
          tags$li("Tab3: Given one or multiple metabolite(s) or gene(s), retrieve all 
		gene(s) or metabolite(s), respectively, that are involved in the same
		reaction. A network of gene-metabolite relationships can be drawn."),
	  tags$li("Tab4: Given one or multiple metabolite(s), retrieve the ontology
		(e.g. biofluid location, cellular location, etc.) they belong in.")
          ),
	HTML("<p>More details can be found in <a href='http://www.mdpi.com/2218-1989/8/1/16' target='_blank'>our manuscript</a></p>"),
        p("The Venn diagram below shows the overlap in metabolites and genes from each database source (as of 01/24/2017). "), 
        HTML(
          "
          <div id = 'figure-container' >
          <a href ='#'>
          <figure>
          <img src = 'fourVennCompound.jpg' width ='200px' height = '200px' title = 'The purpose of plot is to compare information among the four databases that make up the RaMP databases -- wikipathways, reactome, kegg, hmdb. There may be overlaps between two or more of the databases. For example, one metabolite may be present in two databaase. By linking RAMPID to each compound's ID in the databases, the overlap in different databases is shown in this plot.'/>
          <figcaption>Compound Overlap in four databases</figcaption>
          </figure>
          </a>
          <a href ='#'>
          <figure>
          <img src = 'fourVennGene.jpg' width ='200px' height = '200px'title = 'The plot identifies the overlap of gene among the four databases that make up the RaMP databases -- wikipathways, reactome, kegg, hmdb.'/>
          <figcaption>Gene Overlap in four databases</figcaption>
          </figure>
          </a>
          </div>
          "
        )
        )
    )
  )
