tabItem_about <- shinydashboard::tabItem(
  tabName = "About",
  shinydashboard::box(width = 12,
      title = tags$h1(strong("RaMP: Relational Database of Metabolic Pathways")),
      footer = p(strong("Contact: Bofei Zhang zhang.5675@osu.edu")),
      solidHeader = T,
      div(
        style ="margin:25px;",
        hr(),
        p("For usage of this software, you must have RaMP databases imported on
          your local computer. For importing database, you should refer to the 
          repository ",
          tags$a(
            href = "https://bitbucket.org/baskine/mathelabramp",
            rel="noopener noreferrer",
            target = "_blank",
            "mathelabramp"
          ), " to build local RaMP databases."),
        p(" RaMP is a conglomerate of 4 different databases:"),
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
        p("The imported databases pulls together the information scattered across them 
          -- mainly compounds, genes, pathways, and ontology. The relationship between 
          metabolites, pathways and ontology are connected by their RaMP ID. Then, user
          can explore the end from given origins.For example, the current functions
          for this software include: "),
        tags$ul(
          id = "ramp-function",
          tags$li("Tab1: Given a metabolite's synonym, it returns all metabolites that 
                  are involved in same pathway that the metabolite has."),
          tags$li("Tab2: Given a pathway name, it returns all metabolites' synonyms that have this
                  pathway"),
          tags$li("Tab3: Given a metabolite synonym, it returns all pathways that contain this 
                  metabolite"),
          tags$li("Tab4: Given a synonym of metabolite, it returns metabolites based on catalyzation
                  relationship"),
          tags$li("Tab5: Given a synonym of metabolite or ontology, it returns ontology or metabolite
                  respectively")
          ),
        p("Based on content of RaMP databases, the overlap under specific condition could be determined 
          and shown in Venn Diagram."),
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