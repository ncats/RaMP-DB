sideBar <-shinydashboard::dashboardSidebar(
  shinydashboard::sidebarMenu(
    shinydashboard::menuItem(
      "About",
      tabName = "About",
      icon = icon("fa fa-superpowers")
    ),
    # Given synonym, it returns compound or genes based on choices of user which have 
    # same pahtway involved.
    shinydashboard::menuItem(
      HTML("<p>Return analytes in pathways containing a given<br> analyte</p>"),
      tabName = "genesFromCompounds",
      icon = icon("dashboard")
    ),
    # Given pathway name, it returns genes or compound based on user choice which are involved
    # in that pathway.
    shinydashboard::menuItem(
      "Return analyte from given pathway name",
      tabName = "metaFromPath",
      icon = icon("dashboard")
    ),
    # Given metabolites' synonym, it returns pathway name in which the metabolites are 
    # involved in.
    shinydashboard::menuItem(
      "Return pathway from given analytes",
      tabName = "pathFromMeta",
      icon = icon("dashboard")
    ),
    shinydashboard::menuItem(
      "Return metabolites or genes based on catalyzation",
      tabName = "geneCataComp",
      icon = icon("dashboard")
    ),
    shinydashboard::menuItem(
      "Return metabolites or biofluid location",
      tabName = "geneCompOnto",
      icon = icon("dashboard")
    ),
    shinydashboard::menuItem(
      "Customize Query",
      tabName = "customizeQuery",
      icon = icon("dashboard")
    ),
    shinydashboard::menuItem(
        tags$button(
          id = "buttonstop",
          type = "button",
          class = "btn action-button",
          onclick = "setTimeout(function()
          {window.close();},50);",
          "Click to Exit RaMP"
        ),
        icon = icon("sign-out")
        )
    )
  ,width = 350
)
