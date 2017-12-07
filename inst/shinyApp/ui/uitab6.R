
tabItem6 <- tabItem(
  tabName = "customizeQuery",
  shinydashboard::box(width= 12,
      title = strong("Provide Search Condition:"),
      selectInput("givenDataType","Define Type of Searching Condition"
                  ,choices = c("metabolites",
                              "genes",
                              "pathways",
                              "biofluid"),
                  selected = "metabolites"
                  ),
      textInput("kwInput",label = "Input key word",
                placeholder ="Input kewy word to see if the database has this item."),
      selectInput("kwSearch","",choices = NULL,multiple = T),
      actionButton("submitKW","Submit"),
      textAreaInput("givenData",label = "Search",
                    placeholder = "Input a list of selected type"),
      HTML("<div style='inline-block'>"),
      actionButton("searchGivenData","Search"),
      actionButton("clearGivenData","Clear"),
      HTML("</div>"),
      conditionalPanel(
        "input.givenDataType == 'metabolites' || input.givenDataType == 'genes'",
        selectInput("givenMetabolites",
                    label = "Find ... ",
                    choices = c("Synonym",
                                "Source",
                                "Biofluid",
                                "Catalyzation",
                                "Pathway"))
      ),
      conditionalPanel(
        "input.givenDataType == 'pathways'",
        selectInput("givenPathways",
                    label = "Find ...",
                    c("metabolites","genes",
                      "metabolites+genes","source"))
      ),
      conditionalPanel(
        "input.givenDataType == 'biofluid'",
        selectInput("givenBiofluid",
                    label = "Find ...",
                    c("metabolites","genes","metabolites+genes")
                    )
      ),
      DT::dataTableOutput("searchingResult")
  )
)