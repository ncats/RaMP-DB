
source(file.path("ui","header.R"),local = T)$value

source(file.path("ui","side.R"),local = T)$value

source(file.path("ui","uitababout.R"),local = T)$value

#source(file.path("ui","uitab1.R"),local = T)$value

source(file.path("ui","uitab2.R"),local = T)$value

source(file.path("ui","uitab3.R"),local = T)$value

source(file.path("ui","uitab4.R"),local = T)$value

#source(file.path("ui","uitab5.R"),local = T)$value

#source(file.path("ui","uitab6.R"),local = T)$value
body <- shinydashboard::dashboardBody(
  tags$head(
    tags$link(
      rel = "stylesheet", type = "text/css", href ="rampstyle.css"
    )
  ),
  shinydashboard::tabItems(
    tabItem_about,
    #tabItem1,
    tabItem2,
    tabItem3,
    tabItem4
    #tabItem5,
    #tabItem6
  )
) 




ui <- shinydashboard::dashboardPage(header, sideBar, body)


