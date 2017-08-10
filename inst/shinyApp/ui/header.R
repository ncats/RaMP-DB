
header <- dashboardHeader(title = "RaMP", 
                          titleWidth = 350,
                          dropdownMenu(type = "messages",
                                       messageItem(from = "Author",
                                                   message = p("Query need to take some time if", 
                                                   br(), " number of elements is too high (>10,000)"), 
                                                   icon = icon(" fa-hand-paper-o"), 
                                                   time = date()), 
                                       messageItem(from = "Author", 
                                                   message = p("Preview of table display maximum 20 items from query"), 
                                                   icon = icon(" fa-hand-paper-o"), 
                                                   time = date())
                                       )
                          )
