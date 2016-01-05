cleanC = function(x) gsub(".__C__", "", x)
cleanT = function(x) gsub(".__T__", "", x)

packageAPI = function(packname="biocMultiAssay") {
  ns = getNamespace(packname)
  nsElements = ls(ns, all.names = TRUE)
  nsNS = get(".__NAMESPACE__.", ns) # probably volatile?
  classes = cleanC(grep(".__C__", nsElements, value=TRUE))
  todrop = sapply(classes, function(x)extends(x, "language"))
  if (any(todrop)) classes = classes[-which(todrop)]
  TmethsWithin = cleanT(grep(packname, grep(".__T__", nsElements, value=TRUE), value=TRUE))
  TmethsWithout = cleanT(grep(packname, grep(".__T__", nsElements, value=TRUE), value=TRUE, invert=TRUE))
  Mcmeths = lapply(classes, function(x) methods(class=x))  # use methods()
  names(Mcmeths) = classes
  #list(nsElements = nsElements, nsNS = nsNS, classes=classes, cmeths=cmeths)
  list(classes=classes, Mcmeths=Mcmeths, TmethsWithin=TmethsWithin,
       TmethsWithout=TmethsWithout)
}

apiDash = function() {
  requireNamespace(shinydashboard)
  ui <- dashboardPage(
    dashboardHeader(title = "Basic dashboard"),
    dashboardSidebar(
      textInput("packname", "Package", "biocMultiAssay")
    ),
    dashboardBody(
      # Boxes need to be put in a row (or column)
      #    fluidRow(
      #      box(plotOutput("plot1", height = 250)),
      #
      #      box(
      #        title = "Controls",
      #        sliderInput("slider", "Number of observations:", 1, 100, 50)
      #      )
      #    )
      fluidRow(
        tabBox(
          title = "API explorer",
          # The id lets us use input$tabset1 on the server to find the current tab
          id = "tabset1", height = "250px",
          tabPanel("Classes", 
                   tableOutput("ls1")), 
          tabPanel("Methods(inside)", tableOutput("ls2")),
          tabPanel("Methods(outside)", tableOutput("ls3")),
          width=9
        )
      )
    )
  )
  
  server <- function(input, output) {
    
    getconts = reactive( {
      papi = packageAPI( input$packname )
      list(classDF = data.frame(classes=papi$classes),
           methsInDF = data.frame(methsIn=papi$TmethsWithin),
           methsOutDF = data.frame(methsOut=papi$TmethsWithout))
    })
    output$ls1 = renderTable( getconts()$classDF )
    output$ls2 = renderTable( getconts()$methsInDF )
    output$ls3 = renderTable( getconts()$methsOutDF )
  }
  
  shinyApp(ui, server)
}
