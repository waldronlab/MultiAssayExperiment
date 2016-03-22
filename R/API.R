.cleanC = function(x) gsub(".__C__", "", x)
.cleanT = function(x) gsub(".__T__", "", x)

#' Refer to the API documentation
#' 
#' \code{API} opens a browser to the API documentation
#' 
#' @param shiny (logical default FALSE) whether to launch the shiny version
#' of the API (experimental)
#' 
#' @return Documentation via the GitHub wiki
#' 
#' @examples
#' API()
#' 
#' @author Vincent J Carey
#' 
#' @import shinydashboard shiny
#' 
#' @export
API <- function(shiny = FALSE) {
  if (shiny) {
    .apiDash()
  } else {
    utils::browseURL(
"https://github.com/vjcitn/MultiAssayExperiment/wiki/MultiAssayExperiment-API"
    )
  }
}

.packageAPI <- function(packname="MultiAssayExperiment") {
  ns = getNamespace(packname)
  nsElements = ls(ns, all.names = TRUE)
  nsNS = get(".__NAMESPACE__.", ns) # probably volatile?
  classes = .cleanC(grep(".__C__", nsElements, value=TRUE))
  todrop = sapply(classes, function(x)extends(x, "language"))
  if (any(todrop)) classes = classes[-which(todrop)]
  TmethsWithin = .cleanT(grep(packname, grep(".__T__", nsElements,
                                             value=TRUE),
                              value=TRUE))
  TmethsWithout = .cleanT(grep(packname, grep(".__T__", nsElements,
                                              value=TRUE),
                               value=TRUE, invert=TRUE))
  Mcmeths = lapply(classes, function(x) utils::methods(class=x))
  names(Mcmeths) = classes
  #list(nsElements = nsElements, nsNS = nsNS, classes=classes, cmeths=cmeths)
  list(classes=classes, Mcmeths=Mcmeths, TmethsWithin=TmethsWithin,
       TmethsWithout=TmethsWithout)
}

.apiDash <- function() {
  ui <- shinydashboard::dashboardPage(
    shinydashboard::dashboardHeader(title = "Basic dashboard"),
    shinydashboard::dashboardSidebar(
      shiny::textInput("packname", "Package", "MultiAssayExperiment")
    ),
    shinydashboard::dashboardBody(
      # Boxes need to be put in a row (or column)
      #    fluidRow(
      #      box(plotOutput("plot1", height = 250)),
      #
      #      box(
      #        title = "Controls",
      #        sliderInput("slider", "Number of observations:", 1, 100, 50)
      #      )
      #    )
      shiny::fluidRow(
        shinydashboard::tabBox(
          title = "API explorer",
          # The id lets us use input$tabset1 on the server to
          # find the current tab
          id = "tabset1", height = "250px",
          shiny::tabPanel("Classes", 
                   shiny::tableOutput("ls1")), 
          shiny::tabPanel("Methods(inside)", shiny::tableOutput("ls2")),
          shiny::tabPanel("Methods(outside)", shiny::tableOutput("ls3")),
          width=9
        )
      )
    )
  )
  
  server <- function(input, output) {
    
    getconts = shiny::reactive( {
      papi = .packageAPI( input$packname )
      list(classDF = data.frame(classes=papi$classes),
           methsInDF = data.frame(methsIn=papi$TmethsWithin),
           methsOutDF = data.frame(methsOut=papi$TmethsWithout))
    })
    output$ls1 = shiny::renderTable( getconts()$classDF )
    output$ls2 = shiny::renderTable( getconts()$methsInDF )
    output$ls3 = shiny::renderTable( getconts()$methsOutDF )
  }
  
  shiny::shinyApp(ui, server)
}
