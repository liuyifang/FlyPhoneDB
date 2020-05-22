# Load R packages
library(shiny)
library(shinythemes)


  # Define UI
  ui <- fluidPage(theme = shinytheme("yeti"),
    navbarPage(
      # theme = "cerulean",  # <--- To use a theme, uncomment this
      "My first app",
      tabPanel("Navbar 1",
               sidebarPanel(
                 tags$h3("Input:"),
                 textInput("txt1", "Given Name:", "hello"), # txt1 will be send to server
                 textInput("txt2", "Surname:", "mao"),    # txt2 will be send to server

               ), # sidebarPanel
               mainPanel(
                            h1("Header 1"),

                            h4("Output 1"),
                            verbatimTextOutput("txtout"), # txtout

               ) # mainPanel

      ), # Navbar 1, tabPanel
      tabPanel("Navbar 2", "This panel is intentionally left blank"),
      tabPanel("Navbar 3", "This panel is intentionally left blank")

    ) # navbarPage
  ) # fluidPage


  # Define server function
  server <- function(input, output) {

    output$txtout <- renderText({
      paste( input$txt1, input$txt2, sep = " " )
    })
  } # server


  # Create Shiny object
  shinyApp(ui = ui, server = server)
