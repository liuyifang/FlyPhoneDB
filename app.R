# Load R packages
library(shiny)
library(shinythemes)
data(airquality)

# Define UI
ui <- fluidPage(theme = shinytheme("yeti"),
                navbarPage(
                  # theme = "cerulean",  # <--- To use a theme, uncomment this
                  "My first app",
                  tabPanel("Ozone level!",
                           # Sidebar panel for inputs ----
                           sidebarPanel(

                             # Input: Slider for the number of bins ----
                             sliderInput(inputId = "bins",
                                         label = "Number of bins:",
                                         min = 0,
                                         max = 50,
                                         value = 30,
                                         step = 2)

                           ),
                           mainPanel(
                             # Output: Histogram ----
                             plotOutput(outputId = "distPlot")

                           ) # mainPanel

                  ), # Navbar 1, tabPanel
                  tabPanel("Navbar 2", "This panel is intentionally left blank"),
                  tabPanel("Navbar 3", "This panel is intentionally left blank")

                ) # navbarPage
) # fluidPage


# Define server function
server <- function(input, output) {

  output$distPlot <- renderPlot({

    x    <- airquality$Ozone
    x    <- na.omit(x)
    bins <- seq(min(x), max(x), length.out = input$bins + 1)

    hist(x, breaks = bins, col = "green", border = "black",
         xlab = "Ozone level",
         main = "Histogram of Ozone level")

  })

} # server

# Create Shiny object
shinyApp(ui = ui, server = server)
