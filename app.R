# Import libraries
library(shiny)
library(shinythemes)
library(data.table)
library(RCurl)
library(randomForest)

# # Read data
# weather <- read.csv(text = getURL("https://raw.githubusercontent.com/dataprofessor/data/master/weather-weka.csv") )
#
# # Build model
# model <- randomForest(play ~ ., data = weather, ntree = 500, mtry = 4, importance = TRUE)

# Save model to RDS file
# saveRDS(model, "model.rds")

# Read in the RF model
model <- readRDS("model.rds")
# print(model)

####################################
# User interface                   #
####################################

# Define UI
ui <- fluidPage(theme = shinytheme("yeti"),
                navbarPage(
                  # theme = "cerulean",  # <--- To use a theme, uncomment this
                  "CLP",
                  tabPanel("Play Golf?",
                           # Input values
                           sidebarPanel(
                             HTML("<h3>Input parameters</h3>"),

                             selectInput("outlook", label = "Outlook:",
                                         choices = list("Sunny" = "sunny", "Overcast" = "overcast", "Rainy" = "rainy"),
                                         selected = "Sunny"),
                             sliderInput("temperature", label = "Temperature:",
                                         min = 64, max = 86,
                                         value = 70),
                             sliderInput("humidity", label = "Humidity:",
                                         min = 65, max = 96,
                                         value = 90),
                             selectInput("windy", label = "Windy:",
                                         choices = list("Yes" = "TRUE", "No" = "FALSE"),
                                         selected = "TRUE"),

                             actionButton("submitbutton", "Submit", class = "btn btn-primary")
                           ),
                           mainPanel(
                             tags$label(h3("Status/Output")), # Status/Output Text Box
                             verbatimTextOutput("contents"),
                             tableOutput("tabledata") # Prediction results table
                           ) # mainPanel

                  ), # Navbar 1, tabPanel
                  tabPanel("Navbar 2", "This panel is intentionally left blank"),
                  tabPanel("Navbar 3", "This panel is intentionally left blank")

                ) # navbarPage
) # fluidPage

####################################
# Server                           #
####################################

# Define server function
server <- function(input, output, session) {

  # Input Data
  datasetInput <- reactive({

    # outlook,temperature,humidity,windy,play
    df <- data.frame(
      Name = c("outlook",
               "temperature",
               "humidity",
               "windy"),
      Value = as.character(c(input$outlook,
                             input$temperature,
                             input$humidity,
                             input$windy)),
      stringsAsFactors = FALSE)

    play <- "play"
    df <- rbind(df, play)
    input <- transpose(df)

    # export to file, this is a good way to debug
    write.table(input,"input.csv", sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE)

    test <- read.csv(paste("input", ".csv", sep=""), header = TRUE)

    test$outlook <- factor(test$outlook, levels = c("overcast", "rainy", "sunny"))


    Output <- data.frame(Prediction=predict(model,test), round(predict(model,test,type="prob"), 3))
    print(Output)

  })

  # Status/Output Text Box
  output$contents <- renderPrint({
    if (input$submitbutton>0) {
      isolate("Calculation complete.")
    } else {
      return("Server is ready for calculation.")
    }
  })

  # Prediction results table
  output$tabledata <- renderTable({
    if (input$submitbutton>0) {
      isolate(datasetInput())
    }
  })

} # server

####################################
# Create the shiny app             #
####################################

# Create Shiny object
shinyApp(ui = ui, server = server)
