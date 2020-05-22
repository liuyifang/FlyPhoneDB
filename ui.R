# Import libraries
library(shiny)
library(shinythemes)
library(data.table)
library(randomForest)

# Read in the RF model
model <- readRDS("model.rds")

# Training set
TrainSet <- read.csv("training.csv", header = TRUE)
TrainSet <- TrainSet[,-1]

####################################
# User interface                   #
####################################

# Define UI
ui <- fluidPage(theme = shinytheme("yeti"),
                navbarPage(
                  # theme = "cerulean",  # <--- To use a theme, uncomment this
                  "CLP",
                  tabPanel("Iris Predictor",

                           # Input values
                           sidebarPanel(
                             #HTML("<h3>Input parameters</h3>"),
                             tags$label(h3("Input parameters")),
                             sliderInput("Sepal.Length", label = "Sepal Length", value = 5.0,
                                         min = min(TrainSet$Sepal.Length),
                                         max = max(TrainSet$Sepal.Length)
                             ),
                             sliderInput("Sepal.Width", label = "Sepal Width", value = 3.6,
                                         min = min(TrainSet$Sepal.Width),
                                         max = max(TrainSet$Sepal.Width)),
                             sliderInput("Petal.Length", label = "Petal Length", value = 1.4,
                                         min = min(TrainSet$Petal.Length),
                                         max = max(TrainSet$Petal.Length)),
                             sliderInput("Petal.Width", label = "Petal Width", value = 0.2,
                                         min = min(TrainSet$Petal.Width),
                                         max = max(TrainSet$Petal.Width)),

                             actionButton("submitbutton", "Submit", class = "btn btn-primary")
                           ),

                           mainPanel(
                             tags$label(h3("Status/Output")), # Status/Output Text Box, same as HTML
                             verbatimTextOutput("contents"),
                             tableOutput("tabledata") # Prediction results table
                           ) # mainPanel

                  ), #tabPanel(), Home

                  tabPanel("About",
                           titlePanel("About"),
                           div(includeMarkdown("about.md"),
                               align="justify")
                  ) #tabPanel(), About

                ) # navbarPage
) # fluidPage

####################################
# Server                           #
####################################

# Define server function
server <- function(input, output, session) {

  # Input Data
  datasetInput <- reactive({

    df <- data.frame(
      Name = c("Sepal Length",
               "Sepal Width",
               "Petal Length",
               "Petal Width"),
      Value = as.character(c(input$Sepal.Length,
                             input$Sepal.Width,
                             input$Petal.Length,
                             input$Petal.Width)),
      stringsAsFactors = FALSE)

    Species <- 0
    df <- rbind(df, Species)
    input <- transpose(df)
    write.table(input,"input.csv", sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE)

    test <- read.csv(paste("input", ".csv", sep=""), header = TRUE)

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







