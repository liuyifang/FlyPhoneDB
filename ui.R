# Import libraries
library(shiny)
library(shinythemes)
library(data.table)

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

