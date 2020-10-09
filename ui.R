# Import libraries
library(shiny)
library(shinythemes)
library(ggplot2)
library(ggiraph)

options(stringsAsFactors = FALSE)

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
                             sliderInput(inputId = "Sepal.Length", label = "Sepal Length", value = 5.0,
                                         min = min(TrainSet$Sepal.Length),
                                         max = max(TrainSet$Sepal.Length)
                             ),
                             sliderInput(inputId = "Sepal.Width", label = "Sepal Width", value = 3.6,
                                         min = min(TrainSet$Sepal.Width),
                                         max = max(TrainSet$Sepal.Width)),
                             sliderInput(inputId = "Petal.Length", label = "Petal Length", value = 1.4,
                                         min = min(TrainSet$Petal.Length),
                                         max = max(TrainSet$Petal.Length)),
                             sliderInput(inputId = "Petal.Width", label = "Petal Width", value = 0.2,
                                         min = min(TrainSet$Petal.Width),
                                         max = max(TrainSet$Petal.Width)),

                             actionButton(inputId = "submitbutton", "Submit", class = "btn btn-primary")


                           ),

                           mainPanel(
                             tags$label(h3("Status/Output")), # Status/Output Text Box, same as HTML
                             verbatimTextOutput("contents"),
                             girafeOutput("dotplot_girafe"),
                             plotOutput("UMAP_plot")
                           ) # mainPanel

                  ), #tabPanel(), Home

                  tabPanel("About",
                           titlePanel("About"),
                           div(includeMarkdown("about.md"),
                               align="justify")
                  ) #tabPanel(), About

                ) # navbarPage
) # fluidPage

