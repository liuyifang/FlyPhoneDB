# Import libraries
library(shiny)
# library(shinythemes)
library(shinydashboard)
library(ggplot2)
library(ggiraph)
library(dplyr)
library(tidyr)
library(ggdendro)

options(stringsAsFactors = FALSE)

####################################
# User interface                   #
####################################

# Define UI
ui <- dashboardPage(
                dashboardHeader(title = "CLP"),
                dashboardSidebar(
                  sidebarMenu(
                    menuItem("MT", tabName = "mt", icon = icon("tree")),
                    menuItem("Heart", tabName = "heart", icon = icon("tree"))
                  )
                ),
                dashboardBody(
                  box(girafeOutput("heatmap2_girafe"), width = 8),
                  box(verbatimTextOutput("heatmap2_choices")),
                  box(girafeOutput("dotplot_girafe")),
                  box(verbatimTextOutput("dotplot_choices")),
                  box(plotOutput("UMAP_plot"))
                )
)

