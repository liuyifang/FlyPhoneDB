# Import libraries
library(shiny)
# library(shinythemes)
# library(shinydashboard)
library(semantic.dashboard)
library(ggplot2)
library(ggiraph)
library(dplyr)
library(tidyr)
library(ggdendro)

options(stringsAsFactors = FALSE)

# Define UI
ui <- dashboardPage(
                dashboardHeader(title = "CLP"),
                dashboardSidebar(
                  sidebarMenu(
                    menuItem("MT", tabName = "mt"),
                    menuItem("Heart", tabName = "heart")
                  )
                ),
                dashboardBody(
                  tabItems(
                    tabItem("mt",
                            box(girafeOutput("heatmap2_girafe"), width = 12),
                            box(selectInput("heatmap_correlation", "Cluster row and column:", c(TRUE, FALSE)), width = 4),
                            # box(verbatimTextOutput("heatmap2_choices")),
                            box(girafeOutput("dotplot_girafe"), width = 16),
                            # box(verbatimTextOutput("dotplot_choices")),
                            box(plotOutput("UMAP_plot"))
                            ),
                    tabItem("heart",
                            h1("Heart")
                    )
                  )
                )
)

