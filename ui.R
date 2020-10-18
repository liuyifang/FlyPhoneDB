# Import libraries
library(shiny)
library(semantic.dashboard)
library(ggplot2)
library(ggiraph)
library(dplyr)
library(tidyr)
library(ggdendro)
library(DT)
library(DiagrammeR)
library(tidyverse)

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
                            box(selectInput("heatmap_correlation",
                                            "Cluster row and column:",
                                            c(TRUE, FALSE)),
                                width = 4),
                            box(girafeOutput("dotplot_girafe"), width = 10),
                            box(dataTableOutput("dotplot_table"), style = "height:500px; overflow-x: scroll;", width = 6), # https://stackoverflow.com/questions/47505893/adding-a-vertical-and-horizontal-scroll-bar-to-the-dt-table-in-r-shiny
                            # box(girafeOutput("UMAP_cluster_girafe"), width = 8),
                            # box(girafeOutput("UMAP_plot"), width = 8, style = "height:500px")
                            box(grVizOutput("networkPlot"), width = 12),
                            box(selectInput("network_layout",
                                            "layout:",
                                            c("nicely", "circle", "tree", "kk", "fr")),
                                width = 4)
                            ),
                    tabItem("heart",
                            h1("Heart")
                    )
                  )
                )
)

