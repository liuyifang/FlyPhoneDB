# Import libraries
library(shiny)
library(data.table)
library(randomForest)

options(stringsAsFactors = FALSE)

####################################
# Server                           #
####################################

# Define server function
server <- function(input, output, session) {

  # dotplot table
  dotplot_table <- reactive({
    dotplot_table_df <- read.csv("LR_pairs3.csv", row.names = 1)
    print(dotplot_table_df)
  })

  # Status/Output Text Box
  output$contents <- renderPrint({
    if (input$submitbutton>0) {
      isolate("Calculation complete.")
    } else {
      return("Server is ready for calculation.")
    }
  })

  # Export dotplot table
  output$dotplot_table_output <- renderTable({
    if (input$submitbutton>0) {
      isolate(dotplot_table())
    }
  })

} # server
