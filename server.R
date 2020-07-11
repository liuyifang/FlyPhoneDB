# Import libraries
library(shiny)
library(data.table)
library(randomForest)

options(stringsAsFactors = FALSE)

# Read in the RF model
model <- readRDS("model.rds")

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

  # Prediction results table
  output$tabledata_output <- renderTable({
    if (input$submitbutton>0) {
      isolate(datasetInput())
    }
  })

  # Export dotplot table
  output$dotplot_table_output <- renderTable({
    if (input$submitbutton>0) {
      isolate(dotplot_table())
    }
  })

} # server
