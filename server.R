# Import libraries
library(shiny)
library(shinythemes)
library(ggplot2)

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

  # dotplot data
  dotplot_data <- reactive({
    data <- read.csv("LR_pairs3.csv", row.names = 1)
    str(data)
    data$interacting_pair <- paste(data$gene_secreted, data$gene_receptor, sep = "_")

    means <- data[ , c("Lcluster0_Acluster0",
                       "Lcluster1_Acluster0",
                       "Lcluster2_Acluster0",
                       "Lcluster3_Acluster0",
                       "Lcluster4_Acluster0",
                       "Lcluster5_Acluster0",
                       "Lcluster6_Acluster0",
                       "Lcluster7_Acluster0",
                       "Lcluster8_Acluster0",
                       "Lcluster9_Acluster0",
                       "Lcluster10_Acluster0",
                       "Lcluster0_Acluster1",
                       "Lcluster1_Acluster1",
                       "Lcluster2_Acluster1",
                       "Lcluster3_Acluster1",
                       "Lcluster4_Acluster1",
                       "Lcluster5_Acluster1",
                       "Lcluster6_Acluster1",
                       "Lcluster7_Acluster1",
                       "Lcluster8_Acluster1",
                       "Lcluster9_Acluster1",
                       "Lcluster10_Acluster1",
                       "Lcluster0_Acluster2",
                       "Lcluster1_Acluster2",
                       "Lcluster2_Acluster2",
                       "Lcluster3_Acluster2",
                       "Lcluster4_Acluster2",
                       "Lcluster5_Acluster2",
                       "Lcluster6_Acluster2",
                       "Lcluster7_Acluster2",
                       "Lcluster8_Acluster2",
                       "Lcluster9_Acluster2",
                       "Lcluster10_Acluster2")]
    colnames(means) <- c("cluster0_Adipose",
                         "cluster1_Adipose",
                         "cluster2_Adipose",
                         "cluster3_Adipose",
                         "cluster4_Adipose",
                         "cluster5_Adipose",
                         "cluster6_Adipose",
                         "cluster7_Adipose",
                         "cluster8_Adipose",
                         "cluster9_Adipose",
                         "cluster10_Adipose",
                         "cluster0_Oenocyte",
                         "cluster1_Oenocyte",
                         "cluster2_Oenocyte",
                         "cluster3_Oenocyte",
                         "cluster4_Oenocyte",
                         "cluster5_Oenocyte",
                         "cluster6_Oenocyte",
                         "cluster7_Oenocyte",
                         "cluster8_Oenocyte",
                         "cluster9_Oenocyte",
                         "cluster10_Oenocyte",
                         "cluster0_Muscle",
                         "cluster1_Muscle",
                         "cluster2_Muscle",
                         "cluster3_Muscle",
                         "cluster4_Muscle",
                         "cluster5_Muscle",
                         "cluster6_Muscle",
                         "cluster7_Muscle",
                         "cluster8_Muscle",
                         "cluster9_Muscle",
                         "cluster10_Muscle")

    pval <- data[ , c("Lcluster0_Acluster0_pvalues",
                      "Lcluster1_Acluster0_pvalues",
                      "Lcluster2_Acluster0_pvalues",
                      "Lcluster3_Acluster0_pvalues",
                      "Lcluster4_Acluster0_pvalues",
                      "Lcluster5_Acluster0_pvalues",
                      "Lcluster6_Acluster0_pvalues",
                      "Lcluster7_Acluster0_pvalues",
                      "Lcluster8_Acluster0_pvalues",
                      "Lcluster9_Acluster0_pvalues",
                      "Lcluster10_Acluster0_pvalues",
                      "Lcluster0_Acluster1_pvalues",
                      "Lcluster1_Acluster1_pvalues",
                      "Lcluster2_Acluster1_pvalues",
                      "Lcluster3_Acluster1_pvalues",
                      "Lcluster4_Acluster1_pvalues",
                      "Lcluster5_Acluster1_pvalues",
                      "Lcluster6_Acluster1_pvalues",
                      "Lcluster7_Acluster1_pvalues",
                      "Lcluster8_Acluster1_pvalues",
                      "Lcluster9_Acluster1_pvalues",
                      "Lcluster10_Acluster1_pvalues",
                      "Lcluster0_Acluster2_pvalues",
                      "Lcluster1_Acluster2_pvalues",
                      "Lcluster2_Acluster2_pvalues",
                      "Lcluster3_Acluster2_pvalues",
                      "Lcluster4_Acluster2_pvalues",
                      "Lcluster5_Acluster2_pvalues",
                      "Lcluster6_Acluster2_pvalues",
                      "Lcluster7_Acluster2_pvalues",
                      "Lcluster8_Acluster2_pvalues",
                      "Lcluster9_Acluster2_pvalues",
                      "Lcluster10_Acluster2_pvalues")]
    colnames(pval) <- c("cluster0_Adipose",
                        "cluster1_Adipose",
                        "cluster2_Adipose",
                        "cluster3_Adipose",
                        "cluster4_Adipose",
                        "cluster5_Adipose",
                        "cluster6_Adipose",
                        "cluster7_Adipose",
                        "cluster8_Adipose",
                        "cluster9_Adipose",
                        "cluster10_Adipose",
                        "cluster0_Oenocyte",
                        "cluster1_Oenocyte",
                        "cluster2_Oenocyte",
                        "cluster3_Oenocyte",
                        "cluster4_Oenocyte",
                        "cluster5_Oenocyte",
                        "cluster6_Oenocyte",
                        "cluster7_Oenocyte",
                        "cluster8_Oenocyte",
                        "cluster9_Oenocyte",
                        "cluster10_Oenocyte",
                        "cluster0_Muscle",
                        "cluster1_Muscle",
                        "cluster2_Muscle",
                        "cluster3_Muscle",
                        "cluster4_Muscle",
                        "cluster5_Muscle",
                        "cluster6_Muscle",
                        "cluster7_Muscle",
                        "cluster8_Muscle",
                        "cluster9_Muscle",
                        "cluster10_Muscle")
    selected_rows = NULL
    selected_columns = NULL

    intr_pairs = data$interacting_pair
    all_pval = pval
    all_means = means

    if(is.null(selected_rows)){
      selected_rows = intr_pairs
    }

    if(is.null(selected_columns)){
      selected_columns = colnames(all_pval)
    }

    # sel_pval = all_pval[match(selected_rows, intr_pairs), selected_columns]
    # sel_means = all_means[match(selected_rows, intr_pairs), selected_columns]

    sel_pval = all_pval
    sel_means = all_means

    df_names = expand.grid(selected_rows, selected_columns)
    pval = unlist(sel_pval)
    pval[pval==0] = 0.0009
    head(pval)
    plot.data = cbind(df_names, pval)
    pr = unlist(sel_means)
    # pr[pr==0] = 0
    # pr[pr>0.1] = 0.1
    # plot.data = cbind(plot.data,log2(pr))
    plot.data = cbind(plot.data, pr)
    colnames(plot.data) = c("pair", "clusters", "pvalue", "mean")
    plot.data
  })


# Output ------------------------------------------------------------------

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

  # Export dotplot
  output$dotplot_output <- renderPlot({
    if (input$submitbutton>0) {
      isolate({
        my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)

        ggplot(dotplot_data(), aes(x=clusters, y=pair)) +
          geom_point(aes(size=-log10(pvalue), color=mean)) +
          scale_color_gradientn("score", colors=my_palette) +
          theme_bw() +
          theme(panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                axis.text=element_text(size=14, colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1),
                axis.text.y = element_text(size=12, colour = "black"),
                axis.title=element_blank(),
                panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
      })
    }
  })



} # server
