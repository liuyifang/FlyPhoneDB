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

  # heatmap data
  # mydata <- cor(mtcars)
  mydata <- matrix(runif(2500, min = -2, max = 2), ncol = 50)
  row.names(mydata) <- paste0("row_", seq_len(nrow(mydata)))
  colnames(mydata) <- paste0("col_", seq_len(ncol(mydata)))

  # dendrogram for rows
  hc <- hclust(dist(mydata), "ave")
  dhr <- as.dendrogram(hc)
  order_r <- rownames(mydata)[hc$order]

  # dendrogram for columns
  hc <- hclust(dist(t(mydata)), "ave")
  dhc <- as.dendrogram(hc)
  order_c <- colnames(mydata)[hc$order]

  # the data
  expr_set <- bind_cols(
    data_frame(rowvar = rownames(mydata)),
    as.data.frame(mydata)
  )
  expr_set <- gather(expr_set, colvar, measure, -rowvar)
  expr_set$rowvar <- factor( expr_set$rowvar, levels = order_r )
  expr_set$colvar <- factor( expr_set$colvar, levels = order_c )
  expr_set <- arrange(expr_set, rowvar, colvar)

  # get data for dendrograms - IMHO, ggdendro is the hero here...
  data_c <- dendro_data(dhc, type = "rectangle")
  data_c <- segment(data_c) %>% mutate(
    y = y + length(order_r) + .5,
    yend = yend + length(order_r) + .5
  )

  data_r <- dendro_data(dhr, type = "rectangle")
  data_r <- segment(data_r)
  data_r <- data_r %>%
    mutate( x_ = y + length(order_c) + .5,
            xend_ = yend + length(order_c) + .5,
            y_ = x,
            yend_ = xend )

  expr_set <- expr_set %>%
    mutate(
      tooltip = sprintf("Row: %s<br/>Col: %s<br/>measure: %.02f",
                        rowvar, colvar, measure) ,
      data_id = sprintf("%s_%s", rowvar, colvar)
    )

  # heatmap data2
  count_network <- read.csv("count_network.csv")
  str(count_network)
  count_network$tooltip <- paste0("source: ", count_network$SOURCE, "<br/>",
                                  "target: ", count_network$TARGET, "<br/>",
                                  "count: ", count_network$count, "<br/>")
  count_network$id <- paste0(count_network$SOURCE, ">", count_network$TARGET)

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
    plot.data$id <- row.names(plot.data)
    plot.data
  })


# Output ------------------------------------------------------------------

  # Export heatmap
  output$heatmap_girafe <- renderGirafe({
    if (input$submitbutton>0) {
      isolate({
        # all data are tidy and can be now used with ggplot
        gg_heatmap <- ggplot(data = expr_set, aes(x = colvar, y = rowvar) ) +
          geom_tile_interactive(aes(fill = measure, tooltip = tooltip, data_id = data_id), colour = "white") +
          scale_fill_gradient(low = "white", high = "#BC120A") +
          geom_segment(
            data = data_c,
            mapping = aes(x = x, y = yend, xend = xend, yend = y),
            colour = "gray20", size = .2) +
          geom_segment(
            data = data_r,
            mapping = aes(x = x_, y = y_, xend = xend_, yend = yend_),
            colour = "gray20", size = .2) +
          coord_equal()

        # cosmetics
        gg_heatmap <- gg_heatmap + theme_minimal() +
          theme(
            legend.position = "right",
            panel.grid.minor = element_line(color = "transparent"),
            panel.grid.major = element_line(color = "transparent"),
            axis.ticks.length   = unit(2, units = "mm"),
            plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
            axis.title = element_text(size = 9, colour = "gray30"),
            axis.text.y = element_text(hjust = 1, size = 5, colour = "gray40"),
            axis.text.x = element_text(angle = 90, hjust = 1, size = 5, colour = "gray40"),
            legend.title=element_text(face = "bold", hjust = 0.5, size=8),
            legend.text=element_text(size=6)
          )

        girafe(ggobj = gg_heatmap,
               options = list(opts_selection(type = "single")) )
      })
    }
  })

  # Export heatmap2
  output$heatmap2_girafe <- renderGirafe({
    if (input$submitbutton>0) {
      isolate({
        # all data are tidy and can be now used with ggplot
        gg_heatmap2 <- ggplot(data = count_network, aes(x = SOURCE, y = TARGET) ) +
          geom_tile_interactive(aes(fill = count, tooltip = tooltip, data_id = id), colour = "white") +
          scale_fill_gradient(low = "white", high = "#BC120A") +
          # geom_segment(
          #   data = data_c,
          #   mapping = aes(x = x, y = yend, xend = xend, yend = y),
          #   colour = "gray20", size = .2) +
          # geom_segment(
          #   data = data_r,
          #   mapping = aes(x = x_, y = y_, xend = xend_, yend = yend_),
          #   colour = "gray20", size = .2) +
          coord_equal()

        # cosmetics
        gg_heatmap2 <- gg_heatmap2 + theme_minimal() +
          theme(
            legend.position = "right",
            panel.grid.minor = element_line(color = "transparent"),
            panel.grid.major = element_line(color = "transparent"),
            axis.ticks.length   = unit(2, units = "mm"),
            plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
            axis.title = element_text(size = 9, colour = "gray30"),
            axis.text.y = element_text(hjust = 1, size = 5, colour = "gray40"),
            axis.text.x = element_text(angle = 90, hjust = 1, size = 5, colour = "gray40"),
            legend.title=element_text(face = "bold", hjust = 0.5, size=8),
            legend.text=element_text(size=6)
          )

        girafe(ggobj = gg_heatmap2,
               options = list(opts_selection(type = "single")) )
      })
    }
  })

  # Status/Output Text Box
  output$contents <- renderPrint({
    if (input$submitbutton>0) {
      isolate("Calculation complete.")
    } else {
      return("Server is ready for calculation.")
    }
  })

  # Export dotplot
  output$dotplot_girafe <- renderGirafe({
    if (input$submitbutton>0) {
      isolate({
        my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)

        gg_dot <- ggplot(dotplot_data(), aes(x=clusters, y=pair, tooltip = id, data_id = pair)) +
          geom_point_interactive(aes(size=-log10(pvalue), color=mean)) +
          scale_color_gradientn("score", colors=my_palette) +
          theme_bw() +
          theme(panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                axis.text=element_text(size=14, colour = "black"),
                axis.text.x = element_text(angle = 90, hjust = 1),
                axis.text.y = element_text(size=12, colour = "black"),
                axis.title=element_blank(),
                panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
        girafe(ggobj = gg_dot,
               options = list(opts_selection(type = "single")) )
      })
    }
  })

  output$choices <- renderPrint({
    input$dotplot_girafe_selected
  })

  # Export UMAP_plot
  output$UMAP_plot <- renderPlot({


    UMAP <- read.csv("2020-10-08_MT_UMAP.csv", row.names = 1)
    matrix <- readRDS("2020-10-08_MT_matrix.Rds")
    UMAP <- UMAP[colnames(matrix), ]
    # plot(UMAP, pch=16)

    id <- "128up"
    l <- apply(as.data.frame(matrix[id, ]) - .1, 1, sum) + .1
    f <- l == 0
    l <- log2(l)
    l[f] <- NA
    mi <- min(l,na.rm=TRUE)
    ma <- max(l,na.rm=TRUE)
    ColorRamp <- colorRampPalette(c("#F2F2F2", "#F1E9E9", "#EEDBDB", "#F96060", "#FF0000"))(256)
    ColorLevels <- seq(mi, ma, length=length(ColorRamp))
    v <- round((l - mi)/(ma - mi)*255 + 1,0)
    par(mfrow=c(1,1))
    layout(matrix(data=c(1,3,2,4), nrow=2, ncol=2), widths=c(5,1,5,1), heights=c(5,1,1,1))
    plot(UMAP, main=input$dotplot_girafe_selected,pch=20,cex=0.3,col="#F2F2F2",frame.plot=F, axes=FALSE, xlab="", ylab="")
    rand_v <- 1:length(v)
    rand_v <- sample(rand_v)
    for ( k in rand_v ){
      points(UMAP[k,1],UMAP[k,2],col=ColorRamp[v[k]],pch=20,cex=0.3)
    }
  })


} # server
