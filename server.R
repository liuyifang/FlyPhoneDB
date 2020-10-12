# Import libraries
library(shiny)
library(shinythemes)
library(ggplot2)

options(stringsAsFactors = FALSE)

####################################
# Server                           #
####################################

# Define server function
server <- function(input, output) {

  # heatmap data2
  count_network <- read.csv("data/MT/2020-10-10_count_network_MT.csv")

  data_wide <- pivot_wider(count_network, names_from = TARGET, values_from = count) %>% as.data.frame()
  str(data_wide)
  row.names(data_wide) <- data_wide$SOURCE
  data_wide$SOURCE <- NULL

  # dendrogram for rows
  hc <- hclust(dist(data_wide), "ave")
  dhr <- as.dendrogram(hc)
  order_r <- rownames(data_wide)[hc$order]

  # dendrogram for columns
  hc <- hclust(dist(t(data_wide)), "ave")
  dhc <- as.dendrogram(hc)
  order_c <- colnames(data_wide)[hc$order]

  # # get data for dendrograms - IMHO, ggdendro is the hero here...
  data_c <- dendro_data(dhc, type = "rectangle")
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  # data_c <- segment(data_c) %>% mutate(
  #   y = range01(y) + length(order_r) + .5,
  #   yend = range01(yend) + length(order_r) + .5
  # )
  data_c <- segment(data_c)
  data_c_y <- data_c[ , c("y", "yend")]
  data_c_y01 <- range01(data_c_y)
  data_c$y <- data_c_y01$y + length(order_r) + .5
  data_c$yend <- data_c_y01$yend + length(order_r) + .5

  data_r <- dendro_data(dhr, type = "rectangle")
  data_r <- segment(data_r)

  data_r_y <- data_r[ , c("y", "yend")]
  data_r_y01 <- range01(data_r_y)
  data_r$y <- data_r_y01$y + length(order_c) + .5
  data_r$yend <- data_r_y01$yend + length(order_c) + .5

  data_r <- data_r %>%
    mutate( x_ = y,
            xend_ = yend,
            y_ = x,
            yend_ = xend )

  # str(count_network)
  count_network$tooltip <- paste0("source: ", count_network$SOURCE, "<br/>",
                                  "target: ", count_network$TARGET, "<br/>",
                                  "count: ", count_network$count, "<br/>")
  count_network$id <- paste0(count_network$SOURCE, ">", count_network$TARGET)


# Output ------------------------------------------------------------------

  # Export heatmap2
  output$heatmap2_girafe <- renderGirafe({
        # all data are tidy and can be now used with ggplot
        if(input$heatmap_correlation){
          gg_heatmap2 <- ggplot(data = count_network, aes(x = SOURCE, y = TARGET) ) +
            geom_tile_interactive(aes(fill = count, tooltip = tooltip, data_id = id), colour = "white") +
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
        }else{
          gg_heatmap2 <- ggplot(data = count_network, aes(x = SOURCE, y = TARGET) ) +
            geom_tile_interactive(aes(fill = count, tooltip = tooltip, data_id = id), colour = "white") +
            scale_fill_gradient(low = "white", high = "#BC120A") +
            coord_equal()
        }

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
               options = list(opts_selection(type = "multiple")) )
  })

  output$heatmap2_choices <- renderPrint({
    input$heatmap2_girafe_selected
  })

  # Export dotplot
  output$dotplot_girafe <- renderGirafe({
    my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=255)

    dotplot_data <- read.csv("data/MT/2020-10-10_MT_dotplot_data.csv",
                             row.names = 1)

    if (is.null(input$heatmap2_girafe_selected)) {
      print(paste0("input$heatmap2_girafe_selected: is null"))

      dotplot_data_order_by_score <- dotplot_data[order(-dotplot_data$score), ]
      pair_top <- head(unique(dotplot_data_order_by_score$pair), 4)
      clusters_top <- head(unique(dotplot_data_order_by_score$clusters), 10)

      dotplot_data_order_by_score_top <- subset(dotplot_data, pair %in% pair_top & clusters %in% clusters_top)

      gg_dot <- ggplot(dotplot_data_order_by_score_top, aes(x=clusters, y=pair, tooltip = id, data_id = pair)) +
        geom_point_interactive(aes(size=-log10(pvalue), color=score)) +
        scale_color_gradientn("score", colors=my_palette) +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              axis.text=element_text(size=14, colour = "black"),
              axis.text.x = element_text(angle = 45, hjust = 1),
              axis.text.y = element_text(size=12, colour = "black"),
              axis.title=element_blank())

    }else{
      dotplot_data_selected <- subset(dotplot_data, clusters %in% input$heatmap2_girafe_selected)

      dotplot_data_order_by_score <- dotplot_data_selected[order(-dotplot_data_selected$score), ]
      pair_top <- head(unique(dotplot_data_order_by_score$pair), 5)
      dotplot_data_selected_order_by_score_top <- subset(dotplot_data_selected, pair %in% pair_top)

      gg_dot <- ggplot(dotplot_data_selected_order_by_score_top, aes(x=clusters, y=pair, tooltip = id, data_id = pair)) +
        geom_point_interactive(aes(size=-log10(pvalue), color=score)) +
        scale_color_gradientn("score", colors=my_palette) +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              axis.text=element_text(size=14, colour = "black"),
              axis.text.x = element_text(angle = 45, hjust = 1),
              axis.text.y = element_text(size=12, colour = "black"),
              axis.title=element_blank())

    }

    girafe(ggobj = gg_dot,
           options = list(opts_selection(type = "single")) )

  })

  output$dotplot_table <- renderDataTable({
    if (is.null(input$heatmap2_girafe_selected)) {
      dotplot_data
    }else{
      dotplot_data_selected <- subset(dotplot_data, clusters %in% input$heatmap2_girafe_selected)
      dotplot_data_selected
    }
  })

  output$dotplot_choices <- renderPrint({
    input$dotplot_girafe_selected
  })

  # Export UMAP_plot
  output$UMAP_plot <- renderPlot({


    UMAP <- read.csv("data/MT/2020-10-08_MT_UMAP.csv", row.names = 1)
    # str(UMAP)
    matrix <- readRDS("data/MT/2020-10-10_exprMat_norm_filter.Rds")
    UMAP <- UMAP[colnames(matrix), ]
    # plot(UMAP, pch=16)

    id <- "Sema2b"
    l <- apply(as.data.frame(matrix[id, ]) - .1, 1, sum) + .1
    f <- l == 0
    # l <- log2(l)
    l[f] <- NA
    mi <- min(l,na.rm=TRUE)
    ma <- max(l,na.rm=TRUE)
    ColorRamp <- colorRampPalette(c("#F2F2F2", "#F1E9E9", "#EEDBDB", "#F96060", "#FF0000"))(256)
    ColorLevels <- seq(mi, ma, length=length(ColorRamp))
    v <- round((l - mi)/(ma - mi)*255 + 1,0)
    par(mfrow=c(1,1))
    # layout(matrix(data=c(1,3,2,4), nrow=2, ncol=2), widths=c(5,1,5,1), heights=c(5,1,1,1))
    plot(UMAP, main=input$dotplot_girafe_selected,pch=20,cex=0.3,col="#F2F2F2",frame.plot=F, axes=FALSE, xlab="", ylab="")
    rand_v <- 1:length(v)
    rand_v <- sample(rand_v)
    for ( k in rand_v ){
      points(UMAP[k,1],UMAP[k,2],col=ColorRamp[v[k]],pch=20,cex=0.3)
    }
  })


} # server
