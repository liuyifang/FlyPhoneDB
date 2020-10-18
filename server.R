# Import libraries
# library(shiny)
# library(shinythemes)
# library(ggplot2)

options(stringsAsFactors = FALSE)

####################################
# Server                           #
####################################

# Define server function
server <- function(input, output) {

  # read heatmap data2
  count_network <- read.csv("data/MT/2020-10-10_count_network_MT.csv")

  data_wide <- pivot_wider(count_network, names_from = TARGET, values_from = count) %>% as.data.frame()
  # str(data_wide)
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

    # read
    dotplot_data <- read.csv("data/MT/2020-10-12_dotplot_data_MT.csv",
                             row.names = 1)

    if (is.null(input$heatmap2_girafe_selected)) {
      # print(paste0("input$heatmap2_girafe_selected: is null"))

      dotplot_data_order_by_score <- dotplot_data[order(-dotplot_data$score), ]
      pair_top <- head(unique(dotplot_data_order_by_score$pair), 4)
      clusters_top <- head(unique(dotplot_data_order_by_score$clusters), 10)

      dotplot_data_order_by_score_top <- subset(dotplot_data, pair %in% pair_top & clusters %in% clusters_top)

      gg_dot <- ggplot(dotplot_data_order_by_score_top, aes(x=clusters, y=pair, tooltip = id, data_id = id)) +
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

      gg_dot <- ggplot(dotplot_data_selected_order_by_score_top, aes(x=clusters, y=pair, tooltip = id, data_id = id)) +
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
      # read
      dotplot_data_table <- read.csv("data/MT/2020-10-12_dotplot_data_MT.csv",
                               row.names = 1)
      dotplot_data_table$id <- NULL
      dotplot_data_table
    }else{
      # read
      dotplot_data_table <- read.csv("data/MT/2020-10-12_dotplot_data_MT.csv",
                                     row.names = 1)
      dotplot_data_table$id <- NULL
      dotplot_data_selected_table <- subset(dotplot_data_table, clusters %in% input$heatmap2_girafe_selected)
      dotplot_data_selected_table
    }
  })

  # # UMAP cluster
  # output$UMAP_cluster_girafe <- renderGirafe({
  #
  #   UMAP_cluster_annotation <- readRDS("data/MT/2020-10-12_UMAP_matrix_down_sample_to_3000.Rds")
  #
  #   if (is.null(input$dotplot_girafe_selected)) {
  #
  #     gg_UMAP_cluster <- ggplot(UMAP_cluster_annotation, aes(x = UMAP_1, y = UMAP_2,
  #                           colour = annotation, group = annotation)) +
  #       geom_point_interactive(aes(tooltip = annotation, data_id = annotation)) +
  #       scale_color_viridis_d() +
  #       theme_minimal() +
  #       theme(legend.position = "none")
  #
  #   }else{
  #     # print(input$dotplot_girafe_selected)
  #     split.var <- strsplit(input$dotplot_girafe_selected, "\\|")[[1]][1]
  #     split.var <- strsplit(split.var, ">")[[1]]
  #     # print(split.var)
  #     UMAP_cluster_annotation_filter <- subset(UMAP_cluster_annotation, annotation %in% split.var)
  #     UMAP_1_min <- min(UMAP_cluster_annotation$UMAP_1)
  #     UMAP_1_max <- max(UMAP_cluster_annotation$UMAP_1)
  #     UMAP_2_min <- min(UMAP_cluster_annotation$UMAP_2)
  #     UMAP_2_max <- max(UMAP_cluster_annotation$UMAP_2)
  #
  #     gg_UMAP_cluster <- ggplot(UMAP_cluster_annotation_filter, aes(x = UMAP_1, y = UMAP_2,
  #                                                            colour = annotation, group = annotation)) +
  #       geom_point_interactive(aes(tooltip = annotation, data_id = annotation)) +
  #       scale_color_viridis_d() +
  #       theme_minimal() +
  #       theme(legend.position = "none") +
  #       xlim(UMAP_1_min, UMAP_1_max) +
  #       ylim(UMAP_2_min, UMAP_2_max)
  #   }
  #
  #   girafe(ggobj = gg_UMAP_cluster,
  #          options = list(
  #            opts_hover_inv(css = "opacity:0.1;"),
  #            opts_hover(css = "stroke-width:2;")
  #          )
  #   )
  #
  # })
  #
  # # UMAP gene
  # output$UMAP_plot <- renderGirafe({
  #
  #
  #   # UMAP <- read.csv("data/MT/2020-10-08_MT_UMAP.csv", row.names = 1)
  #   # print(UMAP)
  #   matrix <- readRDS("data/MT/2020-10-12_exprMat_norm_filter.Rds") %>% t()
  #   UMAP <- read.csv("data/MT/2020-10-08_MT_UMAP.csv", row.names = 1)
  #   UMAP <- UMAP[row.names(matrix), ]
  #   UMAP_matrix <- cbind(UMAP, matrix)
  #   # plot(UMAP, pch=16)
  #
  #   # print(input$dotplot_girafe_selected)
  #   # split.var <- strsplit("apolpp_dally", "_")
  #   # class(split.var)
  #   # print(split.var[[1]][1])
  #   # split.var <- c("Sema2b" , "dally" )
  #   # print(UMAP_matrix[ , "Sema2b"])
  #
  #   if (is.null(input$dotplot_girafe_selected)) {
  #     # change to total UMI
  #     gg_UMAP_gene <- ggplot(UMAP_matrix, aes(x = UMAP_1, y = UMAP_2,
  #                             colour = log10UMI)) + # https://stackoverflow.com/questions/22309285/how-to-use-a-variable-to-specify-column-name-in-ggplot
  #       geom_point() +
  #       # scale_color_viridis_d() +
  #       theme_minimal() +
  #       theme(legend.position = "none")
  #   }else{
  #     # print(input$dotplot_girafe_selected)
  #     # strsplit("hemocyte>? Crop or hindgut ?|Hml_dally", "\\|")[[1]][2]
  #     split.var <- strsplit(input$dotplot_girafe_selected, "\\|")[[1]][2]
  #
  #     # print(split.var)
  #     # strsplit("Hml_dally", "_")[[1]]
  #     split.var <- strsplit(split.var, "_")[[1]]
  #
  #     gg_UMAP_gene <- ggplot(UMAP_matrix, aes(x = UMAP_1, y = UMAP_2,
  #                             colour = get(split.var[1]))) + # https://stackoverflow.com/questions/22309285/how-to-use-a-variable-to-specify-column-name-in-ggplot
  #       geom_point() +
  #       # scale_color_viridis_d() +
  #       theme_minimal() +
  #       theme(legend.position = "none")
  #   }
  #
  #   girafe(ggobj = gg_UMAP_gene)
  #
  # })


} # server
