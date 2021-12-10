start_time <- Sys.time()
set.seed(123)

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(future.apply)
  library(Seurat)
  library(RColorBrewer)
  library(reshape2)
  library(network)
  library(igraph)
})

option_list = list(
  make_option(c("-i", "--matrix"), type="character", default=NULL,
              help="input matrix", metavar="character"),
  make_option(c("-a", "--metadata"), type="character", default=NULL,
              help="input metadata", metavar="character"),
  make_option(c("-p", "--lrpair"), type="character", default=NULL,
              help="annotation ligand receptor", metavar="character"),
  make_option(c("-s", "--corecomponents"), type="character", default=NULL,
              help="annotation core components", metavar="character"),
  make_option(c("-c", "--cores"), type="character", default=NULL,
              help="number of cores", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="output directory name", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(paste0("matrix input: ", opt$matrix))
print(paste0("metadata input: ", opt$metadata))
print(paste0("L-R pairs input: ", opt$lrpair))
print(paste0("core components input: ", opt$corecomponents))
print(paste0("number of cores: ", opt$cores))
print(paste0("output directory: ", opt$output))

output_dir <- opt$output
if (!dir.exists(output_dir)) {dir.create(output_dir)}

if (!dir.exists(paste0(output_dir, "/heatmap"))) {dir.create(paste0(output_dir, "/heatmap"), recursive = TRUE)}
if (!dir.exists(paste0(output_dir, "/dotplot"))) {dir.create(paste0(output_dir, "/dotplot"), recursive = TRUE)}
if (!dir.exists(paste0(output_dir, "/circleplot"))) {dir.create(paste0(output_dir, "/circleplot"), recursive = TRUE)}

# plan(multiprocess, workers = 8) ## => parallelize on your local computer
plan(multiprocess, workers = as.numeric(opt$cores)) ## => parallelize on your local computer

####################################
# Input and cluster means
####################################

# seuratObj <- readRDS("../2021-02-15_FlyPhone_midgut/Data/2021-02-15_midgut_seuratObj.Rds")

# metadata <- seuratObj@meta.data
# Idents(seuratObj) <- "celltype"
# Idents(seuratObj)
#
# exprMat_original <- GetAssayData(seuratObj, assay = "RNA", slot = "counts")
# exprMat_original[1:3, 1:3]

# matrix <- as.matrix(exprMat_original)
# write.csv(matrix, file = "2021-02-25_matrix_midgut.csv")

# exprMat <- log1p(sweep(exprMat_original, 2, Matrix::colSums(exprMat_original), FUN = "/") * 10000)

# exprMat <- read.csv("./test_dataset_60M/2021-03-18_matrix_midgut.csv", row.names = 1, check.names = FALSE)
exprMat <- read.csv(opt$matrix, row.names = 1, check.names = FALSE)
# cellInfo <- read.csv("./test_dataset_60M/2021-03-18_metadata_midgut.csv", row.names = 1)
cellInfo <- read.csv(opt$metadata, row.names = 1)
cellInfo$celltype <- as.character(cellInfo$celltype)
# str(cellInfo)

exprMat <- exprMat[ , row.names(cellInfo)]

# create seuratObj
seuratObj <- CreateSeuratObject(counts = exprMat)
seuratObj <- NormalizeData(seuratObj)
seuratObj <- FindVariableFeatures(seuratObj, selection.method = "vst", nfeatures = 2000)
all_genes <- rownames(seuratObj)
seuratObj <- ScaleData(seuratObj, features = all_genes)
seuratObj$celltype <- as.factor(cellInfo$celltype)
Idents(seuratObj) <- "celltype"
# clusterMetadataTable <- table(seuratObj@meta.data[ , "celltype"]) %>% as.data.frame()
# colnames(clusterMetadataTable) <- c("celltype", "count")
# print(clusterMetadataTable$celltype)

exprMat <- sweep(exprMat, 2, Matrix::colSums(exprMat), FUN = "/") * 10000

# cellInfo <- metadata[ , c("celltype"), drop = FALSE]
# row.names(cellInfo) <- cellInfo$barcode
# colnames(cellInfo) <- c("ID", "res", "celltype")
# str(cellInfo)

# cellInfo_test <- subset(cellInfo, celltype %in% c("aEC1", "aEC2", "aEC3"))
# exprMat_test <- exprMat[ , row.names(cellInfo_test)]
#
# write.csv(cellInfo_test, file = "./test_dataset/2021-02-25_cellInfo_midgut_test.csv")
# write.csv(exprMat_test, file = "./test_dataset/2021-02-25_matrix_midgut_test.csv")

# write.csv(cellInfo, file = "2021-02-25_cellInfo_midgut.csv")

# saveRDS(exprMat, file = "../Data/2021-02-15_exprMat_midgut.Rds")
# saveRDS(cellInfo, file = "../Data/2021-02-15_cellInfo_midgut.Rds")

# exprMat <- readRDS("../Data/2021-02-15_exprMat_midgut.Rds")
# cellInfo <- readRDS("../Data/2021-02-15_cellInfo_midgut.Rds")

# LR_pairs <- read.csv(file = "./annotation/Ligand_receptor_pair_high_confident_2021vs1_clean.txt", sep = "\t")
LR_pairs <- read.csv(file = opt$lrpair, sep = "\t")
# str(LR_pairs)
# LR_pairs <- interaction_input[ , c("Gene_secreted", "Gene_receptor", "pathway_receptor")]
# colnames(LR_pairs) <- c("Gene_secreted", "Gene_receptor", "pathway_receptor")

# start real score
gene_list <- unique(c(LR_pairs$Gene_secreted, LR_pairs$Gene_receptor))
common_genes <- intersect(gene_list, row.names(exprMat))

LR_pairs <- subset(LR_pairs, Gene_secreted %in% common_genes & Gene_receptor %in% common_genes)

exprMat <- as.matrix(exprMat)
exprMat <- t(exprMat)
df_Ligand <- exprMat[ , unique(LR_pairs$Gene_secreted)]
df_Receptor <- exprMat[ , unique(LR_pairs$Gene_receptor)]

celltype_df_Ligand <- cbind(cellInfo[ , c("celltype"), drop = FALSE], df_Ligand)
celltype_df_Receptor <- cbind(cellInfo[ , c("celltype"), drop = FALSE], df_Receptor)
# celltype_df_Ligand[1:3, 1:3]
# celltype_df_Receptor[1:3, 1:3]

# average Ligand counts by each celltype
df_group_by_celltype_Ligand <- celltype_df_Ligand %>%
  group_by(celltype) %>%
  summarise_all(mean) %>%
  as.data.frame()
# df_group_by_celltype_Ligand[1:3, 1:3]
row.names(df_group_by_celltype_Ligand) <- df_group_by_celltype_Ligand$celltype
df_group_by_celltype_Ligand$celltype <- NULL
df_group_by_celltype_Ligand <- t(df_group_by_celltype_Ligand)
# df_group_by_celltype_Ligand[1:3, 1:3]
# str(df_group_by_celltype_Ligand)
# write.csv(df_group_by_celltype_Ligand,
#           file = "../Data/2021-02-15_df_group_by_celltype_Ligand.csv")

# average Receptor counts by each celltype
df_group_by_celltype_Receptor <- celltype_df_Receptor %>%
  group_by(celltype) %>%
  summarise_all(mean) %>%
  as.data.frame()
# df_group_by_celltype_Receptor[1:3, 1:3]
row.names(df_group_by_celltype_Receptor) <- df_group_by_celltype_Receptor$celltype
df_group_by_celltype_Receptor$celltype <- NULL
df_group_by_celltype_Receptor <- t(df_group_by_celltype_Receptor)
# df_group_by_celltype_Receptor[1:3, 1:3]
# str(df_group_by_celltype_Receptor)
# write.csv(df_group_by_celltype_Receptor,
#           file = "../Data/2021-02-15_df_group_by_celltype_Receptor.csv")

####################################
# Interaction score
####################################

ligand_avg <- df_group_by_celltype_Ligand[LR_pairs$Gene_secreted, ] %>% as.data.frame()
receptor_avg <- df_group_by_celltype_Receptor[LR_pairs$Gene_receptor, ] %>% as.data.frame()
# write.csv(ligand_avg, file = "../Data/2021-02-15_ligand_avg.csv")
# write.csv(receptor_avg, file = "../Data/2021-02-15_receptor_avg.csv")

# x <- sort(unique(cellInfo$celltype)) %>% as.data.frame()

interaction_list <- list()
LR_pairs_one <- LR_pairs # combine

for (i in sort(unique(cellInfo$celltype)) ) {
# for (i in c("aEC1", "aEC2")) {
  # print("")
  # print(i)
  # print("")
  LR_pairs_combine <- LR_pairs # combine
  for (j in sort(unique(cellInfo$celltype)) ) {
  # for (j in c("aEC3", "aEC4") ) {
    print(paste0(i, ">", j))
    LR_pairs_tmp <- LR_pairs
    LR_pairs_tmp[[paste0(i, ">", j, "_score")]] <- log1p(ligand_avg[[i]]) * log1p(receptor_avg[[j]])

    # permutation -------------------------------------------------------------
    # start permutatioin
    permutation_times <- 1000
    y <- future_lapply(1:permutation_times, function(ii) {
    # for (ii in 1:permutation_times) {
    #   LR_pairs_tmp[[paste0("permute", ii)]] <- local({
      # print(i)
      # sample Ligand
      cellInfo_sample_Ligand <- cellInfo
      # str(cellInfo_sample_Ligand)
      cellInfo_sample_Ligand$celltype <- sample(cellInfo_sample_Ligand$celltype)

      celltype_df_sample_Ligand <- cbind(cellInfo_sample_Ligand[ , c("celltype"), drop = FALSE], df_Ligand)

      df_group_by_celltype_sample_Ligand <- celltype_df_sample_Ligand %>%
        group_by(celltype) %>%
        summarise_all(mean) %>%
        as.data.frame()
      # str(df_group_by_celltype_sample_Ligand)
      row.names(df_group_by_celltype_sample_Ligand) <- df_group_by_celltype_sample_Ligand$celltype
      df_group_by_celltype_sample_Ligand$celltype <- NULL
      df_group_by_celltype_sample_Ligand <- t(df_group_by_celltype_sample_Ligand)
      # str(df_group_by_celltype_sample_Ligand)

      # sample Receptor
      cellInfo_sample_Receptor <- cellInfo
      cellInfo_sample_Receptor$celltype <- sample(cellInfo_sample_Receptor$celltype)

      celltype_df_sample_Receptor <- cbind(cellInfo_sample_Receptor[ , c("celltype"), drop = FALSE], df_Receptor)

      df_group_by_celltype_sample_Receptor <- celltype_df_sample_Receptor %>%
        group_by(celltype) %>%
        summarise_all(mean) %>%
        as.data.frame()
      # str(df_group_by_celltype_sample_Receptor)
      row.names(df_group_by_celltype_sample_Receptor) <- df_group_by_celltype_sample_Receptor$celltype
      df_group_by_celltype_sample_Receptor$celltype <- NULL
      df_group_by_celltype_sample_Receptor <- t(df_group_by_celltype_sample_Receptor)
      # df_group_by_celltype_sample_Receptor[1:2, 1:2]
      # str(df_group_by_celltype_sample_Receptor)

      ####################################
      # Interaction score
      ####################################
      ligand_avg_tmp <- df_group_by_celltype_sample_Ligand[LR_pairs$Gene_secreted, ] %>% as.data.frame()
      receptor_avg_tmp <- df_group_by_celltype_sample_Receptor[LR_pairs$Gene_receptor, ] %>% as.data.frame()
      # LR_pairs_tmp <- LR_pairs
      # score <- ligand_avg[[paste0("Ligand_cluster", i)]] * receptor_avg[[paste0("Receptor_cluster", j)]]
      # colnames(LR_pairs_tmp)[ncol(LR_pairs_tmp)] <- paste0("permute", i)
      tmp <- log1p(ligand_avg_tmp[[i]]) * log1p(receptor_avg_tmp[[j]])
      tmp
    }, future.seed = TRUE)
    # }

    df <- data.frame(matrix(unlist(y), nrow=length(y), byrow=TRUE))
    df <- t(df)
    LR_pairs_tmp <- cbind(LR_pairs_tmp, df)
    # dim(LR_pairs_tmp)
    # dim(df)

    LR_pairs_tmp$result <- rowSums(sapply(LR_pairs_tmp[, 13:ncol(LR_pairs_tmp)], function(x) x > LR_pairs_tmp[[paste0(i, ">", j, "_score")]]))
    # head(LR_pairs_tmp[ , c("PM1_nonhemo_interaction_score", "result")])
    LR_pairs_tmp[[paste0(i, ">", j, "_pvalues")]] <- LR_pairs_tmp$result / permutation_times
    # write.csv(LR_pairs_tmp, file = "LR_pairs_tmp.csv")

    LR_pairs_tmp <- LR_pairs_tmp[ , c(1:12, ncol(LR_pairs_tmp))]
    # LR_pairs_tmp[LR_pairs_tmp$PM1_nonhemo_interaction_score]
    LR_pairs_tmp[LR_pairs_tmp[[paste0(i, ">", j, "_score")]] == 0, paste0(i, ">", j, "_pvalues")] <- 1

    LR_pairs_combine <- cbind(LR_pairs_combine,
                       LR_pairs_tmp[ , c(paste0(i, ">", j, "_score"), paste0(i, ">", j, "_pvalues"))]
                       )
    LR_pairs_one <- cbind(LR_pairs_one,
                              LR_pairs_tmp[ , c(paste0(i, ">", j, "_score"), paste0(i, ">", j, "_pvalues"))]
    )
  }

  interaction_list[[i]] <- LR_pairs_combine

}

write.csv(LR_pairs_one, file = paste0(output_dir, "/", "interaction_list.csv"))

end_time <- Sys.time()

end_time - start_time

# heatmap
avgexp <- AverageExpression(seuratObj, assay = "RNA", return.seurat = TRUE)

# Pathway_core_components <- read.table("./annotation/Pathway_core_components_2021vs1_clean.txt", sep = "\t", header = TRUE)
Pathway_core_components <- read.table(file = opt$corecomponents, sep = "\t", header = TRUE)

for(i in unique(Pathway_core_components$pathway) ){
  cat("\n")
  cat("## ", i, " {.tabset} \n")
  df <- subset(Pathway_core_components, pathway == i)
  genes <- df$gene

  cat("\n")
  p <- DoHeatmap(avgexp, features = genes, label = TRUE ,draw.lines = FALSE, raster = FALSE, angle = 90) +
    scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n =4, name = "RdBu"))) # & NoLegend()
  print(p)
  cat("\n")
  ggsave(p, file = paste0(output_dir, "/heatmap/heatmap_", i, ".png"),   # The directory you want to save the file in
         width = 8, # The width of the plot in inches
         height = 12)

}

# dotplot
data_original <- LR_pairs_one

data_original$interacting_pair <- paste(data_original$Gene_secreted, data_original$Gene_receptor, sep = "_")

pathways <- unique(data_original$pathway_receptor)
pathways <- pathways[-length(pathways)]

celltypes <- sort(unique(cellInfo$celltype))

for (p in pathways) {
  for (c in celltypes) {
    data_pathway <- subset(data_original, pathway_receptor == p)

    data <- data_pathway[ , grepl(paste0(c, ">"), colnames(data_pathway)), drop = FALSE]
    data$interacting_pair <- data_pathway$interacting_pair

    score <- data[ , grepl("_score", colnames(data)), drop = FALSE]
    colnames(score) <- str_replace(colnames(score), "_score", "")

    # 下次 input 改成 pvalue
    pvalue <- data[ , grepl("_pvalues", colnames(data)), drop = FALSE]
    colnames(pvalue) <- str_replace(colnames(pvalue), "_pvalues", "")

    selected_rows = NULL
    selected_columns = NULL

    intr_pairs = data$interacting_pair
    all_pvalue = pvalue
    all_score = score

    if(is.null(selected_rows)){
      selected_rows = intr_pairs
    }

    if(is.null(selected_columns)){
      selected_columns = colnames(all_pvalue)
    }

    sel_pvalue = all_pvalue
    sel_score = all_score

    df_names = expand.grid(selected_rows, selected_columns)
    pvalue = unlist(sel_pvalue)
    pvalue[pvalue == 0] <- 0.0009
    head(pvalue)
    plot.data = cbind(df_names, pvalue)
    pr = unlist(sel_score)
    # pr[pr==0] = 0
    # pr[pr>0.1] = 0.1
    # plot.data = cbind(plot.data,log2(pr))
    plot.data = cbind(plot.data, pr)
    colnames(plot.data) = c("pair", "clusters", "pvalue", "score")
    plot.data$id <- paste(plot.data$clusters, plot.data$pair, sep = "|")
    # plot.data

    # write.csv(plot.data, file = "2021-01-25_dotplot_data_Abdomen.csv")

    # my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)
    # my_palette <- wes_palette("Zissou1", 10, type = "continuous")
    # my_palette <- c("#A6A6A6", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B")
    # my_palette <- c("lightgrey", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B")
    my_palette <- colorRampPalette(brewer.pal(9, "Blues"))(100)
    my_palette_white <- rep("white", 100)

    dotplot_data <- plot.data

    # print(paste0("input$heatmap2_girafe_selected: is null"))

    # dotplot_data <- subset(dotplot_data, clusters == "main segment stellate cell>main segment PC")

    # The height of the plot in inches

    if(sum(dotplot_data$score) == 0){
      temp_plot <- ggplot(dotplot_data, aes(x=clusters, y=pair)) +
        geom_point(aes(size=-log10(pvalue), color=score)) +
        scale_color_gradientn("score", colors=my_palette_white) +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              axis.text=element_text(size=14, colour = "black"),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              axis.text.y = element_text(size=10, colour = "black"),
              axis.title=element_blank())
    }else{
      temp_plot <- ggplot(dotplot_data, aes(x=clusters, y=pair)) +
        geom_point(aes(size=-log10(pvalue), color=score)) +
        scale_color_gradientn("score", colors=my_palette) +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              axis.text=element_text(size=14, colour = "black"),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
              axis.text.y = element_text(size=10, colour = "black"),
              axis.title=element_blank())
    }

    ggsave(temp_plot, file = paste0(output_dir, "/dotplot/dotplot_", c, "_", p, ".png"),   # The directory you want to save the file in
           width = 8, # The width of the plot in inches
           height = 9)

  }
}

####
# circle plot
####

# LR_pairs_one <- read.csv("../../2021-05-17_FlyPhone_Blood_injured/Data/2021-05-17_LR_pairs_one.csv", row.names = 1, check.names = FALSE)
# LR_pairs_one <- LR_pairs_one[ , -grep("18", colnames(LR_pairs_one))]

pvalues <- colnames(LR_pairs_one)[grepl("_pvalues", colnames(LR_pairs_one))]
scores <- colnames(LR_pairs_one)[grepl("_score", colnames(LR_pairs_one))]
# celltypes <- c("0", "2", "4", "5", "6", "10", "11", "8", "9", "13", "16", "17", "15", "7", "3", "12")
celltypes <- sort(unique(cellInfo$celltype))

interaction_pvalues <- LR_pairs_one[ , c("Gene_secreted", "Gene_receptor", "pathway_receptor", pvalues)]
interaction_scores <- LR_pairs_one[ , c("Gene_secreted", "Gene_receptor", "pathway_receptor", scores)]

# pathway <- "EGFR signaling pathway"

pathways <- c("EGFR signaling pathway", "PVR RTK signaling pathway", "FGFR signaling pathway", "HEDGEHOG signaling pathway", "HIPPO signaling pathway", "INSULIN signaling pathway", "JAK-STAT signaling pathway", "NOTCH signaling pathway", "TGF beta signaling pathway", "TNF alpha signaling pathway", "WNT signaling pathway", "Toll signaling pathway", "Torso signaling pathway")

for (pathway in pathways) {
  print(pathway)
  interaction_pathway_pvalues <- subset(interaction_pvalues, pathway_receptor == pathway)
  interaction_pathway_scores <- subset(interaction_scores, pathway_receptor == pathway)

  interaction_pathway_pvalues$pathway_receptor <- NULL
  interaction_pathway_long_pvalues <- melt(interaction_pathway_pvalues, id.vars = c("Gene_secreted", "Gene_receptor"))
  interaction_pathway_long_pvalues$variable <- str_replace(interaction_pathway_long_pvalues$variable, "_pvalues", "")
  colnames(interaction_pathway_long_pvalues)[4] <- "pvalue"
  row.names(interaction_pathway_long_pvalues) <- paste(interaction_pathway_long_pvalues$Gene_secreted, interaction_pathway_long_pvalues$Gene_receptor, interaction_pathway_long_pvalues$variable, sep="_")

  interaction_pathway_scores$pathway_receptor <- NULL
  interaction_pathway_long_scores <- melt(interaction_pathway_scores, id.vars = c("Gene_secreted", "Gene_receptor"))
  interaction_pathway_long_scores$variable <- str_replace(interaction_pathway_long_scores$variable, "_score", "")
  colnames(interaction_pathway_long_scores)[4] <- "score"
  row.names(interaction_pathway_long_scores) <- paste(interaction_pathway_long_scores$Gene_secreted, interaction_pathway_long_scores$Gene_receptor, interaction_pathway_long_scores$variable, sep="_")

  all.equal(row.names(interaction_pathway_long_pvalues), row.names(interaction_pathway_long_scores))
  interaction_pathway_long_pvalues$Gene_secreted <- NULL
  interaction_pathway_long_pvalues$Gene_receptor <- NULL
  interaction_pathway_long_pvalues$variable <- NULL
  interaction_pathway_long <- cbind(interaction_pathway_long_scores, interaction_pathway_long_pvalues)

  # interaction_pathway_filter <- subset(interaction_pathway_long, score >= 0.05 & pvalue <= 0.05)
  interaction_pathway_filter <- subset(interaction_pathway_long, pvalue < 0.05)

  data <- interaction_pathway_filter

  # interaction_TSC1 <- LR_pairs_one[ , c("Gene_secreted", "Gene_receptor", "pathway_receptor","0_TSC1>0_TSC1_pvalues", "0_TSC1>1_TSC1_pvalues", "0_TSC1>2_TSC1_pvalues", "1_TSC1>0_TSC1_pvalues", "1_TSC1>1_TSC1_pvalues", "1_TSC1>2_TSC1_pvalues", "2_TSC1>0_TSC1_pvalues", "2_TSC1>1_TSC1_pvalues", "2_TSC1>2_TSC1_pvalues")]
  # interaction_TSC1_pathway <- subset(interaction_TSC1, pathway_receptor == "PVR RTK signaling pathway")
  #
  # interaction_TSC1_pathway$pathway_receptor <- NULL
  # interaction_TSC1_pathway_long <- melt(interaction_TSC1_pathway, id.vars = c("Gene_secreted", "Gene_receptor"))
  # interaction_TSC1_pathway_filter <- subset(interaction_TSC1_pathway_long, value < 0.1)
  # data <- interaction_TSC1_pathway_filter

  # cell_col<-structure(c("#E5E5AF", "#FF66A1", "#BC85A9", "#E76172", "#CCD7D7",
  #                       "#FFA667", "#BC5DBB", "#76EA8E", "#90559F", "#5F9858",
  #                       "#B494D0", "#D5C766", "#959592", "#7CCFF9", "#AF876D",
  #                       "#F9CCDF", "#939DD1", "#4B5FF5", "#6BAFAE"), names= celltypes)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  myColors <- sample(col_vector, length(celltypes))
  # myColors <- brewer.pal(length(celltypes), "Set1")
  cell_col<-structure(myColors, names= celltypes)

  col <- cell_col
  label=FALSE
  edge.curved=0.5
  shape='circle'
  layout=in_circle()
  vertex.size=20
  margin=0.2
  vertex.label.cex=0.8
  vertex.label.color='black'
  arrow.width=3
  edge.label.color='black'
  edge.label.cex=1
  edge.max.width=4

  # net <- data %>% group_by(variable) %>% dplyr::summarize(n=n())
  net <- data %>% group_by(variable) %>% summarize(mean_score = mean(score))

  net <- net %>%
    separate(variable, c("sender", "receiver"), ">")

  # net$sender <- str_replace(net$sender, "_EGFP", "")
  # net$receiver <- str_replace(net$receiver, "_pvalues", "")

  empty_celltype <- setdiff(celltypes, unique(c(net$sender, net$receiver)))

  for(ct in empty_celltype) {
    # print(ct)
    line <- c(ct, ct, 0)
    # print(line)
    net <- rbind(net, line)
  }

  colnames(net) <- c("sender", "receiver", "n")
  net$n <- as.numeric(net$n)
  # net$sender <- str_replace(net$sender, "_TSC1", "")
  # net$receiver <- str_replace(net$receiver, "_TSC1_pvalues", "")

  net<-as.data.frame(net,stringsAsFactors=FALSE)
  g<-graph.data.frame(net,directed=TRUE)
  x <- get.adjacency(g, attr="n", sparse=FALSE)
  x <- x[celltypes, celltypes]
  # x[is.na(x)] <- 0
  # print(x)
  g <- graph_from_adjacency_matrix(x, mode = "directed", weighted = T)

  edge.start <- ends(g, es=E(g), names=FALSE)
  coords<-layout_(g,layout)
  if(nrow(coords)!=1){
    coords_scale=scale(coords)
  }else{
    coords_scale<-coords
  }
  loop.angle<-ifelse(coords_scale[V(g),1]>0,-atan(coords_scale[V(g),2]/coords_scale[V(g),1]),pi-atan(coords_scale[V(g),2]/coords_scale[V(g),1]))
  V(g)$size<-vertex.size
  V(g)$color<-col[V(g)]
  V(g)$label.color<-vertex.label.color
  V(g)$label.cex<-vertex.label.cex
  if(label){
    E(g)$label<-E(g)$n
  }
  if(max(E(g)$weight)==min(E(g)$weight)){
    E(g)$width<-1
  }else{
    # E(g)$width<-0.1 + edge.max.width/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
    E(g)$width <- 0.1 + E(g)$weight/max(E(g)$weight)*edge.max.width
  }
  E(g)$arrow.width<-arrow.width
  E(g)$label.color<-edge.label.color
  # E(g)$label.cex<-edge.label.cex
  E(g)$color<-V(g)$color[edge.start[,1]]
  if(sum(edge.start[,2]==edge.start[,1])!=0){
    E(g)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
  }

  png(file=paste0(output_dir, "/circleplot/circleplot_", pathway, ".png"),
      width     = 7,
      height    = 7,
      units     = "in",
      res       = 300,
  )
  plot(g,edge.curved=0.2,vertex.shape=shape,
       layout=coords_scale,margin=margin,edge.arrow.size=0.5, vertex.frame.color="white"
       , label=FALSE)
  dev.off()

  # png(file=paste0(output_dir, "/dotplot/dotplot_", pathway, ".png"), res=300, width = 1000, height = 1000)
  # plot(g,edge.curved=0.2,vertex.shape=shape,
  #      layout=coords_scale,margin=margin,edge.arrow.size=0.5, vertex.frame.color="white"
  #      , label=FALSE)
  # dev.off()

}




