#load libraries
library(pheatmap)
library(ComplexHeatmap)
library(dplyr)
library(igraph)
library(InteractiveComplexHeatmap)

options(future.globals.maxSize =  4000 * 1024^2)# 40G
setwd("F://human_embryonic_TMJ")

factorize_matrix <- function(x, labels) {
  mtx <- factor(x, labels = labels)
  class(mtx) <- c("matrix", "factor")
  dim(mtx) <- dim(x)
  dimnames(mtx) <- dimnames(x)
  invisible(mtx)
}

# Load the data
binarized_mtx <- read.csv("binarized_joint.txt",
                          row.names=1, check.names = F)
seurat_clusters <- readRDS("seurat_clusters.rds")
rss_matrix <- read.csv("regulon_specificity_bootstrap.csv",
                       row.names=1, check.names = F)
regulon_sizes <- read.table("regulon_sizes.txt",
                             row.names = 1, header = F, col.names = c(NA, "Counts"))

#Pre-process the data
get_top_n <- function(x, n) names(x)[sort.list(x, decreasing = TRUE)[1:n]]
top_regulons <- asplit(rss_matrix, 1) %>% lapply(get_top_n, n = 5)%>%
  `[`(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15))
regulons <- unlist(top_regulons, use.names = FALSE)
regulons <- intersect(regulons, colnames(binarized_mtx))
regulons_group <- top_regulons %>%
  {structure(rep(names(.), lengths(.)), names = unlist(.))} %>%
  `[`(., regulons) %>%
  as.numeric()
barcodes <- split(names(seurat_clusters), seurat_clusters) %>%
  `[`(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)) %>%
  unlist(use.names = FALSE)
sub_matrix <- t(binarized_mtx[barcodes, regulons])
clusters <- seurat_clusters[colnames(sub_matrix)]
colors = list("Regulon OFF" = "#F8F8F8FF", "Regulon ON" = "#000000FF")
colors_clusters <- structure(paste0(scales::hue_pal()(15), "FF"),
                             names=levels(seurat_clusters))
name_fix <- function(text) sprintf("%s (%ig)", text, regulon_sizes[text, 1])
top_annot <- HeatmapAnnotation(Clusters = clusters, show_annotation_name = FALSE,
                               col = list(Clusters = colors_clusters))

#Heatmap
pdf("top_5_regulons_Heatmap.pdf", width = 12, height = 9)
Heatmap(factorize_matrix(sub_matrix, c("Regulon OFF", "Regulon ON")),
        show_column_names = FALSE, col = colors, name = "Regulon Activity",
        row_labels = name_fix(regulons),
        row_names_gp = gpar(fontsize = 6), column_title = NULL, row_title = NULL,
        top_annotation = top_annot, # right_annotation = right_annot,
        cluster_rows = FALSE, cluster_columns = FALSE,
        row_gap = unit(0, "mm"), column_gap = unit(0, "mm"),
        #row_split = regulons_group, column_split = seurat_clusters[barcodes],
        clustering_distance_columns = "binary", show_column_dend = FALSE,
        cluster_column_slices = FALSE) %>%
  draw(merge_legend = TRUE, use_raster = TRUE)
dev.off()

#module analysis
compute_pcc <- function(df) {
  df <- as.matrix(1 * (df > 0))
  a <- df %*% t(df)
  b <- df %*% (1 - t(df))
  c <- (1 - df) %*% t(df)
  d <- ncol(df) - a - b - c
  return((a * d - b * c)/sqrt((a + b) * (a + c) * (b + d) * (d + c)))
}
compute_csi <- function(mat){
  apply(mat, 1, function(x) sapply(colnames(mat), function(y) {
    (sum((x < x[y] - 0.05) | (mat[y,] < x[y] - 0.05)) + 1) / length(x)
  }))
}
pcc <- compute_pcc(t(binarized_mtx))
csi <- compute_csi(pcc)
g = graph_from_adjacency_matrix(csi >= 0.7)
write.graph(g, "regulon.gml", format = "gml")
MCL_cluster <- read.csv("MCL_cluster.csv") %>%
  with(structure(Cluster, names = label))
colors_clusters <- structure(c("#FFA431FF","#00C5E5FF","#00AAB8FF","#A5C900FF",
                               "#D8818FFF","#D6C93AFF","#FF2D43FF","#009F84FF","#00BFFF","#7FFFAA",
                               "red", "blue"),
                             names=1:12)
top_annot_csi <- HeatmapAnnotation(`MCL Cluster` = as.factor(MCL_cluster),
                                   show_annotation_name = FALSE,
                                   col = list(`MCL Cluster` =colors_clusters))
ht_module <- Heatmap(csi, name = "CSI",
                     show_column_names = FALSE, show_row_names = FALSE, 
                     col = circlize::colorRamp2(c(0.2, 1), c("#FCFCFF", "#0000FF")),
                     column_title = NULL, row_title = NULL,
                     clustering_distance_columns = function(x) dist(x >= 0.7, "binary"),
                     clustering_distance_rows = function(x) dist(x >= 0.7, "binary"),
                     row_gap = unit(0, "mm"), column_gap = unit(0, "mm"),
                     top_annotation = top_annot_csi)
pdf("module.pdf", width = 12, height = 12)
p <- draw(ht_module, merge_legend = TRUE, use_raster = TRUE)
dev.off()

