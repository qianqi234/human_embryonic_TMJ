#load libraries
library(ggplot2)
library(dplyr)
library(monocle3)
library(Seurat)
library(ggforce)
options(future.globals.maxSize =  4000 * 1024^2)
setwd("F://human_embryonic_TMJ")

ggsave_sep <- function(plot_obj, filename, width = NA, height = NA,
                       use_cairo = FALSE, ...) {
  require(ggpubr)
  ifelse(use_cairo, cairo_pdf, pdf)(filename, width = width, height = height, ...)
  plot(plot_obj)
  plot(plot_obj + NoLegend())
  tryCatch(plot(as_ggplot(get_legend(plot_obj))),
           error = function(e) print(e),
           finally = dev.off())
}

# Load the data
object <- readRDS("Clustering.RDS")
tmp <- as.CellDataSet(object,assay = "RNA")
cds <- new_cell_data_set(as(object@assays$RNA@counts, "sparseMatrix"),
                         cell_metadata = tmp@phenoData@data,
                         gene_metadata = tmp@featureData@data)

#Pre-process the data
cds <- preprocess_cds(cds, num_dim = 40)
plot_pc_variance_explained(cds)
cds <- align_cds(cds, alignment_group = "orig.ident",
                 residual_model_formula_str = ~ S.Score + G2M.Score)

#Reduce dimensionality and visualize the results
cds <- reduce_dimension(cds, umap.n_neighbors = 60L)
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "seurat_clusters",
           label_branch_points=FALSE,label_leaves=FALSE, label_roots=FALSE,group_label_size = 8)

#Cluster your cells
cds <- cluster_cells(cds) 
plot_cells(cds, color_cells_by = "partition",group_label_size = 8,label_roots=FALSE)
plot_cells(cds, color_cells_by = "cluster",group_label_size = 8,label_roots=FALSE)

#Learn the trajectory graph
cds <- learn_graph(cds) 
p1<- plot_cells(cds,            
           color_cells_by = "seurat_clusters",            
           label_groups_by_cluster=FALSE,            
           label_leaves=FALSE,            
           label_branch_points=FALSE,group_label_size = 8)
p1
ggsave_sep(p1, "all_clusters_trajectories.pdf", width = 9, height = 9)
dev.off()
saveRDS(cds, "all_clusters_trajectories.RDS")

# Load the data
object2 <- subset(object, idents = c("Mesenchymal stem cells","Transition state","Hypertrophic chondrocytes","Osteoblasts","Tenocytes"))
tmp <- as.CellDataSet(object2,assay = "RNA")
cds_2 <- new_cell_data_set(as(object2@assays$RNA@counts, "sparseMatrix"),
                         cell_metadata = tmp@phenoData@data,
                         gene_metadata = tmp@featureData@data)

#Pre-process the data
cds_2 <- preprocess_cds(cds_2, num_dim = 70)
plot_pc_variance_explained(cds_2)
cds_2 <- align_cds(cds_2, alignment_group = "orig.ident",
                 residual_model_formula_str = ~ S.Score + G2M.Score)

#Reduce dimensionality and visualize the results
cds_2 <- reduce_dimension(cds_2, umap.n_neighbors = 110L)
plot_cells(cds_2, label_groups_by_cluster=FALSE,  color_cells_by = "seurat_clusters",
           label_branch_points=FALSE,label_leaves=FALSE, label_roots=FALSE,group_label_size = 8)

#Cluster your cells
cds_2 <- cluster_cells(cds_2) 
plot_cells(cds_2, color_cells_by = "partition",group_label_size = 8,label_roots=FALSE)
plot_cells(cds_2, color_cells_by = "cluster",group_label_size = 8,label_roots=FALSE)

#Learn the trajectory graph
cds_2 <- learn_graph(cds_2) 
p2<- plot_cells(cds_2,            
                color_cells_by = "seurat_clusters",            
                label_groups_by_cluster=FALSE,            
                label_leaves=FALSE,            
                label_branch_points=FALSE,group_label_size = 8)
p2
ggsave_sep(p2, "joint.pdf", width = 9, height = 9)
dev.off()
saveRDS(cds_2, "joint.rds")

#Rename the clusters
cds_2 <- readRDS("joint.rds")
ident.names <- list(`1`="Mesenchymal stem cells", `2`="Transition state cells",
                    `3`="Tenocytes",`6`="Hypertrophic chondrocytes",`13`="Osteoblasts")
cds_2@colData$label <- factor(ident.names[as.character(cds_2@colData$seurat_clusters)],
                            levels = ident.names)
p3 <- plot_cells(cds_2,            
                color_cells_by = "label",            
                label_groups_by_cluster=FALSE,            
                label_leaves=FALSE,            
                label_branch_points=FALSE,group_label_size = 8)
p3
ggsave_sep(p3, "joint_rename.pdf", width = 9, height = 9)
dev.off()

cds_3M <- cds_2[, colData(cds_2)$orig.ident == "3M_joint"]
p4 <- plot_cells(cds_3M,            
           color_cells_by = "seurat_clusters",            
           label_groups_by_cluster=FALSE,            
           label_leaves=FALSE,            
           label_branch_points=FALSE,group_label_size = 8)
p4
ggsave_sep(p4, "3M_joint.pdf", width = 9, height = 9)
dev.off()

cds_4M <- cds_2[, colData(cds_2)$orig.ident == "4M_joint"]
p5 <- plot_cells(cds_4M,            
           color_cells_by = "seurat_clusters",            
           label_groups_by_cluster=FALSE,            
           label_leaves=FALSE,            
           label_branch_points=FALSE,group_label_size = 8)
p5
ggsave_sep(p5, "4M_joint.pdf", width = 9, height = 9)
dev.off()

#gene expressing across the clusters
p6 <- plot_cells(cds_2,
                genes = "RUNX2",
                label_cell_groups=FALSE,
                group_label_size = 20,
                show_trajectory_graph=FALSE,
                min_expr = 0.1,
                cell_size = 1) +
  scale_color_gradient2(low = "grey", high = "#FF0000", mid = "#E9AAAA")
p6
ggsave("RUNX2.pdf", plot = p6, width = 8, height = 6)
dev.off()
p7 <- plot_cells(cds_2,
                genes = "SOX9",
                label_cell_groups=FALSE,
                group_label_size = 20,
                show_trajectory_graph=FALSE,
                min_expr = 0.1,
                cell_size = 1) +
  scale_color_gradient2(low = "grey", high = "#FF0000", mid = "#E9AAAA")
p7
ggsave("SOX9.pdf", plot = p7, width = 8, height = 6)
dev.off()
p8 <- plot_cells(cds_2,
                genes = "MMP13",
                label_cell_groups=FALSE,
                group_label_size = 20,
                show_trajectory_graph=FALSE,
                min_expr = 0.1,
                cell_size = 1) +
  scale_color_gradient2(low = "grey", high = "#FF0000", mid = "#E9AAAA")
p8
ggsave("MMP13.pdf", plot = p8, width = 8, height = 6)
dev.off()
p9 <- plot_cells(cds_2,
                genes = "SOST",
                label_cell_groups=FALSE,
                group_label_size = 20,
                show_trajectory_graph=FALSE,
                min_expr = 0.1,
                cell_size = 1) +
  scale_color_gradient2(low = "grey", high = "#FF0000", mid = "#E9AAAA")
p9
ggsave("SOST.pdf", plot = p9, width = 8, height = 6)
dev.off()
p10 <- plot_cells(cds_2,
                genes = "VEGFA",
                label_cell_groups=FALSE,
                group_label_size = 20,
                show_trajectory_graph=FALSE,
                min_expr = 0.1,
                cell_size = 1) +
  scale_color_gradient2(low = "grey", high = "#FF0000", mid = "#E9AAAA")
p10
ggsave("VEGFA.pdf", plot = p10, width = 8, height = 6)
dev.off()
p11 <- plot_cells(cds_2,
                genes = "COL10A1",
                label_cell_groups=FALSE,
                group_label_size = 20,
                show_trajectory_graph=FALSE,
                min_expr = 0.1,
                cell_size = 1) +
  scale_color_gradient2(low = "grey", high = "#FF0000", mid = "#E9AAAA")
p11
ggsave("COL10A1.pdf", plot = p11, width = 8, height = 6)
dev.off()
p12 <- plot_cells(cds_2,
                genes = "MAF",
                label_cell_groups=FALSE,
                group_label_size = 20,
                show_trajectory_graph=FALSE,
                min_expr = 0.1,
                cell_size = 1) +
  scale_color_gradient2(low = "grey", high = "#FF0000", mid = "#E9AAAA")
p12
ggsave("MAF.pdf", plot = p12, width = 8, height = 6)
dev.off()
p13 <- plot_cells(cds_2,
                genes = "CCDC80",
                label_cell_groups=FALSE,
                group_label_size = 20,
                show_trajectory_graph=FALSE,
                min_expr = 0.1,
                cell_size = 1) +
  scale_color_gradient2(low = "grey", high = "#FF0000", mid = "#E9AAAA")
p13
ggsave("CCDC80.pdf", plot = p13, width = 8, height = 6)
dev.off()
p14 <- plot_cells(cds_2,
                genes = "SYNE2",
                label_cell_groups=FALSE,
                group_label_size = 20,
                show_trajectory_graph=FALSE,
                min_expr = 0.1,
                cell_size = 1) +
  scale_color_gradient2(low = "grey", high = "#FF0000", mid = "#E9AAAA")
p14
ggsave("SYNE2.pdf", plot = p14, width = 8, height = 6)
dev.off()
