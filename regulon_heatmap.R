#load libraries
library(Seurat)
library(SeuratObject)
library(clusterProfiler)
library(ggplot2)
library(pheatmap)
library(ggpubr)
options(future.globals.maxSize =  4000 * 1024^2)# 40G
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

object <- readRDS("object_with_score.rds")
binarized_mtx <- read.csv("binarized_joint.txt",
                          header = T, check.names = F, row.names = 1)
colnames(binarized_mtx) <- gsub("[(+)]", "", colnames(binarized_mtx))
auc_mtx <- object@meta.data[, startsWith(colnames(object@meta.data), "Score_")]
colnames(auc_mtx) <- sub("Score_", "", colnames(auc_mtx))
object[["regulons"]] <- CreateAssayObject(data = t(as.matrix(auc_mtx)))
object[["regulons_binary"]] <- CreateAssayObject(data = t(as.matrix(binarized_mtx)))
object <- ScaleData(object, assay = "regulons", do.scale = FALSE, do.center = TRUE)
object <- RunPCA(object, assay = "regulons", features = colnames(auc_mtx), npcs = 50, reduction.name = "PCA_regulons")
p1<- DoHeatmap(object, c("ZSCAN32", "NRF1", "LBX2", "MAFA", "PPARA", "ZNF354B", "PAX1", "TFEB", "ETV4", "CEBPZ", "SALL2", "ZSCAN29",
                    "ZNF287", "MEIS1", "ZNF48", "NFYA", "YY2", "SRY", "BRF1", "ZDHHC15", "TBX15", "SHOX2", "PLAGL1", "TCF7L2", 
                    "JUN", "PBX1", "TBX18", "NFIA", "ARNT", "NR1I2", "TWIST1", "DLX1", "ALX4", "PAX9", "DLX4", "ZNF148", "GLI2",
                    "DLX6", "ZNF607", "MSX1", "ZNF157", "RARA", "DLX2", "JUNB", "HIC1", "CPEB1", "SREBF1", "UBTF", "GTF3C2", "ZNF658", 
                    "POU2F3", "ZNF672", "POU3F1", "NR2F2", "FOXD3", "POU3F2", "NR6A1"), assay = "regulons",
          disp.min = 0, disp.max = 0.1, label = F) + 
scale_fill_gradient(low = "#F5F5F577", high = "blue",name = "AUC")
ggsave("heatmap_regulons.pdf", plot = p1, width = 15, height = 9)
dev.off()

