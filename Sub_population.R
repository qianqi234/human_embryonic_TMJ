#load libraries
library(sctransform)
library(cowplot)
library(Seurat)
library(SeuratObject)
library(GSEABase)
library(enrichplot)
library(clusterProfiler)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(ggpubr)
options(future.globals.maxSize =  4000 * 1024^2)# 40G

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

setwd("F://human_embryonic_TMJ")
object.new <- readRDS("Clustering.RDS")

#load cell cycle genes
s.genes=as.character(cc.genes.updated.2019$s.genes)
g2m.genes=as.character(cc.genes.updated.2019$g2m.genes)

#SCT normalize for each sample, regress by cell cycle scores and others
cluster2 <- subset(object.new, idents = "Transition state")
cluster2.list <- SplitObject(cluster2, split.by = "orig.ident")

#first normalized by 'percent.mt', 'nFeature_RNA', 'nCount_RNA'
cluster2.list <-
  lapply(cluster2.list,
         SCTransform,
         assay = 'RNA',
         new.assay.name = 'SCT',
         do.scale = T,
         do.center = T,
         vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA'))
#then cacluate cell cycle scores
cluster2.list <-
  lapply(cluster2.list,
         CellCycleScoring,
         s.features = s.genes,
         g2m.features = g2m.genes,
         assay = 'SCT',
         set.ident = TRUE)
#SCTranscform again by cell cycle scores and other three
cluster2.list <-
  lapply(cluster2.list,
         SCTransform,
         assay = 'RNA',
         new.assay.name = 'SCT',
         vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA',
                             'S.Score', 'G2M.Score'))

 cluster2.features <-
  SelectIntegrationFeatures(object.list = cluster2.list,
                            nfeatures = 3000)
cluster2.list <-
  PrepSCTIntegration(object.list = cluster2.list,
                     anchor.features = cluster2.features, 
                     verbose = FALSE)
cluster2.anchors <-
  FindIntegrationAnchors(object.list = cluster2.list,
                         normalization.method = "SCT", 
                         anchor.features = cluster2.features,
                         verbose = FALSE)
cluster2 <- IntegrateData(anchorset = cluster2.anchors,
                          normalization.method = "SCT", 
                          verbose = FALSE)

#Clustering
cluster2 <- RunPCA(cluster2, npcs = 50, verbose = FALSE)
cluster2 <- FindNeighbors(cluster2, dims = 1:40)
cluster2 <- FindClusters(cluster2, dims = 1:40, resolution = 0.1)
cluster2 <- RunUMAP(cluster2, dims = 1:50)
saveRDS(cluster2, "trantion_cluster.rds")

Cluster2.subcluster.joint.markers <- FindAllMarkers(cluster2,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Cluster2.subcluster.joint.markers, file = "Cluster2_markers.csv")
#Rename each cluster
ident.names <- list(`0`="Osteoblasts",`1`="Preosteoblasts",`2`="Hypertrophic chondrocytes",
                    `3`="Chondrocytes",
                    `4`="Mesenchymal stem cells")
cluster2.new <- RenameIdents(cluster2, ident.names)
p1 <- DimPlot(cluster2.new, reduction = "umap") %>% LabelClusters("ident")
p1
ggsave_sep(p1, "Cluster2.pdf", width = 9, height = 9)
dev.off()

#The Percentage of each cluster in Cluster2
col <- scales::hue_pal()(nlevels(object.new))
p2 <- (table(cluster2.new@active.ident) %>%
          prop.table() %>%
          as.data.frame %>%
          ggplot(aes(x=Var1,y=Freq)))+
  scale_fill_manual(values = col, name = "Cell Types") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = c(0, 0, 0.1, 0),
                     name = "Percentage") +
  stat_summary(geom = "text", position = position_dodge(width = 1),
               mapping = aes(label = scales::percent_format(accuracy = 0.1)(Freq),
                             y = Freq + 0.01)) +
  geom_col(aes(fill = Var1), width = 0.8, size = 0) +
  cowplot::theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1),
        axis.title.x = element_blank())
p2
ggsave("clusters2_Percentage.pdf", p2, width = 12, height = 12)
dev.off()

# The Heatmap of marker genes for annotation 
dotplot_markers <- c("OSTN","COL14A1","IBSP", "COL9A2", "GFRA1")
p3 <- DotPlot(cluster2.new, features = dotplot_markers, ) +
  scale_y_discrete(limits = rev) +
  viridis::scale_color_viridis() +
  guides(size = guide_legend(title = "Percentage"),
         color = guide_colorbar(title = "Expression")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.1)) +
  ylab("Cluster")
p3
ggsave_sep(p3, "clusters2_marker_genes_annotation.pdf", width =13, height = 13)
dev.off()