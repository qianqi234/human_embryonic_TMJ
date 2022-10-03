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

#load human embryonic temporomandibular joint (TMJ) matrix data
data1_1<- Read10X(data.dir = "RTMJ_matrix")
object1_1 <- CreateSeuratObject(data1_1, project = "3M_joint", min.cells = 3, min.features =200,assay = "RNA")
object1_1 <- PercentageFeatureSet(object1_1, pattern = "^MT-", col.name = "percent.mt")
data1_2=read.table("G5_matrix.tsv",sep="\t",header=T, row.names = 1)
object1_2 <- CreateSeuratObject(data1_2, project = "4M_joint", min.cells = 3, min.features =200,assay = "RNA")
object1_2 <- PercentageFeatureSet(object1_2, pattern = "^MT-", col.name = "percent.mt")
object=merge(object1_1,object1_2, c("3M_joint", "4M_joint"))

#filter
object <- subset(object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000& percent.mt < 50& nCount_RNA<20000)

#QC plot
p1<- VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p1
ggsave_sep(p1, "quality_control.pdf", width = 9, height = 9)
dev.off()

#load cell cycle genes
s.genes=as.character(cc.genes.updated.2019$s.genes)
g2m.genes=as.character(cc.genes.updated.2019$g2m.genes)

#SCT normalize for each sample, regress by cell cycle scores and others
object.list = SplitObject(object)

#first normalized by 'percent.mt', 'nFeature_RNA', 'nCount_RNA'
for (i in 1:2) {
    print(object.list[[i]])
  object.list[[i]]<- SCTransform(object.list[[i]],
    assay = 'RNA',
    new.assay.name = 'SCT',
    do.scale=T,
    do.center=T,
    vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA')
  )
#then cacluate cell cycle scores
  object.list[[i]]<- CellCycleScoring(
    object.list[[i]],
    s.features = s.genes,
    g2m.features = g2m.genes,
    assay = 'SCT',
    set.ident = TRUE
  )
#SCTranscform again by cell cycle scores and other three
  object.list[[i]] <- SCTransform(
    object.list[[i]],
    assay = 'RNA',
    new.assay.name = 'SCT',
    vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'S.Score', 'G2M.Score')
  )
}
object.features <- SelectIntegrationFeatures(object.list = object.list, nfeatures = 3000)
object.list <- PrepSCTIntegration(object.list = object.list, anchor.features = object.features, 
                                  verbose = FALSE)
object.anchors <- FindIntegrationAnchors(object.list = object.list, normalization.method = "SCT", 
                                  anchor.features = object.features, verbose = FALSE)
object.integrated <- IntegrateData(anchorset = object.anchors, normalization.method = "SCT", 
                                   verbose = FALSE)

#Clustering
pc.num = 50
object.integrated <- RunPCA(object.integrated, npcs = pc.num, verbose = FALSE)
object.integrated <- RunUMAP(object.integrated, reduction = "pca", dims = 1:50)
object.integrated <- FindNeighbors(object.integrated, dims = 1:pc.num, verbose = FALSE)
object.integrated <- FindClusters(object.integrated, resolution =0.25,  verbose = FALSE)

#The Heatmap of top five DEGs of each cluster
object.integrated.markers <- FindAllMarkers(object.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
object.integrated.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
object.integrated.markers %>% write.csv("marker_genes.csv", quote = F)
top_5_markers <- object.integrated.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
p2<- DoHeatmap(object.integrated, features = top_5_markers$gene) +
  scale_fill_gradient2(low = "#92C4EA", high = "#DC5F44", mid = "white")
p2
ggsave_sep(p2, "top_5_markers_heatmap.pdf", width = 9, height = 9)
dev.off()

#Rename each cluster
ident.names <- list(`0`="Satellite cells",`1`="Mesenchymal stem cells", `2`="Transition state",
                    `3`="Tenocytes",`4`="Myoblasts",`5`="Endothelial cells",`6`="Hypertrophic chondrocytes",
                    `7`="Erythrocytes",`8`="Proliferating cells",`9`="Leukocytes",`10`="Pericytes",
                    `11`="Chondrocytes",`12`="Schwann cells",`13`="Osteoblasts",
                      `14`="Osteoclasts")
object.new <- RenameIdents(object.integrated, ident.names)
p3<- DimPlot(object.new, reduction = "umap") %>% LabelClusters("ident")
p3
ggsave_sep(p3, "cluster_names.pdf", width = 9, height = 9)
dev.off()
p4<- DimPlot(object.new, reduction = "umap", split.by  = "orig.ident") %>% LabelClusters("ident")
p4
ggsave_sep(p4, "cluster_split_ident.pdf", width = 9, height = 9)
dev.off()
saveRDS(object.new, "Clustering.RDS")

# The Heatmap of marker genes for annotation 
dotplot_markers <- c("PAX7", "MYF5","PDGFRA","NRK","THBS2","PTN", "KERA", "TNMD","MYOG", "KLHL41","CDH5", "PECAM1", 
                     "LOXL4", "COL10A1","HBA1", "HBA2","SGO1","HIST1H4C","PTPRC", "CXCL8", "RGS5", "ABCC9",
                     "CYTL1", "ITGA10","PLP1", "MPZ", "IBSP", "ALPL","RGS10", "MMP9","MKI67","TOP2A" )
p5<- DotPlot(object.new, features = dotplot_markers, ) +
  scale_y_discrete(limits = rev) +
  viridis::scale_color_viridis() +
  guides(size = guide_legend(title = "Percentage"),
         color = guide_colorbar(title = "Expression")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1.1)) +
  ylab("Cluster")
 p5
 ggsave_sep(p5, "marker_genes_annotation.pdf", width = 9, height = 9)
 dev.off()
 
#The Number of cells of each cluster
col2 <- c("3M_joint" = "#006400", "4M_joint" = "#FA8072")
p6 <- (table(object.new@active.ident, object.new$orig.ident) %>%
          as.data.frame %>%
          ggplot(aes(x=Var1,y=Freq, fill=Var2)))+
  scale_fill_manual(values = col2, name = "Sample", 
                    labels = c("3M_Joint", "4M_Joint")) +
  scale_y_continuous(expand = c(0, 0, 0.1, 0),
                     name = "Number of cells") +
  stat_summary(geom = "text", position = position_dodge(width = 1),
               mapping = aes(y = Freq + 0.01, label=Freq)) +
  geom_col(width = 0.85, size = 0, position = "dodge") +
  cowplot::theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1),
        axis.title.x = element_blank())
p6
ggsave_sep(p6, "Number_of_cells_of_each_cluster.pdf", width = 9, height = 9)
dev.off()

#The GO enrichment analysis
library("clusterProfiler")
library(stringr)
genes=as.character(object.integrated.markers$gene)
cid2name=str_pad(levels(object.integrated$seurat_clusters), 2, side = "left", pad = "0")
names(cid2name)=levels(object.integrated$seurat_clusters)
cid2name=as.list(cid2name)
eg = bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
name2id=list()
for(i in c(1:nrow(eg))){
  name=as.character(eg[i,1])
  id=as.character(eg[i,2])
  name2id[[name]]=id
}
gene.f=genes %in% eg$SYMBOL
genes_keep=genes[gene.f]
gene_id=sapply(genes_keep,function(x){name2id[[x]]})
gene_id=as.character(gene_id)
clusters=object.integrated.markers$cluster
my.df=data.frame(Entrez=gene_id,cluster=clusters[gene.f])
ont_type="BP"
fun_type="enrichGO"
formula_res <- compareCluster(Entrez~cluster, data=my.df, fun=fun_type ,OrgDb="org.Hs.eg.db",ont=ont_type,pvalueCutoff=0.01,pAdjustMethod="BH")
formula_res@compareClusterResult$cluster <-
  factor(formula_res@compareClusterResult$cluster, levels = levels(object.integrated))
p7<- dotplot(formula_res, x=~cluster, showCategory=5)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
p7
ggsave_sep(p7, "clusters_GO.pdf", width = 9, height = 9)
dev.off()

# The GSEA enrichment analysis
log2fc <- FoldChange(object.new,
                     ident.1 = "Mesenchymal stem cells",
                     ident.2 = "Transition state")
geneList <- sort(setNames(log2fc$avg_log2FC, rownames(log2fc)), decreasing = T)
geneset <- read.gmt("genesets_human_differentiation_proliferation.gmt")
egmt <- GSEA(geneList, TERM2GENE = geneset, verbose = T, pvalueCutoff = 1)

p8<- gseaplot2(egmt, geneSetID = "WP_OSTEOBLAST_DIFFERENTIATION",
          title = "WP_OSTEOBLAST_DIFFERENTIATION")
p8
ggsave_sep(p8, "WP_OSTEOBLAST_DIFFERENTIATION.pdf",width = 13, height = 14)
dev.off()

log2fc <- FoldChange(object.new,
                     ident.1 = "Transition state",
                     ident.2 = "Osteoblasts")
geneList <- sort(setNames(log2fc$avg_log2FC, rownames(log2fc)), decreasing = T)
geneset <- read.gmt("genesets_bone.gmt")
egmt <- GSEA(geneList, TERM2GENE = geneset, verbose = T, pvalueCutoff = 1)

p9<- gseaplot2(egmt, geneSetID = "GOBP_BONE_DEVELOPMENT",
               title = "GOBP_BONE_DEVELOPMENT")
p9
ggsave_sep(p9, "GOBP_BONE_DEVELOPMENT.pdf",width = 13, height = 14)
dev.off()
p10<- gseaplot2(egmt, geneSetID = "GOBP_OSSIFICATION",
               title = "GOBP_OSSIFICATION")
p10
ggsave_sep(p10, "GOBP_OSSIFICATION.pdf",width = 13, height = 14)
dev.off()
#The expression of PTN, THBS1 and CAPN6 in each cluster 
p11=FeaturePlot(object.new, feature="PTN" ,
               min.cutoff = 0.1,
               max.cutoff = 6
)
p11
ggsave_sep(p11, "PTN.pdf", width = 9, height = 9)
dev.off()

p12=FeaturePlot(object.new, feature="THBS1" ,
               min.cutoff = 0.1,
               max.cutoff = 12)
p12
ggsave_sep(p12, "THBS1.pdf", width = 9, height = 9)
dev.off()

p13=FeaturePlot(object.new, feature="CAPN6" ,
               min.cutoff = 0.1,
               max.cutoff = 8)
p13
ggsave_sep(p13, "CAPN6.pdf", width = 9, height = 9)
dev.off()