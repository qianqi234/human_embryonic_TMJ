# Loading required package
library(CellChat)
library(Seurat)
library(NMF)
library(ggalluvial)
library(tidyverse)
library(ggalluvial)
library(ggpubr)
library(ggalluvial)
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

rds = readRDS('Clustering.RDS')
rds <- NormalizeData(rds,verbose = FALSE, normalization.method = "LogNormalize", assay="RNA",scale.factor = 1e6) #log1p(RPM) 
rds <- ScaleData(rds,do.scale = TRUE , do.center = TRUE, assay="RNA") 
levels(rds)
assignInNamespace("scPalette", scales::hue_pal(), "CellChat")

# Create a CellChat object
meta = GetAssayData(rds, slot='data', assay='RNA')
meta[1:5,1:5]
idents = data.frame(row.names=rownames(rds@meta.data), celltype=Idents(rds)) 
head(idents)
cellchat <- createCellChat( meta , meta=idents , group.by='celltype' )
cellchat
levels( cellchat@idents )

# Load and set the required CellChatDB database
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
head(CellChatDB$interaction)
dplyr::glimpse(CellChatDB$interaction)
unique(CellChatDB$interaction$pathway_name)
CellChatDB.ss <- subsetDB(CellChatDB,  key='annotation') 
cellchat@DB <- CellChatDB.ss
cellchat@DB <- CellChatDB

# Preprocessed expression data for cell-to-cell interaction analysis
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

# Inferring cellular interaction networks
cellchat <- computeCommunProb(cellchat, raw.use = FALSE, do.fast = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
df.net
head(df.net)
cellchat <- aggregateNet(cellchat) 
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,1))

# Number of interactions of clusters
opar <- par()
par(mar = c(2, 2, 2, 2))
p1<- netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Number of interactions")
pdf(file = "Number of interactions.pdf", width = 9, height = 9)
p1
dev.off()
do.call(par, opar)

# Interaction strength of clusters
p2<- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Interaction strength")
pdf(file = "Interaction strength.pdf", width=9, height=9)
p2
dev.off()
do.call(par,opar)

# Extracting communication relationships among specific cell types
df.net2 <- subsetCommunication(cellchat, sources.use=c('Mesenchymal stem cells','Transition state','Hypertrophic chondrocytes',
                                                       'Osteoblasts','Tenocytes'),
                               targets.use= c('Mesenchymal stem cells','Transition state','Hypertrophic chondrocytes',
                                              'Osteoblasts','Tenocytes'))
df.net2

# Extracting the signaling pathway of FGF among cell types
df.net3 <- subsetCommunication(cellchat, signaling = c("FGF"))
head(df.net3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat@netP

# Integrate communication network of the signaling pathway of FGF
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mar = c(0, 0, 0, 0), mai = c(0.1, 0, 0, 0.1), mfrow = c(1,1))
p3<- netVisual_aggregate(cellchat, signaling = "FGF", vertex.receiver = c(1, 3, 4, 11, 10) )
pdf(file = "FGF_netVisual_aggregate.pdf", width = 20, height = 10)
p3
dev.off()
do.call(par, opar)
par(mar = c(2, 2, 2, 2), mai = c(0.5, 0.5, 0.5, 0.5), mfrow = c(1,1))
p4<- netVisual_aggregate(cellchat, signaling = "FGF", layout = "circle")
pdf(file = "FGF_circle.pdf", width = 9, height = 9)
p4
dev.off()
do.call(par, opar)

#The degree of interaction contribution of the signaling pathway of FGF
p <- netAnalysis_contribution(cellchat, signaling = "FGF")
p$data <- p$data[p$data$contribution > 0, ]
p5<- p + scale_x_discrete(limit=rev)
pdf(file = " p + scale_x_discrete.pdf", width = 9, height = 9)
p5
dev.off()

# Integrate communication network of the signaling pathway of FGF7_FGFR1
par(mar = c(0.2, 0.2, 0.2, 0.2), mai = c(0.1, 0.1, 0.1, 0.1), mfrow = c(1,1))
p6<- netVisual_individual(
  cellchat,
  signaling = "FGF",
  pairLR.use = "FGF7_FGFR1",
  vertex.receiver = c(1, 3, 4, 11, 10)
)
pdf(file = "FGF7_FGFR1_pathway.pdf", width = 30, height = 10)
p6
dev.off()
do.call(par, opar)
par(mfrow=c(1,1)) 
levels(cellchat@idents)
ct_group = c('SC','MSC','TS','TES',
             'MYS','ENC','HC','ERS',
             'PC','LTS','PYS','CYS',
             'SCS','OBTS','OC') 
names(ct_group) = levels(cellchat@idents) 
ct_group
strwidth <- function(x){0.3}
p7<- netVisual_individual(cellchat, signaling= "FGF", pairLR.use= "FGF7_FGFR1" , layout = "chord" )
pdf(file = "FGF7_FGFR1_pathway_chord.pdf", width = 9, height = 9)
p7
dev.off()
do.call(par, opar)

# The expression of FGF signaling molecules
plotGeneExpression(cellchat, signaling = "FGF")
p8 <- plotGeneExpression(cellchat, signaling = "FGF", enriched.only = FALSE)
p8$patches <- within(p8$patches, {
  plots[[length(plots)]] <- plots[[length(plots)]] + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
})
p8
pdf(file = "FGF_expression_molecules.pdf", width = 9, height = 9)
p8
dev.off()
do.call(par, opar)

# Computing network centrality
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Heatmap displaying
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", width = 40, height = 40, font.size = 10)
pdf(file = "outgoing netAnalysis.pdf", width = 20, height = 20)
ht1
dev.off()
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", width = 40, height = 40, font.size = 10)
pdf(file = "incoming netAnalysis.pdf", width = 20, height = 20)
ht2
dev.off()
selectK(cellchat, pattern = "outgoing")
dev.off()
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
dev.off()
p9<- netAnalysis_river(cellchat, pattern = "outgoing")
pdf(file = "outgoing communication patterns.pdf", width = 18, height = 18)
p9
dev.off()
do.call(par, opar)
selectK(cellchat, pattern = "incoming")
dev.off()
nPatterns = 5
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
dev.off()
p10<- netAnalysis_river(cellchat, pattern = "incoming")
pdf(file="incoming communcation patterns.pdf",width=18, height=18)
p10
dev.off()
do.call(par,opar)