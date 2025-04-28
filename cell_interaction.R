
##### Luo spatial ref mapping
library(spatstat.explore)
library(tidyverse)
library(anndata)
library(Seurat)
packageVersion("Seurat") ## 5.1.0
library(SeuratObject)
packageVersion("SeuratObject") ## 5.0.2
library(SeuratData)
library(tools) ## for toTitleCase

rm(list=ls())
gc()


setwd('/ix/sliu/Luo/LoopSpatial/cpdb/')
getwd()

options(future.globals.maxSize= 1e12)



#### cell chat
library(CellChat)
library(patchwork)

library(Seurat)
library(ggplot2)
library(GenomicRanges)
library(cowplot)
library(dplyr)
library(ggnewscale)

table(sp_aviti$library)
table(sp_nextseq$library)

sp_aviti$newCluster <- paste(sp_aviti$TumorLabel, sp_aviti$ImmuneLabel, sep = '_')
sp_nextseq$newCluster <- paste(sp_nextseq$TumorLabel, sp_nextseq$ImmuneLabel, sep = '_')

sp_aviti$newCluster[sp_aviti$newCluster == "ADJ_Immune"] <- "Benign\nImmune cells"
sp_aviti$newCluster[sp_aviti$newCluster == "ADJ_NonImmune"] <- "Benign\naway from immune cells"
sp_aviti$newCluster[sp_aviti$newCluster == "ADJ_SurroundImmune"] <- "Benign\nAdjacent to immune cells"
sp_aviti$newCluster[sp_aviti$newCluster == "TUMOR_Immune"] <- "HCC\nImmune cells"
sp_aviti$newCluster[sp_aviti$newCluster == "TUMOR_NonImmune"] <- "HCC\naway from immune cells"
sp_aviti$newCluster[sp_aviti$newCluster == "TUMOR_SurroundImmune"] <- "HCC\nAdjacent to immune cells"

table(sp_aviti$newCluster)

#### NEXTSEQ
sp_nextseq$newCluster[sp_nextseq$newCluster == "ADJ_Immune"] <- "Benign\nImmune cells"
sp_nextseq$newCluster[sp_nextseq$newCluster == "ADJ_NonImmune"] <- "Benign\naway from immune cells"
sp_nextseq$newCluster[sp_nextseq$newCluster == "ADJ_SurroundImmune"] <- "Benign\nAdjacent to immune cells"
sp_nextseq$newCluster[sp_nextseq$newCluster == "TUMOR_Immune"] <- "HCC\nImmune cells"
sp_nextseq$newCluster[sp_nextseq$newCluster == "TUMOR_NonImmune"] <- "HCC\naway from immune cells"
sp_nextseq$newCluster[sp_nextseq$newCluster == "TUMOR_SurroundImmune"] <- "HCC\nAdjacent to immune cells"

table(sp_nextseq$newCluster)

#J259B <- subset(sp_aviti, subset = library == "J259B")
#J259X <- subset(sp_aviti, subset = library == "J259X")

table(sp_nextseq$newCluster)

##### creat cellchat obj
cellchat=createCellChat(sp_nextseq, group.by="newCluster")
head(cellchat@meta)

CellChatDB <- CellChatDB.human   # use CellChatDB.mouse if running on mouse data; CellChatDB.human for human
showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# project gene expression data onto PPI (Optional: when running it, 
# USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
#cellchat <- projectData(cellchat, PPI.human) # PPI.human for human


cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)

#df.net_aviti <- subsetCommunication(cellchat) 
df.net_nextseq <- subsetCommunication(cellchat) 

getwd()
library(openxlsx)
library(readxl)
#write.xlsx(df.net_aviti, "Ligand.Receptor.all.aviti.xlsx")
write.xlsx(df.net_nextseq, "Ligand.Receptor.all.nextseq.xlsx")

df_aviti <- read_excel("Ligand.Receptor.all.aviti.xlsx")
df_nextseq <- read_excel("Ligand.Receptor.all.nextseq.xlsx")

length(unique(df_aviti$pathway_name))



cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

pdf("Number.of.interaction.all.aviti.pdf", width = 6.5, height = 6)
p <- netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 #color.use = c("Immune" = "#F8766D", "NonImmune" = "#00BA38", "SurroundImmune" = "#619CFF"),
                 label.edge= F, title.name = "Number of interactions")
print(p)
dev.off()

pdf("Number.of.interaction.all.nextseq.empty.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 #color.use = c("Immune" = "#F8766D", "NonImmune" = "#00BA38", "SurroundImmune" = "#619CFF"),
                 label.edge= F, vertex.label.cex = "none", title.name = "Number of interactions")
dev.off()

 

pdf("interaction.per.celltype.all.nextseq.pdf")
mat <- cellchat@net$weight
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), 
                   #color.use = c("Immune" = "#F8766D", "NonImmune" = "#00BA38", "SurroundImmune" = "#619CFF"),
                   title.name = rownames(mat)[i])
}
dev.off()


#### pathway
head(cellchat@meta)

# Signaling analysis: CCL, CXCL, TNF Interferons (IFNs), Interleukins CTLA-4/B7 PD-1/PD-L1

PATH_aviti=sort(unique(df.net_aviti$pathway_name))
PATH_nextseq=sort(unique(df.net_nextseq$pathway_name))

diff <- setdiff(PATH_aviti, PATH_nextseq)

diff <- setdiff(PATH_nextseq, PATH_aviti)

library(eulerr)
### venn plot
fit=euler(c("a" = length(PATH_aviti), "b" = length(PATH_nextseq),
            "a&b" = length(intersect(PATH_aviti, PATH_nextseq))), 
          input="union")

pdf("Venn.cellchat.label.pdf",height=6,width=6)
plot(fit,edges=FALSE,labels=c(paste("AVITI n=",length(PATH_aviti),sep=""),
                              paste("NextSeq n=",length(PATH_nextseq),sep="")),
     quantities=T,fills=c("#dec6ee","#ffe796"))
dev.off()

pdf("Venn.cellchat.pdf",height=6,width=6)
plot(fit,edges=FALSE,labels=c(" ", " "),
     quantities=T,fills=c("#dec6ee", "#ffe796"))
dev.off()


fit=euler(c("a" = length(unique(df_aviti$pathway_name)), "b" = length(unique(df_nextseq$pathway_name)),
            "a&b" = length(intersect(unique(df_aviti$pathway_name), unique(df_nextseq$pathway_name)))), 
          input="union")

pdf("Venn.cellchat.empty.pdf",height=6,width=6)
plot(fit,edges=FALSE,labels=c(paste(""),
                              paste("")),
     quantities=F,fills=c("#dec6ee","#ffe796"))
dev.off()


for(pathways.show in PATH_nextseq){
  print(pathways.show)
  pdf(paste0("Pathway.",pathways.show,"_nextseq.pdf"))
  print(netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord"))
  dev.off()
}



netVisual_chord_cell(cellchat, signaling = "ADGRA", 
                     title.name = paste0(pathways.show, "signaling network"))

netVisual_aggregate(cellchat, signaling = "ADGRA", layout = "chord")




#### NOT RUN
#### cellphoneDB
library(SeuratDisk)

getwd()
SaveH5Seurat(sp_aviti, filename = "sp_aviti.h5Seurat")

Convert("sp_aviti.h5Seurat", dest = "h5ad")




######## NOT RUN mapping to ref sc dataset

##### ref 
sc <- readRDS('/ix/sliu/jun/project/Sungjin_SALL4_datamining/data/Human_GSE189903/OBJ.GSE189903.hcc_6samples.RDS')
sc

DimPlot(sc, split.by = "lib")

sc <- SCTransform(sc)


##### query 
sp_aviti <- readRDS("/ix/sliu/Luo/LoopSpatial/Seurat.yes/OBJ.merge.immune.RDS")
sp_aviti

sp_nextseq <- readRDS("/ix/sliu/Luo/Spatial2023/Seurat/OBJ.merge.immune.RDS")
sp_nextseq

DimPlot(sp_aviti)

#list <- toupper(c(#"TF", "TTR", "CYP3A4", "APOA1",
#  "ALB",
#                  "CD68", "CD163", "LYZ", "C1QA", "AIF1",
#           "CD3D", "CD2", "IL7R", "TRBC2", "CD69",
#           "CD79A", "CD27", "CD74",
#           "CD4", "CD8", "CD11",
#           "Nkg7", "Klrd1", "Prf1", "Cd7", 'Trdc'))
#
#
#list <- list[list%in%rownames(sc@assays$SCT)] 
#list <- list[list%in%rownames(sp_aviti@assays$SCT)] 

#### mapping
anchors <- FindTransferAnchors(reference = sc, query = sp_aviti, dims = 1:30, reference.reduction = "pca", 
                               features = list,
                               normalization.method = 'SCT')
predictions <- TransferData(anchorset = anchors, refdata = sc$newCluster, dims = 1:30) #for my data column labels in scRNA had celltype annotation, you can change it as per your data 
sp_aviti <- AddMetaData(sp_aviti, metadata = predictions)
sp_aviti

DimPlot(sp_aviti, group.by = "predicted.id")
DimPlot(sp_aviti, group.by = "ImmuneLabel")







