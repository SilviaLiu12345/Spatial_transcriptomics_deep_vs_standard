## new UMAP shape by the DEGs

library(spatstat.explore)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

rm(list=ls())

options(future.globals.maxSize= 891289600000)

SAMPLENAME=c("J259B","J259X")


data=list() # data list


### read in the data
for(sInd in 1:length(SAMPLENAME)){
  print(c(sInd, SAMPLENAME[sInd]))
  ## read in data
  data_dir <- paste0("../code/",SAMPLENAME[sInd],"/outs/")
  list.files(data_dir) # Should show filtered_feature_bc_matrix.h5
  data[[sInd]] <- Load10X_Spatial(data.dir = data_dir)
  
  data[[sInd]]$library=SAMPLENAME[sInd]
  
  print(data[[sInd]])
}
names(data)=SAMPLENAME

# [1] "1"     "J259B"
# An object of class Seurat
# 18085 features across 2397 samples within 1 assay
# Active assay: Spatial (18085 features, 0 variable features)
# 1 image present: slice1
# [1] "2"     "J259X"
# An object of class Seurat
# 18085 features across 3923 samples within 1 assay
# Active assay: Spatial (18085 features, 0 variable features)
# 1 image present: slice1


####### Data pre-processing:
for(sInd in 1:length(data)){
  print(c(sInd, SAMPLENAME[sInd]))
  
  ## data visualization
  pdf(paste0("Vln.",sInd,".",SAMPLENAME[sInd],".pdf"), width=10, height=5)
  plot1 <- VlnPlot(data[[sInd]], features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
  plot2 <- SpatialFeaturePlot(data[[sInd]], features = "nCount_Spatial") + theme(legend.position = "right")
  print(wrap_plots(plot1, plot2))
  dev.off()
  
  ## filter out zero cells (add this step, incase SCTransform will show error)
  # data[[i]] <- subset(x = data[[i]], subset = nCount_Spatial > 0)
  # print(data[[i]])
  
  ## transformation
  data[[sInd]] <- SCTransform(data[[sInd]], assay = "Spatial", verbose = FALSE)
  
  ## Gene expression visualization
  # ALLGENE=rownames(data[[i]]@assays$SCT)
  # plotGene=c("Alb","Gapdh")
  # plotGene%in%ALLGENE
  #
  # pdf(paste0("Feature.",i,".pdf"), width=10, height=5)
  # print(SpatialFeaturePlot(data[[i]], features = plotGene))
  # dev.off()
  
}

# save(data, file="data.rdata")
# 
# load("data.rdata")


#### integration, pipeline 1, as suggested by ST pipeline, but no normalization
# OBJ.merge=merge(x=data[[1]], y=data[[2]], ## list all
#                 add.cell.ids=SAMPLENAME,
#                 project="Spatial")
# 
# OBJ.merge
# # An object of class Seurat 
# # 36126 features across 6320 samples within 2 assays 
# # Active assay: SCT (18041 features, 0 variable features)
# # 1 other assay present: Spatial
# # 2 images present: slice1, slice1_J259X
# 
# DefaultAssay(OBJ.merge) <- "SCT"
# 
# tmp=c()
# for(i in 1:length(data)){
#   tmp=c(tmp, VariableFeatures(data[[i]]))
# }
# VariableFeatures(OBJ.merge) <- tmp
# 
# OBJ.merge <- RunPCA(OBJ.merge, verbose = FALSE)
# OBJ.merge <- FindNeighbors(OBJ.merge, dims = 1:30)
# OBJ.merge <- FindClusters(OBJ.merge, verbose = FALSE, resolution=0.8)
# table(OBJ.merge@meta.data$seurat_clusters)
# 
# OBJ.merge <- RunUMAP(OBJ.merge, dims = 1:30)
# 
# ## add in the library information
# LIB=sapply(strsplit(rownames(OBJ.merge@meta.data), split="_"), function(x) return(x[1]))
# OBJ.merge@meta.data$library=LIB
# 
# table(OBJ.merge@meta.data$library)
# # J259B J259X 
# # 2397  3923 
# 
# saveRDS(file = "OBJ.merge.RDS", OBJ.merge)
# 
# OBJ.merge <- readRDS("OBJ.merge.RDS")
# 
# 
# table(OBJ.merge@meta.data$seurat_clusters, OBJ.merge@meta.data$library )
# 
# # J259B J259X
# # 0      0  1107
# # 1      9   899
# # 2      6   822
# # 3    685     2
# # 4    656    14
# # 5    385    85
# # 6    405     0
# # 7      5   378
# # 8      0   325
# # 9    246     0
# # 10     0   150
# # 11     0   141


## additional pre-processing for Seurat V5
for(sInd in 1:length(data)){
  print(c(sInd, SAMPLENAME[sInd]))
  data[[sInd]]@assays[["RNA"]] <- data[[sInd]]@assays[["Spatial"]]
}


############# pipeline 2: SCTransform, as scRNA-seq data

### Preform the Data integration between two samples using DEGs 

## read in DEGs

DEfile1=read.table("DE.adj_surround.vs.adj_nonImmune.txt", header=T, sep="\t")
gene1=unique(DEfile1$Gene[DEfile1$p_val<0.01])
length(gene1) # 424

DEfile2=read.table("DE.tumor_surround.vs.tumor_nonImmune.txt", header=T, sep="\t")
gene2=unique(DEfile2$Gene[DEfile2$p_val<0.01])
length(gene2) # 1202

#select.features=unique(c(gene1, gene2))
select.features=gene1
select.features=intersect(select.features, rownames(data[[1]]@assays$SCT))
length(select.features) # 424

select.features2 <- SelectIntegrationFeatures(object.list = data, nfeatures = 3000)
library(openxlsx)
write.xlsx(data.frame(Gene=select.features2), file="AVITI_full_3000.xlsx")

select.features=intersect(select.features, select.features2)
length(select.features) # 182

data <- PrepSCTIntegration(object.list = data,
                           anchor.features = select.features, 
                           verbose = TRUE)


select.anchors <- FindIntegrationAnchors(object.list = data, 
                                         normalization.method = "SCT", 
                                         anchor.features = select.features,
                                         verbose = T,reference = 1,)

OBJ.merge <- IntegrateData(anchorset = select.anchors, normalization.method = "SCT", 
                           verbose = FALSE,)


## Generate the PCA, UMAP and visualizations
OBJ.merge <- RunPCA(OBJ.merge, verbose = FALSE)
OBJ.merge <- RunUMAP(OBJ.merge, dims = 1:20,n.neighbors = 30)

DimPlot(OBJ.merge, label=T,group.by = c("library"), combine = FALSE)


##Generate initial cluster calls
OBJ.merge <- FindNeighbors(OBJ.merge, reduction = "pca", dims = 1:20)
OBJ.merge <- FindClusters(OBJ.merge, resolution = 1) ## old 9 clusters with resolution=0.5

### swap and scale to RNA assay for visuals
DefaultAssay(OBJ.merge) <- "Spatial"
OBJ.merge <- FindVariableFeatures(OBJ.merge)
OBJ.merge <- NormalizeData(OBJ.merge)
OBJ.merge <- ScaleData(OBJ.merge)

FeaturePlot(OBJ.merge,features=c("ALB","GAPDH"),split.by="library")

names(OBJ.merge@images)=SAMPLENAME

nrow(OBJ.merge@meta.data) #  6320


## load in the Immune, tumor label
OBJ.old=readRDS("OBJ.merge.immune.RDS")
nrow(OBJ.old@meta.data) # 6320
sum(rownames(OBJ.merge@meta.data)==rownames(OBJ.old@meta.data)) # 6320

OBJ.merge@meta.data$ImmuneLabel=OBJ.old@meta.data$ImmuneLabel
OBJ.merge@meta.data$TumorLabel=OBJ.old@meta.data$TumorLabel


saveRDS(OBJ.merge, file = "OBJ.merge.DEV2.RDS")
# OBJ.merge=readRDS("OBJ.merge.DEV2.RDS")

### plot

pdf(paste0("DimPlot.newDEV2.pdf"), width=10, height=5)
DimPlot(OBJ.merge, reduction = "umap", group.by = c("seurat_clusters", "library"))
dev.off()

pdf(paste0("DimPlot.newDEV2.ImmuneLabel.pdf"), width=6, height=5)
DimPlot(OBJ.merge, reduction = "umap", group.by = c("ImmuneLabel")) +
  scale_color_manual(values=c("Immune"="red", "SurroundImmune"="blue", "NonImmune"="lightgrey"))
dev.off()

pdf(paste0("DimPlot.newDEV2.TumorLabel.pdf"), width=6, height=5)
DimPlot(OBJ.merge, reduction = "umap", group.by = c("TumorLabel"))
dev.off()


### dot plot of 182 genes

# expression across 8 clusters
pdf(paste0("DotPlot.182gene.cluster.pdf"), width=8, height=30)
DotPlot(object=OBJ.merge, features=select.features, group.by = "seurat_clusters") +
  coord_flip()
dev.off()


# expression across immune and tumor categories
Idents(OBJ.merge)=paste0(OBJ.merge@meta.data$TumorLabel,"__", OBJ.merge@meta.data$ImmuneLabel)
pdf(paste0("DotPlot.182gene.tumor_immune.pdf"), width=8, height=30)
DotPlot(object=OBJ.merge, features=select.features) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 30, hjust=1))
dev.off()



####### check gene expression, spatial

DefaultAssay(OBJ.merge) <- "Spatial"

ALLGENE=rownames(OBJ.merge@assays$SCT)


# Hepatocytes
# ALB, TF, TTR, CYP3A4, APOA1
plotGene=c("ALB", "TF", "TTR", "CYP3A4", "APOA1") 
plotGene[!plotGene%in%ALLGENE]

pdf(paste0("FeaturePlot.newDEV2.Hepatocytes.pdf"), width=10, height=10)
FeaturePlot(OBJ.merge, features = plotGene, ncol=2)
dev.off()


# Liver sinusoidal ECs
# LYVE1, ENG, STAB1, PECAM1, CD36
plotGene=c("LYVE1", "ENG", "STAB1", "PECAM1", "CD36") 
plotGene[!plotGene%in%ALLGENE]


pdf(paste0("FeaturePlot.newDEV2.SinusoidalEC.pdf"), width=10, height=10)
FeaturePlot(OBJ.merge, features = plotGene)
dev.off()


# Vascular ECs
# PECAM1, CCL21, MMRN1, GNG11, FLT4
plotGene=c("PECAM1", "CCL21", "MMRN1", "GNG11", "FLT4") 
plotGene[!plotGene%in%ALLGENE]


pdf(paste0("FeaturePlot.newDEV2.VascularEC.pdf"), width=10, height=10)
FeaturePlot(OBJ.merge, features = plotGene)
dev.off()


# Cholangiocytes
# EPCAM, KRT7, KRT19, CLDN4, TACSTD2
plotGene=c("EPCAM", "KRT7", "KRT19", "CLDN4", "TACSTD2") 
plotGene[!plotGene%in%ALLGENE]


pdf(paste0("FeaturePlot.newDEV2.Cholangiocytes.pdf"), width=10, height=10)
FeaturePlot(OBJ.merge, features = plotGene)
dev.off()


# Stellate cells
# ACTA2, DCN, BGN, COLEC11, IGFBP3
plotGene=c("ACTA2", "DCN", "BGN", "COLEC11", "IGFBP3") 
plotGene[!plotGene%in%ALLGENE]


pdf(paste0("FeaturePlot.newDEV2.StellateCells.pdf"), width=10, height=10)
FeaturePlot(OBJ.merge, features = plotGene)
dev.off()


# Kupffer cells
# CD68, CD163, LYZ, C1QA, AIF1
plotGene=c("CD68", "CD163", "LYZ", "C1QA", "AIF1") 
plotGene[!plotGene%in%ALLGENE]

pdf(paste0("FeaturePlot.newDEV2.KupfferCells.pdf"), width=10, height=10)
FeaturePlot(OBJ.merge, features = plotGene)
dev.off()



# T cells
# CD3D, CD2, IL7R, TRBC2, CD69
plotGene=c("CD3D", "CD2", "IL7R", "TRBC2", "CD69") 
plotGene[!plotGene%in%ALLGENE]

pdf(paste0("FeaturePlot.newDEV2.Tcells.pdf"), width=10, height=10)
FeaturePlot(OBJ.merge, features = plotGene)
dev.off()


# B cells / Plasma cells
# IGKC, JCHAIN, CD79A, CD27, CD74
plotGene=c("IGKC", "JCHAIN", "CD79A", "CD27", "CD74") 
plotGene[!plotGene%in%ALLGENE]

pdf(paste0("FeaturePlot.newDEV2.Bcells.pdf"), width=10, height=10)
FeaturePlot(OBJ.merge, features = plotGene)
dev.off()


# CD genes and NK cell markers
# Check CD4, CD8, CD11
# NK cell markers: Nkg7, Klrd1, Prf1, Cd7, Trdc

## CD11 -> ITGAM

plotGene=c("CD4","CD8A","ITGAM",
           "NKG7", "KLRD1", "PRF1", "CD7", "TRDC") 
plotGene[!plotGene%in%ALLGENE]


pdf(paste0("FeaturePlot.newDEV2.CDgene.NKcell.pdf"), width=10, height=10)
FeaturePlot(OBJ.merge, features = plotGene)
dev.off()





