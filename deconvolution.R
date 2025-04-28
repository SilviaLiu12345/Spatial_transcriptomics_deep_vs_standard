
##### Luo deconvolution CARD
# Spatially informed cell-type deconvolution for spatial transcriptomics
# https://pubmed.ncbi.nlm.nih.gov/35501392/
# https://yma-lab.github.io/CARD/documentation/04_CARD_Example.html


#### spatial
# The AVITI-sequencing data:
#   /ix/sliu/Luo/LoopSpatial/Seurat.yes/OBJ.merge.immune.RDS
# 
# The NextSeq data
# /ix/sliu/Luo/Spatial2023/Seurat/OBJ.merge.immune.RDS
# 
# My meta.data contains the immune labelling and tumor/normal labelling.

### sc reference needed annotation
### 


#devtools::install_github('YMa-lab/CARD')
library(CARD)
packageVersion('CARD') ## ‘1.1’
library(spatstat.explore)
library(Seurat)
packageVersion("Seurat") ## ‘5.1.0’
library(SeuratData)
library(tidyverse)
library(patchwork)
packageVersion("patchwork") ## ‘1.2.0’
library(ggpubr)

options(future.globals.maxSize= 1e12)

getwd()
setwd('/ix/sliu/Luo/LoopSpatial/deconvolution_s.yes/')

#### load sp mtx
#load("./spatial_count.RData")
#spatial_count[1:4,1:4]
sp_aviti <- readRDS("/ix/sliu/Luo/LoopSpatial/Seurat.yes/OBJ.merge.immune.RDS")
sp_aviti
# An object of class Seurat 
# 39126 features across 6320 samples within 3 assays 
# Active assay: Spatial (18085 features, 2000 variable features)
# 3 layers present: counts, data, scale.data
# 2 other assays present: SCT, integrated
# 2 dimensional reductions calculated: pca, umap
# 2 images present: J259B, J259X

sp_nextseq <- readRDS("/ix/sliu/Luo/Spatial2023/Seurat/OBJ.merge.immune.RDS")
sp_nextseq
# An object of class Seurat 
# 39122 features across 6320 samples within 3 assays 
# Active assay: Spatial (18085 features, 2000 variable features)
# 3 layers present: counts, data, scale.data
# 2 other assays present: SCT, integrated
# 2 dimensional reductions calculated: pca, umap
# 2 images present: J259B, J259X

table(sp_aviti$ImmuneLabel)
table(sp_nextseq$ImmuneLabel)

getwd()
colour = c(ADJ = "deepskyblue", TUMOR = "firebrick4")

DimPlot(sp_aviti)
pdf("SpatialDimPlot_aviti.pdf", height = 5, width = 8)
SpatialDimPlot(sp_aviti, pt.size.factor = 80, group.by = "TumorLabel", cols = colour, alpha = .5)
dev.off()

DimPlot(sp_nextseq)
pdf("SpatialDimPlot_nextseq.pdf", height = 5, width = 8)
SpatialDimPlot(sp_nextseq, pt.size.factor = 80, group.by = "TumorLabel", cols = colour, alpha = .5)
dev.off()

pdf("SpatialDimPlot_aviti.immune.pdf", height = 5, width = 8)
SpatialDimPlot(sp_aviti, pt.size.factor = 80, group.by = "ImmuneLabel")
dev.off()

#sp_aviti <- IntegrateLayers(object = sp_aviti, method = HarmonyIntegration)

table(sp_aviti@meta.data[["library"]])
# J259B J259X 
# 2397  3923 

sub_J259B <- subset(sp_aviti, subset = library == "J259B")
sub_J259X <- subset(sp_aviti, subset = library == "J259X")

########## nextseq
table(sp_nextseq$library)

sub_J259B <- subset(sp_nextseq, subset = library == "J259B")
sub_J259X <- subset(sp_nextseq, subset = library == "J259X")

sp_count <- sub_J259B@assays[["Spatial"]]@data
sp_count[1:4,1:4]

sp_count <- sub_J259X@assays[["Spatial"]]@data
sp_count[1:4,1:4]


#### load sp location x,y coord
#load("./spatial_location.RData")
#spatial_location[1:4,]

spatial_location <- data.frame(cell = sp_aviti@images$J259B@coordinates)%>%select(cell.imagerow, cell.imagecol)
spatial_location <- data.frame(cell = sp_aviti@images$J259X@coordinates)%>%select(cell.imagerow, cell.imagecol)


spatial_location <- data.frame(cell = sp_nextseq@images$J259B@coordinates)%>%select(cell.imagerow, cell.imagecol)
spatial_location <- data.frame(cell = sp_nextseq@images$J259X@coordinates)%>%select(cell.imagerow, cell.imagecol)

#spatial_location <- rbind(spatial_location_J259B, spatial_location_J259X)%>%
#  select(cell.imagerow, cell.imagecol)

colnames(spatial_location)[1] <- 'x'
colnames(spatial_location)[2] <- 'y'
head(spatial_location)



### sc count mtx
#load("./sc_count.RData")
#sc_count[1:4,1:4]

sc <- readRDS('/ix/sliu/jun/project/Sungjin_SALL4_datamining/data/Human_GSE189903/OBJ.GSE189903.hcc_6samples.RDS')
sc
# An object of class Seurat 
# 55354 features across 61297 samples within 3 assays 
# Active assay: RNA (33538 features, 2000 variable features)
# 3 layers present: data, counts, scale.data
# 2 other assays present: SCT, integrated
# 2 dimensional reductions calculated: pca, umap

DimPlot(sc, group.by = "newCluster", label = T)

#table(sc$lib)

sc_sub <- subset(sc, subset = lib %in% c("HB", "HT"))

sc_count <- sc_sub@assays[["RNA"]]$data
sc_count[1:4,1:4]

ncol(sc_count) ## 61297
nrow(sc_count) ## 33538

### sc cell type
#load("./sc_meta.RData")
#sc_meta[1:4,]

sc_meta <- data.frame(sc_sub$newCluster)%>%
  rownames_to_column(var = "cell")
colnames(sc_meta)[2] <- 'cellType'
rownames(sc_meta) <- sc_meta$cell
sc_meta$sampleInfo <- "sample1"
head(sc_meta)


print(nrow(sc_meta)) ## 35653 



### check input
Dim(sp_count) ## 18085  6320
Dim(spatial_location) ## 6320    2
Dim(sc_count) ## 33538 61297
Dim(sc_meta) ## 61297     3

### debug
#
#load("./spatial_count.RData")
#load("./spatial_location.RData")
#load("./sc_count.RData")
#load("./sc_meta.RData")
#


CARD_object = createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = sp_count,
  spatial_location = spatial_location,
  ct.varname = "cellType",
  ct.select = unique(sc_meta$cellType),
  sample.varname = "sampleInfo") 
#CARD_object

CARD_obj = CARD_deconvolution(CARD_object = CARD_object)

## error
# Basis_ref = createscRef(sc_eset, ct.select, ct.varname, sample.varname)


print(CARD_obj@Proportion_CARD[1:2,])


# res_CARD = as.data.frame(CARD_obj@Proportion_CARD)
# location = CARD_obj@spatial_location
# data = cbind(res_CARD,location)
# 
# radius = (max(data$x) - min(data$x)) * (max(data$y) - min(data$y))
# radius = radius / nrow(data)
# radius = radius / pi
# radius = sqrt(radius) * 0.85 >> 269

table(sub_J259B@meta.data[["ImmuneLabel"]])

df_immune <- data.frame(sub_J259B$ImmuneLabel)%>%filter(sub_J259B.ImmuneLabel == "Immune")
df_immune <- data.frame(sub_J259X$ImmuneLabel)%>%filter(sub_J259X.ImmuneLabel == "Immune")

head(df_immune)
head(CARD_obj@Proportion_CARD[, 'Hepatocytes'], 20)
head(CARD_obj@Proportion_CARD[, 'Bcells'], 20)
head(sub_J259X$ImmuneLabel, 20)


CARD_obj@Proportion_CARD[rownames(CARD_obj@Proportion_CARD) %in% setdiff(rownames(CARD_obj@Proportion_CARD), rownames(df_immune)), 'Hepatocytes'] <- 1
CARD_obj@Proportion_CARD[, 'Bcells'][CARD_obj@Proportion_CARD[, 'Hepatocytes'] == 1] <- 0
CARD_obj@Proportion_CARD[, 'Tcells'][CARD_obj@Proportion_CARD[, 'Hepatocytes'] == 1] <- 0
CARD_obj@Proportion_CARD[, 'Myeloid'][CARD_obj@Proportion_CARD[, 'Hepatocytes'] == 1] <- 0
CARD_obj@Proportion_CARD[, 'NKcells'][CARD_obj@Proportion_CARD[, 'Hepatocytes'] == 1] <- 0
CARD_obj@Proportion_CARD[, 'Proliferative'][CARD_obj@Proportion_CARD[, 'Hepatocytes'] == 1] <- 0


#### output the immune cell composition
immune_J259B <- data.frame(CARD_obj@Proportion_CARD)%>%rownames_to_column(var = "cells")%>%
  filter(cells %in% rownames(df_immune))

immune_J259X <- data.frame(CARD_obj@Proportion_CARD)%>%rownames_to_column(var = "cells")%>%
  filter(cells %in% rownames(df_immune))


aviti <- rbind(immune_J259B, immune_J259X)
nextseq <- rbind(immune_J259B, immune_J259X)

getwd()
write.table(av_tumor, "aviti_immune_prop.csv", quote = F, sep = ',', row.names = F)
write.table(nextseq_tumor, "nextseq_immune_prop.csv", quote = F, sep = ',', row.names = F)

#### add tumor annotation
aviti <- read.delim("aviti_immune_prop.csv", header = T, sep = ',')
nextseq <- read.delim("nextseq_immune_prop.csv", header = T, sep = ',')

av_tumor <- data.frame(sp_aviti$TumorLabel)%>%rownames_to_column(var = "cells")%>%
  right_join(.,aviti, by = "cells")

nextseq_tumor <- data.frame(sp_nextseq$TumorLabel)%>%rownames_to_column(var = "cells")%>%
  right_join(.,nextseq, by = "cells")



colors = c("chartreuse3","gray90","orangered2","#FFD92F","#377EB8","#F0027F") ##A
# colors = c("chartreuse3","floralwhite","firebrick2","#FFD92F","#377EB8","#F0027F") ##B 
# colors = c("springgreen4","snow1","red1","#FFD92F","#377EB8","#F0027F") ##C
# colors = c("limegreen","#F0F8FF","red1","#FFD92F","dodgerblue","palevioletred1") ##D

getwd()
pdf("card_nextseq.J259X.immune.pdf", height = 6, width = 7)
CARD.visualize.pie(
  proportion = CARD_obj@Proportion_CARD,
  spatial_location = CARD_obj@spatial_location, 
  colors = colors,
  radius = 10)+ ### You can choose radius = NULL or your own radius number
  theme_void()
dev.off()

tmp <- data.frame(CARD_obj@Proportion_CARD)


### output df
##### AVITI
output <- data.frame(sp_aviti$ImmuneLabel, sp_aviti$TumorLabel)%>%
  rownames_to_column(var = "cells")
head(output)

aviti <- data.frame(CARD_obj@Proportion_CARD)%>%rownames_to_column(var = "cells")
aviti1 <- data.frame(CARD_obj@Proportion_CARD)%>%rownames_to_column(var = "cells")

aviti_prop <- rbind(aviti, aviti1)
aviti_prop <- inner_join(output, aviti_prop, by = "cells")

aviti_prop_immune <- aviti_prop%>%
  filter(sp_aviti.ImmuneLabel %in% c("Immune"))%>%
  gather("celltype", "prop", -c("cells", "sp_aviti.ImmuneLabel", "sp_aviti.TumorLabel"))%>%
  group_by(celltype)


ggplot(aviti_prop_immune, mapping = aes(x = celltype, y = prop, fill = sp_aviti.TumorLabel))+
  geom_violin(alpha = .5)+
  theme_classic()+
  facet_wrap(~celltype, scale = "free")+
  stat_compare_means(method = "wilcox.test")+
  labs(x = "")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

head(nextseq_prop)

immu_adj <- nextseq_prop%>%filter(sp_nextseq.ImmuneLabel == "Immune")%>%
  filter(sp_nextseq.TumorLabel == "ADJ")

immu_tumor <- nextseq_prop%>%filter(sp_nextseq.ImmuneLabel == "Immune")%>%
  filter(sp_nextseq.TumorLabel == "TUMOR")

head(immu_adj,3)
head(immu_tumor,3)

b <- wilcox.test(x=immu_adj$Bcells, y=immu_tumor$Bcells)$p.value
t <- wilcox.test(x=immu_adj$Tcells, y=immu_tumor$Tcells)$p.value
h <- wilcox.test(x=immu_adj$Hepatocytes, y=immu_tumor$Hepatocytes)$p.value
m <- wilcox.test(x=immu_adj$Myeloid, y=immu_tumor$Myeloid)$p.value
n <- wilcox.test(x=immu_adj$NKcells, y=immu_tumor$NKcells)$p.value
p <- wilcox.test(x=immu_adj$Proliferative, y=immu_tumor$Proliferative)$p.value

list <- c(b,t,h,m,n,p)
list
## 1.762134e-18 4.449185e-29 4.083855e-36 1.866337e-37 1.216984e-34 1.205995e-30


##### NETseq
output <- data.frame(sp_nextseq$ImmuneLabel, sp_nextseq$TumorLabel)%>%
  rownames_to_column(var = "cells")
head(output)

nextseq <- data.frame(CARD_obj@Proportion_CARD)%>%rownames_to_column(var = "cells")
nextseq1 <- data.frame(CARD_obj@Proportion_CARD)%>%rownames_to_column(var = "cells")

nextseq_prop <- rbind(nextseq, nextseq1)
nextseq_prop <- inner_join(output, nextseq_prop, by = "cells")

nextseq_prop_immune <- nextseq_prop%>%
  filter(sp_nextseq.ImmuneLabel %in% c("Immune"))%>%
  gather("celltype", "prop", -c("cells", "sp_nextseq.ImmuneLabel", "sp_nextseq.TumorLabel"))%>%
  group_by(celltype)

head(nextseq_prop_immune)

ggplot(nextseq_prop_immune, mapping = aes(x = celltype, y = prop, fill = sp_nextseq.TumorLabel))+
  geom_violin(alpha = .5)+
  theme_classic()+
  facet_wrap(~celltype, scale = "free")+
  stat_compare_means(method = "wilcox.test")+
  labs(x = "")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


table(sc$newCluster)
ct.visualize = c('NKcells', 'Tcells', 'Bcells',  'Myeloid', 'Hepatocytes', 'Proliferative')

## visualize the spatial distribution of the cell type proportion
#pdf("prop_nc_all.pdf", height = 12, width = 22)
CARD.visualize.prop(
  proportion = CARD_obj@Proportion_CARD,        
  spatial_location = CARD_obj@spatial_location, 
  ct.visualize = ct.visualize,                 ### selected cell types to visualize
  colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
  NumCols = 6,                                 ### number of columns in the figure panel
  pointSize = 1.0)+                            ### point size in ggplot2 scatterplot  
  theme_classic()+
  labs(x = "", y = "")
#dev.off()









######## not run

CARD.visualize.prop.2CT(
  proportion = CARD_obj@Proportion_CARD,                             ### Cell type proportion estimated by CARD
  spatial_location = CARD_obj@spatial_location,                      ### spatial location information
  ct2.visualize = c("PT","Myofib"),              ### two cell types you want to visualize
  colors = list(c("lightblue4","lightyellow3","red4"),c("lightblue","lightyellow","red")))+       ### two color scales      
  theme_void()
print(p3)

CARD.visualize.Cor(CARD_obj@Proportion_CARD,colors = NULL) # if not provide, we will use the default colors
print(p4)




CARD_obj = CARD.imputation(CARD_obj,NumGrids = 2000,ineibor = 10,exclude = NULL)
location_imputation = cbind.data.frame(x=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",1)),
                                       y=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",2)))
rownames(location_imputation) = rownames(CARD_obj@refined_prop)
library(ggplot2)


CARD.visualize.prop(
  proportion = CARD_obj@refined_prop,                         
  spatial_location = location_imputation,            
  ct.visualize = ct.visualize,                    
  colors = c("lightblue","lightyellow","red"),    
  NumCols = 6)


ALLGENE=rownames(sp@assays$RNA)

genelist <- c("Ecm1", "Fbln1", "Fbln2", "Fbln5", "Fbln7", "Mfap2", "Mfap4", "Mfap5", "Nid1", "Nid2", 
              "Tnc", "Vtn", "Col1a1", "Col1a2", "Col3a1")

genelist[!genelist%in%ALLGENE]
## "Cyr61", "Mfap1", "Mfap3"

pdf("card.gene.pdf", height = 11, width = 21)
CARD.visualize.gene(
  spatial_expression = CARD_obj@refined_expression,
  spatial_location = location_imputation,
  gene.visualize = genelist,
  colors = NULL,
  NumCols = 6)+
  theme_classic()+
  labs(x = "", y = "")
dev.off()


CARD.visualize.gene(
  spatial_expression = CARD_obj@spatial_countMat,
  spatial_location = CARD_obj@spatial_location,
  gene.visualize = genelist,
  colors = NULL,
  NumCols = 6)
print(p8)



######## not run 
#### error
scMapping = CARD_SCMapping(CARD_obj,shapeSpot="Square",numCell=20,ncore=10)
library(SingleCellExperiment)
MapCellCords = as.data.frame(colData(scMapping))
count_SC = assays(scMapping)$counts

df = MapCellCords
colors = c("#8DD3C7","#CFECBB","#F4F4B9","#CFCCCF","#D1A7B9","#E9D3DE","#F4867C","#C0979F",
           "#D5CFD6","#86B1CD","#CEB28B","#EDBC63","#C59CC5","#C09CBF","#C2D567","#C9DAC3","#E1EBA0",
           "#FFED6F","#CDD796","#F8CDDE")
ggplot(df, aes(x = x, y = y, colour = CT)) + 
  geom_point(size = 3.0) +
  scale_colour_manual(values =  colors) +
  #facet_wrap(~Method,ncol = 2,nrow = 3) + 
  theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        panel.background = element_rect(colour = "white", fill="white"),
        plot.background = element_rect(colour = "white", fill="white"),
        legend.position="bottom",
        panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.title =element_blank(),
        legend.title=element_text(size = 13,face="bold"),
        legend.text=element_text(size = 12),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.key.size = unit(0.45, 'cm'),
        strip.text = element_text(size = 15,face="bold"))+
  guides(color=guide_legend(title="Cell Type"))
print(p10)




## Spatial transcriptomics deconvolution at single-cell resolution using Redeconve
## https://www.nature.com/articles/s41467-023-43600-9

devtools::install_github("ZxZhou4150/Redeconve", build_vignettes = F)

library(Redeconve)
packageVersion("Redeconve") ## ‘1.1.2’
library(spatstat.explore)
library(Seurat)
packageVersion("Seurat") ## ‘5.1.0’
library(SeuratData)
InstallData("pbmc3k")
InstallData("stxBrain")

## sc
pbmc3k <- UpdateSeuratObject(pbmc3k)
sclist <- extract.seurat.sc(pbmc3k)

## st
stxBrain <- LoadData("stxBrain", type = "anterior1")
stxBrain[['RNA']] <- stxBrain[['Spatial']]

stlist <- extract.seurat.st(stxBrain@assays[["RNA"]]@layers$counts)

