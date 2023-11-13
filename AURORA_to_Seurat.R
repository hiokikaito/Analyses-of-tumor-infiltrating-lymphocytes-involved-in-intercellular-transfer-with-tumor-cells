#################################################
# Name: scRNAseq_CD45+_TILs.R
#
# Description: R codes used for the analysis of multicolor flow cytometry data (from Cytek AURORA) on Tumor Infiltrating Lymphocytes collected from B16 or B16 ZsGreen+ tumors (also for comparison of MC38 and MC38 ZsGreen+ tumors).
#
# Input: Flow cytometry data converted into 10x Genomics file format, using the "csv_to_10x.py" script.
# 
# Date updated: Nov. 13, 2023
# 
# Author: Kaito Hioki 
# Contact: khioki@umass.edu
#################################################



#install.packages("dplyr")
library(dplyr)
#install.packages('Seurat')
library(Seurat)
#install.packages('patchwork')
library(patchwork)
library(ggplot2)



########################################
# Data preparation
########################################
###First the B16 file
B16.data <- Read10X(data.dir = "B16_10x_format/")
B16 <- CreateSeuratObject(counts = B16.data, project = "AURORA", min.cells = 0, min.features = 0, names.field = 2, names.delim = "-")
zsvals <- read.csv(file = "B16_TILs_ZsGreen.csv", header = F)
B16$ZsGreen <- zsvals$V2
dapivals <- read.csv(file = "B16_TILs_DAPI.csv", header = F)
B16$DAPI <- dapivals$V2
B16$orig.ident <- "B16"
head(B16)

###Next the B16ZS file
B16ZS.data <- Read10X(data.dir = "B16ZS_10x_format/")
B16ZS <- CreateSeuratObject(counts = B16ZS.data, project = "AURORA", min.cells = 0, min.features = 0, names.field = 2, names.delim = "-")
zsvals <- read.csv(file = "B16ZS_TILs_ZsGreen.csv", header = F)
B16ZS$ZsGreen <- zsvals$V2
dapivals <- read.csv(file = "B16ZS_TILs_DAPI.csv", header = F)
B16ZS$DAPI <- dapivals$V2
B16ZS$orig.ident <- "B16ZS"
head(B16ZS)






########################################
# Merge B16 and B16ZS samples, annotate together. 
# Integration method to combine two datasets fails because there are not enough variables to compare.
########################################
aurora <- merge(B16, B16ZS)
head(aurora)
table(aurora$orig.ident)

#The number of features is set to the number of relevant variables measured by the flow cytometer.
aurora <- FindVariableFeatures(aurora, selection.method = "vst", nfeatures = 15)
all.genes <- rownames(aurora)
aurora <- ScaleData(aurora, features = all.genes)
aurora <- RunPCA(aurora, features = VariableFeatures(object = aurora), approx = F)
aurora <- FindNeighbors(aurora, dims = 1:15)
aurora <- FindClusters(aurora, resolution = 0.5)
aurora <- RunUMAP(aurora, dims = 1:15)
# aurora <- RunTSNE(object = aurora)
head(aurora)
table(aurora$orig.ident)


#Plot the merged dataset
DimPlot(aurora, reduction = "umap", pt.size = 0.1, raster = T)
DimPlot(subset(aurora, subset = orig.ident == "B16"), reduction = "umap", pt.size = 0.1, raster = T)
DimPlot(subset(aurora, subset = orig.ident == "B16ZS"), reduction = "umap", pt.size = 0.1, raster = T)








########################################
# Visualize the distribution of cell type markers that were measured.
########################################
colors = c("#95e9f0", "#07f71b", "#d11111")

FeaturePlot(aurora, features = c("CD3", "CD4"), blend = TRUE, cols=colors) #CD4 T cells
FeaturePlot(aurora, features = c("CD3", "CD8"), blend = TRUE,cols=colors) #CD8 T cells
FeaturePlot(aurora, features = c("CD3", "gd"), blend = TRUE,cols=colors) #gd
FeaturePlot(aurora, features = c("CD3", "NK1.1"), blend = TRUE,cols=colors) #NKT cells
FeaturePlot(aurora, features = c("B220", "SiglecH"), blend = TRUE,cols=colors) #Plasmacytoid dendritic cells
FeaturePlot(aurora, features = c("MHCII", "CD11c"), blend = TRUE,cols=colors) #Dendritic cells
FeaturePlot(aurora, features = c("MHCII", "CD24"), blend = TRUE,cols=colors) # DC1
FeaturePlot(aurora, features = c("CD11c", "CD24"), blend = TRUE,cols=colors) # DC1
FeaturePlot(aurora, features = c("MHCII", "CD172"), blend = TRUE,cols=colors) #DC2
FeaturePlot(aurora, features = c("CD11c", "CD172"), blend = TRUE,cols=colors) #DC2
FeaturePlot(aurora, features = c("CD11b", "Ly6g"), blend = TRUE,cols=colors) #Neutrophils
FeaturePlot(aurora, features = c("CD11b", "Ly6c"), blend = TRUE,cols=colors) #Monocytes
FeaturePlot(aurora, features = c("CD11b", "F480"), blend = TRUE,cols=colors) #Macrophages

for (n in 0:23){
  img <-DimPlot(aurora, reduction = "umap", pt.size = 0.1, cells.highlight = WhichCells(aurora, idents = n), label = TRUE)
  print(img)
}

###Make violin plot for marker expression per cluster
variables = c("CD3","CD4","CD8","CD11b","NK1.1","B220","Ly6c","Ly6g","F480","SiglecH","MHCII","CD11c","CD24","CD172")
VlnPlot(aurora, features = variables, stack = T, flip = T)




########################################
# Visualize the distribution of ZsGreen intensity across clusters.
########################################
VlnPlot(aurora, features = "ZsGreen", flip = T, pt.size = 0, raster = F)

FeaturePlot(aurora, features = "ZsGreen", pt.size = 0.1, raster = F, cols = c("#e7ebe6", "#47de1d"), min.cutoff = 200, max.cutoff = 400) +
  NoLegend() + labs(title=NULL)
FeaturePlot(aurora, features = "ZsGreen", pt.size = 0.1, raster = F, cols = c("#e7ebe6", "#47de1d"), min.cutoff = 200, max.cutoff = 400)



########################################
# Visualize the distribution of DAPI intensity across clusters.
########################################
FeaturePlot(aurora, features = "DAPI", pt.size = 0.1, raster = F)







########################################
# Annotate cell based on cluster number
########################################

Idents(aurora, cells = WhichCells(aurora, idents = c("15","16","19","22","23"))) <- "Other"
Idents(aurora, cells = WhichCells(aurora, idents = c("2","5","14","20"))) <- "Macrophages"
Idents(aurora, cells = WhichCells(aurora, idents = c("0","3","4","17","18"))) <- "Monocytes"
Idents(aurora, cells = WhichCells(aurora, idents = c("10"))) <- "Neutrophils"
Idents(aurora, cells = WhichCells(aurora, idents = c("11"))) <- "pDC"
Idents(aurora, cells = WhichCells(aurora, idents = c("12"))) <- "DC1"
Idents(aurora, cells = WhichCells(aurora, idents = c("6"))) <- "DC2"
# Idents(aurora, cells = WhichCells(aurora, idents = c("5"))) <- "DC1 + DC2"
Idents(aurora, cells = WhichCells(aurora, idents = c("1"))) <- "NK"
# Idents(aurora, cells = WhichCells(aurora, idents = c("13"))) <- "NKT"
# Idents(aurora, cells = WhichCells(aurora, idents = c("12"))) <- "γδ"
Idents(aurora, cells = WhichCells(aurora, idents = c("9"))) <- "NKT + γδ"
Idents(aurora, cells = WhichCells(aurora, idents = c("13"))) <- "B Cells"
Idents(aurora, cells = WhichCells(aurora, idents = c("7"))) <- "CD8 T Cells"
Idents(aurora, cells = WhichCells(aurora, idents = c("8","21"))) <- "CD4 T Cells"






###fix the order of cell types
celltypes <- c("Other","Macrophages","Monocytes","Neutrophils","pDC","DC2","DC1","NK","NKT","γδ","B Cells","CD8 T Cells", "CD4 T Cells")
for (i in celltypes){
  Idents(aurora, cells = WhichCells(aurora, idents = i)) <- i
}

##########Confirm changes are made
levels(aurora)
# head(aurora)
table(Idents(aurora))
# Idents(aurora, cells = WhichCells(aurora, idents = c("γδ"))) <- "gd"
aurora$celltype <- Idents(aurora)
Idents(aurora) <- aurora$celltype
head(aurora)
DimPlot(aurora, reduction = "umap", pt.size = 0.1, raster = T)



##########Change the color of clusters. Manually enter hex code values in my_cols
require(scales)
# show_col(hue_pal()(12))
# my_cols <- c( 'CD4 T Cells'="#F8766D", 'CD8 T Cells'="#DE8C00", 'B Cells'="#B79F00", 'γδ'="#7CAE00",  'NKT'="#00BA38", 'NK'="#00C08B", 'DC1'="#00BFC4", 'DC2'="#00B4F0", 'pDC'="#619CFF", 'Neutrophils'="#C77CFF", 'Monocytes'="#F564E3", 'Macrophages'="#FF64B0", 'Other'="lightgrey")
my_cols <- c( 'CD4 T Cells'="#F8766D", 'CD8 T Cells'="#DE8C00", 'B Cells'="#B79F00", 'γδ'="#7CAE00",  'NKT'="#00BA38", 'NK'="#00C08B", 'DC1'="#00BFC4", 'DC2'="#00B4F0", 'pDC'="#619CFF", 'Neutrophils'="#C77CFF", 'Monocytes'="#F564E3", 'Macrophages'="#FF64B0", 'Other'="lightgrey")

head(aurora)
table(Idents(aurora))
table(aurora$celltype)
# head(aurora@assays$RNA@counts@Dimnames[[2]])
metadata <- data.frame(aurora@assays$RNA@counts@Dimnames[[2]], aurora$orig.ident, aurora$ZsGreen, aurora$DAPI, aurora@active.ident)
#export relevant metadata



##########Gives number of cells within each cluster.
table(aurora@active.ident)
data.frame(table(aurora@active.ident))
table(aurora$orig.ident)
subdata <- subset(aurora, subset = orig.ident == "B16ZS")
table(subdata$celltype)





##########UMAP plot without labels on clusters.
DimPlot(aurora, reduction = "umap", pt.size = 0.1, raster = F, cols = my_cols, label = F, repel = F)

DimPlot(aurora, reduction = "umap", pt.size = 0.1, raster = F, cols = my_cols, label = F) + NoLegend()
# DimPlot(aurora, reduction = "umap", pt.size = 1, raster = T, cols = my_cols, label = F) + NoLegend()

aurora.downsample = subset(aurora, cells = sample(Cells(aurora), 25000))
DimPlot(aurora.downsample, reduction = "umap", pt.size = 0.1, raster = F, cols = my_cols, label = F) + NoLegend()
FeaturePlot(aurora.downsample, features = "ZsGreen", pt.size = 0.1, raster = F, cols = c("#e7ebe6", "#47de1d"), min.cutoff = 200, max.cutoff = 400) +
  NoLegend() + labs(title=NULL)
FeaturePlot(aurora.downsample, features = "ZsGreen", pt.size = 0.1, raster = F, cols = c("#e7ebe6", "#47de1d"), min.cutoff = 200, max.cutoff = 400)




##########UMAP plot with labels on clusters.
DimPlot(aurora, reduction = "umap", pt.size = 0.1, raster = F, cols = my_cols, label = T, repel = T)





##########Violin plots by cluster for marker variables, and ZsGreen.
VlnPlot(aurora, features = variables, stack = T, flip = T, fill.by = "ident", cols=my_cols)
DotPlot(aurora, features = variables, assay = "RNA") + RotatedAxis() + coord_flip()

VlnPlot(aurora, features = "ZsGreen", flip = T, pt.size = 0, raster = F, cols=my_cols)

max(aurora$ZsGreen)

test <- subset(aurora, subset = orig.ident == "B16")
plots <- VlnPlot(test, features = "ZsGreen", pt.size = 0, raster = F, cols=my_cols, y.max = 1023)
plots


data$types <- test$celltype










