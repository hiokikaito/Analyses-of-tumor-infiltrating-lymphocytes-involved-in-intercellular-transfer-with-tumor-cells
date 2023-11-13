#################################################
# Name: scRNAseq_CD45+_TILs.R
# 
# Description: R codes used for the analysis of scRNAseq data on all CD45+ sorted Tumor Infiltrating Lymphocytes collected from B16 ZsGreen+ tumors.
#
# Input: scRNAseq data in 10x Genomics file format.
# 
# Date updated: Nov. 10, 2023
# 
# Author: Kaito Hioki 
# Contact: khioki@umass.edu
#################################################


library(ggplot2)
#install.packages("dplyr")
library(dplyr)
# install.packages('Seurat')
library(Seurat)
# install.packages('patchwork')
library(patchwork)

library(RColorBrewer)


###Obtain input file and create Seurat object with raw data.
zsme.data <- Read10X(data.dir = "10x-file_directory/")

#Initialize the Seurat object with the raw (non-normalized data).
zsme <- CreateSeuratObject(counts = zsme.data, project = "ZSME", min.cells = 3, min.features = 0, names.field = 2, names.delim = "-")





########################
#Quality control
########################
###calculates the percentage of counts originating from Mitochondrial genes
zsme[["percent.mt"]] <- PercentageFeatureSet(zsme, pattern = "^mt-")

# Visualize QC metrics as a violin plot
VlnPlot(zsme, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

plot1 <- FeatureScatter(zsme, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(zsme, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

###Filter out cells with less than 1000 features and more than 6000 features, or cells with >10% mitochondrial reads
#From the plots, most of cells are between those thresholds.
zsme <- subset(zsme, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 10)







########################
#Normalizing data
########################
zsme <- NormalizeData(zsme, normalization.method = "LogNormalize", scale.factor = 10000)



########################
#Identify highly variable features
########################
zsme <- FindVariableFeatures(zsme, selection.method = "vst", nfeatures = 2000)



########################
#Scaling the data
########################
all.genes <- rownames(zsme)
zsme <- ScaleData(zsme, features = all.genes)




########################
#Perform linear dimensional reduction; PCA
########################
zsme <- RunPCA(zsme, features = VariableFeatures(object = zsme))




########################
#Determine the 'dimensionality' of the dataset
########################
#distribution of p-values for each PC. significant PCs with low p-values look like a solid curve above the dashed line.
zsme <- JackStraw(zsme, num.replicate = 100, dims = 50)
JackStrawPlot(zsme, dims = 1:50)

#ElbowPlot ranks PCs based on percentage of variance explained.
zsme <- ScoreJackStraw(zsme, dims = 1:50)
ElbowPlot(zsme, ndims=50)




########################
#Cluster the cells
########################
zsme <- FindNeighbors(zsme, dims = 1:13)
zsme <- FindClusters(zsme, resolution = 1)





########################
#Run non-linear dimensional reduction 
########################
zsme <- RunUMAP(zsme, dims = 1:13)
DimPlot(zsme, reduction = "umap", pt.size = 1, label = T, raster = T)



###Got 17 clusters (0-16)
#check the location of individual clusters
for (num in seq(0:16)){
  plot <- DimPlot(zsme, reduction = "umap", pt.size = 0.1, cells.highlight = WhichCells(zsme, idents = num), label = TRUE) + NoLegend()
  print(plot)
}  







##############################
###Identify cell types per cluster
##############################

cluster.markers <- FindMarkers(zsme, ident.1 = 12, min.pct = 0.25)



###distribution of cell type marker genes
T_markers <- c("Cd3e", "Cd3d", "Cd3g", "Cd4", "Cd8a", "Cd44")
B_markers <- c("Cd19","Ms4a1", )
gd_markers <-c("Maf","Ccr2","Tmem176a","Zbtb16","Icos")
NKT_markers <- c("Gzma","Gzmb","Klra9","Klra5","Klrc1","Klre1")
NK_markers <- c("Eomes", "Tbx21")
ILC_markers <- c("Kit","Ncr1")

DC1_markers <- c("Xcr1", "Cd24a", "Clec9a","Itgae","Itgax","Ly75","Thbd","Btla","Cadm1","Cd8a","Gpr141","Naaa","Cd36","Rab7b","Irf8")
DC2_markers <- c("Clec4a4","Sirpa","Cd14","Adgre1","Csf1r","Cd163","Clec10a","Notch2","Itgam","Cx3cr1","Cd2","Cd207","Mgl2","Cd209","Mmp12","Cybb","Irf4","Zeb2")


pDC_markers <- c("Spib", "Ccr9", "Siglech", "Batf3", "Itgae")
Macr_markers <- c("Mertk","Fcgr1","Lyz2","Itgam","Itgax")
Neut_markers <- c("Fcgr3","Cd33","Cxcr2","Cxcr4")

markers <- c(T_markers,B_markers,gd_markers,NKT_markers,NK_markers,ILC_markers,DC_markers,pDC_markers,Macr_markers)

# for (gene in T_markers){
#   plot <- FeaturePlot(zsme, features = gene)
#   print(plot)
# }




###########
#Dot plot of marker genes
###########
T_markers <- c("Cd3e", "Cd3d", "Cd3g", "Cd8a", "Cd4", "Cd44")
B_markers <- c("Cd19","Ms4a1","Cd79a","Cd79b")
gd_markers <-c("Maf","Ccr2","Tmem176a","Zbtb16","Icos")
NKT_markers <- c("Gzma","Gzmb","Klra9","Klra5","Klrc1","Klre1")
NK_markers <- c("Eomes", "Tbx21","Itga2","Klrb1c","Ncr1")
ILC_markers <- c("Kit","Ncr1")
DC_markers <- c("Id2","Zbtb46","Xcr1","Cd24a","Clec9a","Tlr3","Ly75","Thbd","Btla","Cadm1","Irf8","Sirpa","Cd14","Clec10a","Cd2","Irf4","Klf4")
DC1_markers <- c("Xcr1", "Cd24a", "Clec9a")
DC2_markers <- c("Mgl2","Sirpa","Cd14")
pDC_markers <- c("Spib", "Ccr9", "Siglech")
Macr_markers <- c("Csf1r","Mertk","Fcgr1","Lyz2","Itgam")

Neut_markers <- c("Fcgr3","Cd33","Cxcr2","Cxcr4")

markers <- c(T_markers,B_markers,NK_markers,DC1_markers,DC2_markers,pDC_markers,Macr_markers)




###show bubble plots for selected marker genes
DotPlot(zsme, features = markers, assay = "RNA", dot.scale = 10) + 
  scale_color_gradient2(low="#4575b4", mid = "#ffffbf",  high="#d73027") +
  RotatedAxis()


###show average expression of variable marker genes
test.markers <- FindAllMarkers(zsme, only.pos = T, min.pct = .25, logfc.threshold = 0.25)
test.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC) -> topmarkers
DoHeatmap(zsme, features = topmarkers$gene)

test2 <- AverageExpression(zsme, features = topmarkers$gene)
test2 <- data.frame(test2$RNA)
library(pheatmap)
cal_z_score <- function(x){
  # x <- log(x*1000)
  (x - mean(x)) / sd(x)
}

data_norm <- t(apply(test2, 1, cal_z_score))
data_norm[data_norm == "NaN"] <- 0
pheatmap(data_norm, 
         color = colorRampPalette(c("purple", "black", "yellow"))(100),
         cluster_cols = F, 
         cluster_rows = F)





########################
###Assemble plot 
###Try to recluster
########################
#Lymphocytes are clusters 0,2,3,6,12
#Myeloid are 1,4,5,7,8,9,10,11
head(zsme)
Idents(zsme) <- zsme$RNA_snn_res.1
table(Idents(zsme))
Idents(zsme, cells = WhichCells(zsme, idents = c(16))) <- "B cells"
Idents(zsme, cells = WhichCells(zsme, idents = c(12))) <- "pDC"
Idents(zsme, cells = WhichCells(zsme, idents = c(3,9,8))) <- "T_cells"
Idents(zsme, cells = WhichCells(zsme, idents = c(1,2,10,14))) <- "NK_cells"
Idents(zsme, cells = WhichCells(zsme, idents = c(7,11,15))) <- "cDCs"
Idents(zsme, cells = WhichCells(zsme, idents = c(0,4,5,6,13))) <- "Macr"

zsme$first_clusters <- zsme$RNA_snn_res.1
zsme$celltype <- Idents(zsme)
table(zsme$celltype)

###updating cell type annotations to full object
table(Idents(zsme))
zsme$group <- Idents(zsme)
table(zsme$group)
table(zsme$celltype)
Idents(zsme) <- zsme$celltype
DimPlot(zsme, reduction = "umap", pt.size = 1, label = T, raster = T)




######################
#subset macrophage cell types to recluster
######################
sub <- subset(zsme, subset = group == "Macr")
head(sub)


###recluster cells in subset 
#figure out dimensions
sub <- JackStraw(sub, num.replicate = 100, dims = 50)
sub <- ScoreJackStraw(sub, dims = 1:50)
JackStrawPlot(sub, dims = 1:50)
ElbowPlot(sub, ndims=50)

#chose dimensions for this cluster
sub <- FindNeighbors(sub, dims = 1:13)
sub <- FindClusters(sub, resolution = 0.5)

sub <- RunUMAP(sub, dims = 1:13)


DimPlot(sub, reduction = "umap", pt.size = 1, label = T, raster = F)
DimPlot(subset(sub, subset = orig.ident == 1), reduction = "umap", pt.size = 1)
DimPlot(subset(sub, subset = orig.ident == 2), reduction = "umap", pt.size = 1)

###Distribution of marker genes
Macr_markers <- c("Mertk","Fcgr1","Lyz2","Itgam","Itgax")
Macro_markers <- c("Cx3cr1", "Nos2", "Xcr1", "Flt3", "Ccr7", "Siglech", "S100a8", "Trgv2", "Mki67", "Cd79a", "Ncr1", "Ptprc")
M1_markers <- c("Cd80","Cd86","Fcgr1","Fcgr2b", "Stat1","Nos2","Il1b","Il6","Il12b","Ccr7","Inhba","Tnf")
M2_markers <- c("Mrc1","Cd163","Msr1","Vegfa","Maf","Arg1","Chil3","Retnla","Egr2","Fn1")

FeaturePlot(sub, features = M2_markers)
# FeaturePlot(sub, features = c("Folr2", "Trem2"))


paper_markerst <- c("Chil1","Mertk","Vcam1","Mrc1","Ccl7","C1qa","Cd81","Pf4","Cd63",
                    "Cx3cr1","Cd72","Fgd2","Ccl2","Mki67","Itga4","Ccr2","Cd83","Il1r2",
                    "Cxcl9","Itgax","Sell","Plac8","Ly6c2","Fas","Rsad2","Isg15","Ifit2",
                    "Ifit3","Ms4a4c","Cd9","Cxcl2","Cd1d1","Nos2","Arg1","Inhba","Il7r",
                    "Hcar2","Pdcd1lg2","Clec4e","Ly6i","Cd274","Sgk1","Nov")
DoHeatmap(sub, features = paper_markerst, size = 1)


MonoMacs_markers <- c("Cx3cr1","Ly6i","Ly6c1", "Ly6c2","Ccr2","Itga4","Sell","Mertk","Cd63","Cd9", "Itgam","Adgre1")

FeaturePlot(sub, features = MonoMacs_markers)
for (marker in MonoMacs_markers){
  plot <- FeaturePlot(sub, features = marker, pt.size = 2) + NoLegend()
  print(plot)
}



phagocyte_markers <- c("Cd14","Mrc1","Cd163","Cd209a")
FeaturePlot(sub, features = phagocyte_markers)



###unbiased way to find subcluster markers
sub.markers <- FindAllMarkers(sub, only.pos = T, min.pct = .25, logfc.threshold = 0)
# write.csv(sub.markers, "G:/.shortcut-targets-by-id/1d1S_Tpnk3gFjNVrvVi2QWzl0o0BtDwzP/KH/KH-22/DEG_analysis/Macrophages/Mac_variable_genes_upregulated-only.csv")
sub.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC) -> topmarkers
topmarkers <- topmarkers$gene
DoHeatmap(sub, features = topmarkers)
DoHeatmap(sub, features = "Cx3cr1")


#get markers of one cluster compared to all other clusters
test <- FindMarkers(sub, ident.1 = 5,only.pos = T, min.pct = .25, logfc.threshold = 0.25)


# topmarkers <- append(topmarkers, MonoMacs_markers)
topmarkers <- MonoMacs_markers

DotPlot(sub, features = topmarkers, assay = "RNA", dot.scale = 10) + 
  scale_color_gradient2(low="#4575b4", mid = "#ffffbf",  high="#d73027") +
  RotatedAxis()





#######make a plot of signature genes
signature_plot <- function (seuratobject,signature_genes){
  sig_data <- as.data.frame(Embeddings(object = seuratobject, reduction = "umap"))
  sig_data$signature <- rowMeans(t(seuratobject@assays$RNA@scale.data[signature_genes,]))
  ggplot(sig_data, aes (UMAP_1, UMAP_2))+
    geom_point(aes (color = signature))+
    scale_color_gradientn(colors = topo.colors(8))+
    theme_classic() +
    labs(title = NULL) +
    NoLegend()
}
signature_plot_2 <- function (seuratobject,signature_genes){
  sig_data <- as.data.frame(Embeddings(object = seuratobject, reduction = "umap"))
  sig_data$signature <- rowMeans(t(seuratobject@assays$RNA@scale.data[signature_genes,]))
  ggplot(sig_data, aes (UMAP_1, UMAP_2))+
    geom_point(aes (color = signature))+
    scale_color_gradientn(colors = topo.colors(8))+
    theme_classic() +
    labs(title = NULL)
}


#monocyte genes
# signature_genes <- c("Cx3cr1","Ly6i","Ly6c1","Ly6c2","Il1r2","Ccr2","Itga4","Sell")
signature_genes <- c("Ly6i","Ly6c1","Ly6c2","Ccr2","Itga4","Sell")
#macrophage genes
signature_genes <- c("Mertk","Cd63","Cd9", "Itgam", "Adgre1")

signature_plot(sub, signature_genes)
signature_plot_2(sub, signature_genes)


#individual genes
for (gene in signature_genes){
  plot <- FeaturePlot(sub, features = gene, pt.size = 2) + labs(title = NULL) + NoLegend()
  print(plot)
}

FeaturePlot(sub, features = "Cx3cr1", pt.size = 2) + labs(title = NULL) + NoLegend()





##############
#Add id of subset macrophages to zsme object
######################
table(Idents(zsme))
head(zsme)
table(zsme$celltype)
Idents(zsme) <- zsme$celltype


Idents(zsme, cells = WhichCells(zsme, idents = c("Macr"))) <- Idents(sub)
zsme$celltype <- Idents(zsme)



#############
#Calculate DEGs
#####################
mac_types <- c("Antiviral_Macs", "Circulating_Macs", "Inflammatory_Macs", "Nos2_Macs", "Cx3cr1_Macs")

head(zsme)
sub2 <- subset(zsme, subset = idents)
head(sub2)
table(Idents(sub2))
Idents(sub2) <- sub2$orig.ident
table(Idents(sub2))

DEGs <- FindMarkers(sub2, ident.1 = 1, ident.2 = 2, min.pct = 0, logfc.threshold = 0, test.use = "wilcox")





######################
#calculate DEGs of all ZsHigh vs ZsLow macrophages 
######################
zsme$group <- Idents(zsme)
table(Idents(zsme))
table(zsme$group)
table(zsme$celltype)

Idents(zsme, cells = WhichCells(zsme, idents = c("Cx3cr1_Macs"))) <- "Macrophages"
zsme$group <- Idents(zsme)

sub2 <- subset(zsme, subset = group == "Macrophages")
Idents(sub2) <- sub2$orig.ident
DEGs <- FindMarkers(sub2, ident.1 = 1, ident.2 = 2, min.pct = 0, logfc.threshold = 0, test.use = "wilcox")




######################
#calculate DEGs of all ZsHigh vs ZsLow B cells 
######################


#B cells were not analyzed, since there were not enough cells in the ZsHigh group to compare to the ZsLow group.








######################
#subset NK cell types to recluster
######################
sub <- subset(zsme, subset = celltype == "NK cells")
head(sub)

###recluster cells in subset 
#figure out dimensions
sub <- JackStraw(sub, num.replicate = 100, dims = 50)
sub <- ScoreJackStraw(sub, dims = 1:50)
JackStrawPlot(sub, dims = 1:50)
ElbowPlot(sub, ndims=50)

#chose dimensions for this cluster
sub <- FindNeighbors(sub, dims = 1:8)
sub <- FindClusters(sub, resolution = 0.5)

sub <- RunUMAP(sub, dims = 1:8)

DimPlot(sub, reduction = "umap", pt.size = 1, label = T, raster = T)
DimPlot(subset(sub, subset = orig.ident == 1), reduction = "umap", pt.size = 1)
DimPlot(subset(sub, subset = orig.ident == 2), reduction = "umap", pt.size = 1)

table(sub$orig.ident)
table(Idents(sub))
table(Idents(subset(sub, subset = orig.ident == 1)))
table(Idents(subset(sub, subset = orig.ident == 2)))


DimPlot(sub, reduction = "umap", pt.size = 1, label = T, raster = T)


sub.markers <- FindAllMarkers(sub, only.pos = T, min.pct = .25, logfc.threshold = 0.25)
sub.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC) -> topmarkers
topmarkers <- topmarkers$gene
DoHeatmap(sub, features = topmarkers)


Idents(sub) <- sub$orig.ident
DEGs <- FindMarkers(sub, ident.1 = 1, ident.2 = 2, min.pct = 0, logfc.threshold = 0, test.use = "wilcox")


NK_markers <- c("Eomes", "Tbx21", "Itga2", "Klrb1c", "Ncr1", "Itgam", "Cd27")

FeaturePlot(sub, features = NK_markers)

VlnPlot(sub, features = NK_markers, stack = T, flip = T, fill.by = "ident", assay = "RNA")





#####################
# DEGs for pDC
#######################
sub <- subset(zsme, subset = celltype == "pDC")
table(Idents(sub))
table(Idents(subset(sub, subset = orig.ident == 1)))
table(Idents(subset(sub, subset = orig.ident == 2)))
Idents(sub) <- sub$orig.ident
DDEGs <- FindMarkers(sub, ident.1 = 1, ident.2 = 2, min.pct = 0, logfc.threshold = 0, test.use = "wilcox")










#####################
# Recluster cDC for DC1 vs DC2, and pDC
#######################
table(Idents(zsme))
head(zsme)
Idents(zsme, cells = WhichCells(zsme, idents = c("pDC"))) <- "DC"
Idents(zsme, cells = WhichCells(zsme, idents = c("DC1"))) <- "DC"
Idents(zsme, cells = WhichCells(zsme, idents = c("DC2"))) <- "DC"
zsme$celltype <- Idents(zsme)


sub <- subset(zsme, subset = celltype == "DC")
table(sub$orig.ident)
table(Idents(sub))
table(Idents(subset(sub, subset = orig.ident == 1)))
table(Idents(subset(sub, subset = orig.ident == 2)))

###recluster cells in subset 
#figure out dimensions
sub <- JackStraw(sub, num.replicate = 100, dims = 50)
sub <- ScoreJackStraw(sub, dims = 1:50)
JackStrawPlot(sub, dims = 1:50)
ElbowPlot(sub, ndims=50)

#chose dimensions for this cluster
sub <- FindNeighbors(sub, dims = 1:10)
sub <- FindClusters(sub, resolution = 1)

sub <- RunUMAP(sub, dims = 1:10)

DimPlot(sub, reduction = "umap", pt.size = 1, label = T, raster = T)
DimPlot(subset(sub, subset = orig.ident == 1), reduction = "umap", pt.size = 1)
DimPlot(subset(sub, subset = orig.ident == 2), reduction = "umap", pt.size = 1)

table(sub$orig.ident)
table(Idents(sub))
table(Idents(subset(sub, subset = orig.ident == 1)))
table(Idents(subset(sub, subset = orig.ident == 2)))



####evaluate markers
DC1_markers <- c("Xcr1", "Cd24a", "Clec9a","Itgae","Itgax","Ly75","Thbd","Btla","Cadm1","Cd8a","Gpr141","Naaa","Cd36","Rab7b","Irf8")
DC2_markers <- c("Clec4a4","Sirpa","Cd14","Adgre1","Csf1r","Cd163","Clec10a","Notch2","Itgam","Cx3cr1","Cd2","Cd207","Mgl2","Cd209","Mmp12","Cybb","Irf4","Zeb2")
pDC_markers <- c("Spib", "Ccr9", "Siglech", "Batf3", "Itgae")

FeaturePlot(sub, features = DC1_markers)
FeaturePlot(sub, features = DC2_markers)
FeaturePlot(sub, features = pDC_markers)

sub.markers <- FindAllMarkers(sub, only.pos = T, min.pct = .25, logfc.threshold = 0.25)
sub.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(sub, features = top10$gene)


########based on markers:
Idents(sub, WhichCells(sub, idents = c(1,2,6,7))) <- "DC1"
Idents(sub, WhichCells(sub, idents = c(0,4,8,9))) <- "DC2"
Idents(sub, WhichCells(sub, idents = c(3,5,10))) <- "pDC"
sub$celltype <- Idents(sub)

DimPlot(sub, reduction = "umap", pt.size = 1, label = T, raster = T)
DimPlot(subset(sub, subset = orig.ident == 1), reduction = "umap", pt.size = 1)
DimPlot(subset(sub, subset = orig.ident == 2), reduction = "umap", pt.size = 1)

table(sub$orig.ident)
table(Idents(sub))
table(Idents(subset(sub, subset = orig.ident == 1)))
table(Idents(subset(sub, subset = orig.ident == 2)))

##############DEGs for DCs
sub2 <- subset(sub, subset = celltype == "pDC")
table(sub2$orig.ident)
table(Idents(sub2))
table(Idents(subset(sub2, subset = orig.ident == 1)))
table(Idents(subset(sub2, subset = orig.ident == 2)))
Idents(sub2) <- sub2$orig.ident
DEGs <- FindMarkers(sub2, ident.1 = 1, ident.2 = 2, min.pct = 0, logfc.threshold = 0, test.use = "wilcox")



Idents(zsme, WhichCells(sub, idents = "DC1")) <- "DC1"
Idents(zsme, WhichCells(sub, idents = "DC2")) <- "DC2"
Idents(zsme, WhichCells(sub, idents = "pDC")) <- "pDC"

table(Idents(zsme))
DimPlot(zsme)









#######################################
# Combine T cells (CD4, CD8, NKT)
#######################################
table(Idents(zsme))
table(zsme$celltype)
table(zsme$group)

sub3 <- zsme
Idents(sub3) <- sub3$celltype

Idents(zsme, WhichCells(sub3, idents = "T_cells")) <- "T cells"

###get DEGs
zsme$celltype <- Idents(zsme)
sub2 <- subset(zsme, subset = celltype == "T cells")
table(sub2$orig.ident)
table(Idents(sub2))
table(Idents(subset(sub2, subset = orig.ident == 1)))
table(Idents(subset(sub2, subset = orig.ident == 2)))
Idents(sub2) <- sub2$orig.ident
DEGs <- FindMarkers(sub2, ident.1 = 1, ident.2 = 2, min.pct = 0, logfc.threshold = 0, test.use = "wilcox")












############################
#Finish annotations, save rds file
##############################
DimPlot(zsme)


celltypes <- c("Macrophages","pDC","DC2","DC1","NK cells","B cells","T cells")
for (i in celltypes){
  Idents(zsme, cells = WhichCells(zsme, idents = i)) <- i
}

DimPlot(zsme, reduction = "umap", label = T, raster = T) 
DimPlot(subset(zsme, subset = orig.ident == 1), reduction = "umap", pt.size = 1)
DimPlot(subset(zsme, subset = orig.ident == 2), reduction = "umap", pt.size = 1)


table(zsme$orig.ident)
table(Idents(zsme))
table(Idents(subset(zsme, subset = orig.ident == 1)))
table(Idents(subset(zsme, subset = orig.ident == 2)))

# saveRDS(zsme, file = "annotated_CD45+_TILs.rds")







##################################
#Make average expression heatmaps for genes selected from GSEA (ones that contributed most to enrichment of pathway terms)
##################################
head(zsme)
table(Idents(zsme))
table(zsme$celltype)

celltypes <-  c("T cells", "B cells", "NK cells", "DC1", "DC2", "pDC", "Macrophages")
ID <- c("T", "B", "NK", "DC1", "DC2", "pDC", "Macs")
for (i in 1:7){
  high <- paste("ZsHigh_", ID[i])
  low <- paste("ZsLow_", ID[i])
  
  sub <- subset(zsme, subset = celltype == celltypes[i])
  sub$cellgroup <- ifelse(sub$orig.ident == 2, high, low)
  # table(sub$cellgroup)
  Idents(zsme, WhichCells(sub)) <- sub$cellgroup
}




library(pheatmap)
cal_z_score <- function(x){
  # x <- log(x*1000)
  (x - mean(x)) / sd(x)
}

GSEA_genes <- read.csv("genelist.csv")


table(Idents(zsme))

test <- AverageExpression(object = zsme)
average_data <- data.frame(test$RNA)

sample <- average_data
sample$gene <- rownames(average_data)

sample2 <- sample %>% 
  filter(gene %in% GSEA_genes$WNT_beta_catenin)
sample2 <- sample2[,-(ncol(sample2))]

data_norm <- t(apply(sample2, 1, cal_z_score))
data_norm[data_norm == "NaN"] <- 0
pheatmap(data_norm, 
         color = colorRampPalette(c("#4575b4", "#ffffbf", "#d73027"))(100),
         cluster_cols = F, 
         cluster_rows = F,
         cellwidth = 15, cellheight = 15)





##############################
# Extract the myeloid cells only and  show UMAP
##############################
table(Idents(zsme))
zsme$celltype <- Idents(zsme)
table(zsme$celltype)
zsme$group <- Idents(zsme)
table(zsme$group)

Idents(zsme, WhichCells(zsme, idents = "DC1")) <- "myel"
Idents(zsme, WhichCells(zsme, idents = "DC2")) <- "myel"
Idents(zsme, WhichCells(zsme, idents = "pDC")) <- "myel"
Idents(zsme, WhichCells(zsme, idents = "Macrophages")) <- "myel"

zsme$group <- Idents(zsme)

Idents(zsme) <- zsme$celltype


myeloids <- subset(zsme, subset = group == "myel")

DimPlot(myeloids)

###recluster cells in subset 
#figure out dimensions
myeloids <- JackStraw(myeloids, num.replicate = 100, dims = 50)
myeloids <- ScoreJackStraw(myeloids, dims = 1:50)
JackStrawPlot(myeloids, dims = 1:50)
ElbowPlot(myeloids, ndims=50)

#chose dimensions for this cluster
myeloids <- FindNeighbors(myeloids, dims = 1:13)
myeloids <- FindClusters(myeloids, resolution = 0.5)

myeloids <- RunUMAP(myeloids, dims = 1:13)

DimPlot(myeloids, reduction = "umap", pt.size = 1, label = T, raster = T)
DimPlot(subset(myeloids, subset = orig.ident == 1), reduction = "umap", pt.size = 1)
DimPlot(subset(myeloids, subset = orig.ident == 2), reduction = "umap", pt.size = 1)

Idents(myeloids) <- myeloids$celltype



