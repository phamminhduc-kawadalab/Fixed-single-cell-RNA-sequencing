library(Seurat)
library(dplyr)
library(patchwork)
library(tidyverse)

markerGenes <- c("Myh11","Actg22", #VSMC
                 "Rgs5", #Pericytes
                 "Lrat","Dcn","Colec11","Pdgfra","Pdgfrb", #HSC
                 "Acta2","Lox","Col3a1","Col1a1", #Act-HSC
                 "Mgp","Eln","Mfap4","Dpt","Gsn", #PF markers
                 "Kdr","Aqp1", #Endothelial
                 "Top2a","Mki67", #Cycling
                 'Upk3b','Msln', #Mesothelial
                 'Epcam', 'Krt7', #Cholagiocytes
                 "Il7r","Cd3d","Trac", #T-cells
                 "Nkg7","Gzma", #NK-T
                 "Clec4f","Visg4",#Kupffer cells
                 "Lgals3","Ccr2", "Lyz2", #Macrophages
                 "Kldr1","Mycl","Cd74","Itgax", # Monocytes
                 "Cxcr2", "S100a9","Ngp", #Neutorphils
                 "Bpgm","Hba-a1", #Erythrocyte
                 "Cd79a","Ebf1","Ms4a1", #B-cell
                 "Jchain",          #Plasma
                 "Siglech","Cox6a2", #Siglech
                 "Alb","Serpina1a", #Hepatocytes
                 "Apoa1","Serpina3k",
                 "Abcc2","Sema4g")
#1, WT1198

wt1198.data <- Read10X_h5("/Users/lethithanhthuy/Desktop/20230410_Single\ fixed\ cell/Fixed\ cell\ TAA/filtered_feature_bc_matrix_WT1198_Normal.h5")
dim(wt1198.data)
sampleName <-  'wt1198'
dataset1 <- CreateSeuratObject(counts = wt1198.data, project = sampleName, min.cells = 3, min.features = 200)
dataset1
dataset1[["percent.mt"]] <- PercentageFeatureSet(dataset1, pattern = "^mt-")
dataset1[["percent.ribo"]] <- PercentageFeatureSet(dataset1, pattern = "^Rpl|^Rps|^Mrps|^Mrpl")
dataset1<- subset(dataset1, subset = percent.mt < 20 & nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 20000 )
VlnPlot(dataset1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4)
FeatureScatter(dataset1, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(dataset1, feature1 = "nFeature_RNA", feature2 = "percent.mt")
FeatureScatter(dataset1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dataset1 <- NormalizeData(object = dataset1)
dataset1 <- FindVariableFeatures(object = dataset1, nfeatures = 3000)
dataset1 <- ScaleData(object = dataset1)
dataset1 <- RunPCA(object = dataset1)
dataset1 <- FindNeighbors(object = dataset1, dims = 1:30)
dataset1 <- RunUMAP(object = dataset1, reduction = "pca", dims = 1:30)
dataset1 <- FindClusters(object = dataset1)
DimPlot(dataset1, label = T)

FeaturePlot(object = dataset1, features = c("percent.mt"))
FeaturePlot(object = dataset1, features = c("percent.ribo"))
FeaturePlot(object = dataset1, features = c("nCount_RNA"))
FeaturePlot(object = dataset1, features = c("nFeature_RNA"))
DimPlot(object = dataset1, label = T)
DotPlot(dataset1, features = markerGenes,cluster.idents = T) + RotatedAxis()
VlnPlot(dataset1, features = "PECAM1", pt.size = 0.25)

#, wt1199

wt1199.data <- Read10X_h5("/Users/lethithanhthuy/Desktop/20230410_Single\ fixed\ cell/Fixed\ cell\ TAA/filtered_feature_bc_matrix_WT1199_Normal.h5")
dim(wt1199.data)
sampleName <-  'wt1199'
dataset2 <- CreateSeuratObject(counts = wt1199.data, project = sampleName, min.cells = 3, min.features = 200)
dataset2
dataset2[["percent.mt"]] <- PercentageFeatureSet(dataset2, pattern = "^mt-")
dataset2[["percent.ribo"]] <- PercentageFeatureSet(dataset2, pattern = "^Rpl|^Rps|^Mrps|^Mrpl")
dataset2<- subset(dataset2, subset = percent.mt < 20 & nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 20000 )
VlnPlot(dataset2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4)
FeatureScatter(dataset2, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(dataset2, feature1 = "nFeature_RNA", feature2 = "percent.mt")
FeatureScatter(dataset2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dataset2 <- NormalizeData(object = dataset2)
dataset2 <- FindVariableFeatures(object = dataset2, nfeatures = 3000)
dataset2 <- ScaleData(object = dataset2)
dataset2 <- RunPCA(object = dataset2)
dataset2 <- FindNeighbors(object = dataset2, dims = 1:30)
dataset2 <- RunUMAP(object = dataset2, reduction = "pca", dims = 1:30)
dataset2 <- FindClusters(object = dataset2)

FeaturePlot(object = dataset2, features = c("percent.mt"))
FeaturePlot(object = dataset2, features = c("percent.ribo"))
FeaturePlot(object = dataset2, features = c("nCount_RNA"))
FeaturePlot(object = dataset2, features = c("nFeature_RNA"))
DimPlot(object = dataset2, label = T)
DotPlot(dataset2, features = markerGenes,cluster.idents = T) + RotatedAxis()
VlnPlot(dataset2, features = "MRC1", pt.size = 0.25)

#, wt1200

wt1200.data <- Read10X_h5("/Users/lethithanhthuy/Desktop/20230410_Single\ fixed\ cell/Fixed\ cell\ TAA/filtered_feature_bc_matrix_WT1200_Normal.h5")
dim(wt1200.data)
sampleName <-  'wt1200'
dataset3 <- CreateSeuratObject(counts = wt1200.data, project = sampleName, min.cells = 3, min.features = 200)
dataset3
dataset3[["percent.mt"]] <- PercentageFeatureSet(dataset3, pattern = "^mt-")
dataset3[["percent.ribo"]] <- PercentageFeatureSet(dataset3, pattern = "^Rpl|^Rps|^Mrps|^Mrpl")
dataset3<- subset(dataset3, subset = percent.mt < 20 & nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 20000 )
VlnPlot(dataset3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4)
FeatureScatter(dataset3, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(dataset3, feature1 = "nFeature_RNA", feature2 = "percent.mt")
FeatureScatter(dataset3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dataset3 <- NormalizeData(object = dataset3)
dataset3 <- FindVariableFeatures(object = dataset3, nfeatures = 3000)
dataset3 <- ScaleData(object = dataset3)
dataset3 <- RunPCA(object = dataset3)
dataset3 <- FindNeighbors(object = dataset3, dims = 1:30)
dataset3 <- RunUMAP(object = dataset3, reduction = "pca", dims = 1:30)
dataset3 <- FindClusters(object = dataset3)

FeaturePlot(object = dataset3, features = c("percent.mt"))
FeaturePlot(object = dataset3, features = c("percent.ribo"))
FeaturePlot(object = dataset3, features = c("nCount_RNA"))
FeaturePlot(object = dataset3, features = c("nFeature_RNA"))
DimPlot(object = dataset3, label = T)
DotPlot(dataset3, features = markerGenes,cluster.idents = T) + RotatedAxis()
VlnPlot(dataset3, features = "LYVE1", pt.size = 0.25)

#1, WT11

wt11.data <- Read10X_h5("/Users/lethithanhthuy/Desktop/20230410_Single\ fixed\ cell/Fixed\ cell\ TAA/filtered_feature_bc_matrix_WT11-TAA10wk.h5")
dim(wt11.data)
sampleName <-  'wt11'
dataset1_TAA <- CreateSeuratObject(counts = wt11.data, project = sampleName, min.cells = 3, min.features = 200)
dataset1_TAA
dataset1_TAA[["percent.mt"]] <- PercentageFeatureSet(dataset1_TAA, pattern = "^mt-")
dataset1_TAA[["percent.ribo"]] <- PercentageFeatureSet(dataset1_TAA, pattern = "^Rpl|^Rps|^Mrps|^Mrpl")
dataset1_TAA<- subset(dataset1_TAA, subset = percent.mt < 20 & nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 20000 )
VlnPlot(dataset1_TAA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4)
FeatureScatter(dataset1_TAA, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(dataset1_TAA, feature1 = "nFeature_RNA", feature2 = "percent.mt")
FeatureScatter(dataset1_TAA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dataset1_TAA <- NormalizeData(object = dataset1_TAA)
dataset1_TAA <- FindVariableFeatures(object = dataset1_TAA, nfeatures = 3000)
dataset1_TAA <- ScaleData(object = dataset1_TAA)
dataset1_TAA <- RunPCA(object = dataset1_TAA)
dataset1_TAA <- FindNeighbors(object = dataset1_TAA, dims = 1:30)
dataset1_TAA <- RunUMAP(object = dataset1_TAA, reduction = "pca", dims = 1:30)
dataset1_TAA <- FindClusters(object = dataset1_TAA)

FeaturePlot(object = dataset1_TAA, features = c("percent.mt"))
FeaturePlot(object = dataset1_TAA, features = c("percent.ribo"))
FeaturePlot(object = dataset1_TAA, features = c("nCount_RNA"))
FeaturePlot(object = dataset1_TAA, features = c("nFeature_RNA"))
DimPlot(object = dataset1_TAA, label = T)
DotPlot(dataset1_TAA, features = markerGenes,cluster.idents = T) + RotatedAxis()
VlnPlot(dataset1_TAA, features = "PECAM1", pt.size = 0.25)

#, wt13

wt13.data <- Read10X_h5("/Users/lethithanhthuy/Desktop/20230410_Single\ fixed\ cell/Fixed\ cell\ TAA/filtered_feature_bc_matrix_WT13-TAA10wk.h5")
dim(wt13.data)
sampleName <-  'wt13'
dataset2_TAA <- CreateSeuratObject(counts = wt13.data, project = sampleName, min.cells = 3, min.features = 200)
dataset2_TAA
dataset2_TAA[["percent.mt"]] <- PercentageFeatureSet(dataset2_TAA, pattern = "^mt-")
dataset2_TAA[["percent.ribo"]] <- PercentageFeatureSet(dataset2_TAA, pattern = "^Rpl|^Rps|^Mrps|^Mrpl")
dataset2_TAA<- subset(dataset2_TAA, subset = percent.mt < 20 & nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 20000 )
VlnPlot(dataset2_TAA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4)
FeatureScatter(dataset2_TAA, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(dataset2_TAA, feature1 = "nFeature_RNA", feature2 = "percent.mt")
FeatureScatter(dataset2_TAA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dataset2_TAA <- NormalizeData(object = dataset2_TAA)
dataset2_TAA <- FindVariableFeatures(object = dataset2_TAA, nfeatures = 3000)
dataset2_TAA <- ScaleData(object = dataset2_TAA)
dataset2_TAA <- RunPCA(object = dataset2_TAA)
dataset2_TAA <- FindNeighbors(object = dataset2_TAA, dims = 1:30)
dataset2_TAA <- RunUMAP(object = dataset2_TAA, reduction = "pca", dims = 1:30)
dataset2_TAA <- FindClusters(object = dataset2_TAA)

FeaturePlot(object = dataset2_TAA, features = c("percent.mt"))
FeaturePlot(object = dataset2_TAA, features = c("percent.ribo"))
FeaturePlot(object = dataset2_TAA, features = c("nCount_RNA"))
FeaturePlot(object = dataset2_TAA, features = c("nFeature_RNA"))
DimPlot(object = dataset2_TAA, label = T)
DotPlot(dataset2_TAA, features = markerGenes,cluster.idents = T) + RotatedAxis()
VlnPlot(dataset2_TAA, features = "MRC1", pt.size = 0.25)

#, wt14

wt14.data <- Read10X_h5("/Users/lethithanhthuy/Desktop/20230410_Single\ fixed\ cell/Fixed\ cell\ TAA/filtered_feature_bc_matrix_WT14-TAA10wk.h5")
dim(wt14.data)
sampleName <-  'wt14'
dataset3_TAA <- CreateSeuratObject(counts = wt14.data, project = sampleName, min.cells = 3, min.features = 200)
dataset3_TAA
dataset3_TAA[["percent.mt"]] <- PercentageFeatureSet(dataset3_TAA, pattern = "^mt-")
dataset3_TAA[["percent.ribo"]] <- PercentageFeatureSet(dataset3_TAA, pattern = "^Rpl|^Rps|^Mrps|^Mrpl")
dataset3_TAA<- subset(dataset3_TAA, subset = percent.mt < 20 & nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 20000 )
VlnPlot(dataset3_TAA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4)
FeatureScatter(dataset3_TAA, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(dataset3_TAA, feature1 = "nFeature_RNA", feature2 = "percent.mt")
FeatureScatter(dataset3_TAA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dataset3_TAA <- NormalizeData(object = dataset3_TAA)
dataset3_TAA <- FindVariableFeatures(object = dataset3_TAA, nfeatures = 3000)
dataset3_TAA <- ScaleData(object = dataset3_TAA)
dataset3_TAA <- RunPCA(object = dataset3_TAA)
dataset3_TAA <- FindNeighbors(object = dataset3_TAA, dims = 1:30)
dataset3_TAA <- RunUMAP(object = dataset3_TAA, reduction = "pca", dims = 1:30)
dataset3_TAA <- FindClusters(object = dataset3_TAA)

FeaturePlot(object = dataset3_TAA, features = c("percent.mt"))
FeaturePlot(object = dataset3_TAA, features = c("percent.ribo"))
FeaturePlot(object = dataset3_TAA, features = c("nCount_RNA"))
FeaturePlot(object = dataset3_TAA, features = c("nFeature_RNA"))
DimPlot(object = dataset3_TAA, label = T)
DotPlot(dataset3_TAA, features = markerGenes,cluster.idents = T) + RotatedAxis()
VlnPlot(dataset3_TAA, features = "LYVE1", pt.size = 0.25)

#1, WT23

wt23.data <- Read10X_h5("/Users/lethithanhthuy/Desktop/20230410_Single\ fixed\ cell/Fixed\ cell\ TAA/filtered_feature_bc_matrix_WT23-TAA10wk-CYGB5wk.h5")
dim(wt23.data)
sampleName <-  'wt23'
dataset1_TAACY <- CreateSeuratObject(counts = wt23.data, project = sampleName, min.cells = 3, min.features = 200)
dataset1_TAACY
dataset1_TAACY[["percent.mt"]] <- PercentageFeatureSet(dataset1_TAACY, pattern = "^mt-")
dataset1_TAACY[["percent.ribo"]] <- PercentageFeatureSet(dataset1_TAACY, pattern = "^Rpl|^Rps|^Mrps|^Mrpl")
dataset1_TAACY<- subset(dataset1_TAACY, subset = percent.mt < 20 & nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 20000 )
VlnPlot(dataset1_TAACY, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4)
FeatureScatter(dataset1_TAACY, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(dataset1_TAACY, feature1 = "nFeature_RNA", feature2 = "percent.mt")
FeatureScatter(dataset1_TAACY, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dataset1_TAACY <- NormalizeData(object = dataset1_TAACY)
dataset1_TAACY <- FindVariableFeatures(object = dataset1_TAACY, nfeatures = 3000)
dataset1_TAACY <- ScaleData(object = dataset1_TAACY)
dataset1_TAACY <- RunPCA(object = dataset1_TAACY)
dataset1_TAACY <- FindNeighbors(object = dataset1_TAACY, dims = 1:30)
dataset1_TAACY <- RunUMAP(object = dataset1_TAACY, reduction = "pca", dims = 1:30)
dataset1_TAACY <- FindClusters(object = dataset1_TAACY)

FeaturePlot(object = dataset1_TAACY, features = c("percent.mt"))
FeaturePlot(object = dataset1_TAACY, features = c("percent.ribo"))
FeaturePlot(object = dataset1_TAACY, features = c("nCount_RNA"))
FeaturePlot(object = dataset1_TAACY, features = c("nFeature_RNA"))
DimPlot(object = dataset1_TAACY, label = T)
DotPlot(dataset1_TAACY, features = markerGenes,cluster.idents = T) + RotatedAxis()
VlnPlot(dataset1_TAACY, features = "PECAM1", pt.size = 0.25)

#, wt27

wt27.data <- Read10X_h5("/Users/lethithanhthuy/Desktop/20230410_Single\ fixed\ cell/Fixed\ cell\ TAA/filtered_feature_bc_matrix_WT27-TAA10wk-CYGB5wk.h5")
dim(wt27.data)
sampleName <-  'wt27'
dataset2_TAACY <- CreateSeuratObject(counts = wt27.data, project = sampleName, min.cells = 3, min.features = 200)
dataset2_TAACY
dataset2_TAACY[["percent.mt"]] <- PercentageFeatureSet(dataset2_TAACY, pattern = "^mt-")
dataset2_TAACY[["percent.ribo"]] <- PercentageFeatureSet(dataset2_TAACY, pattern = "^Rpl|^Rps|^Mrps|^Mrpl")
dataset2_TAACY<- subset(dataset2_TAACY, subset = percent.mt < 20 & nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 20000 )
VlnPlot(dataset2_TAACY, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4)
FeatureScatter(dataset2_TAACY, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(dataset2_TAACY, feature1 = "nFeature_RNA", feature2 = "percent.mt")
FeatureScatter(dataset2_TAACY, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dataset2_TAACY <- NormalizeData(object = dataset2_TAACY)
dataset2_TAACY <- FindVariableFeatures(object = dataset2_TAACY, nfeatures = 3000)
dataset2_TAACY <- ScaleData(object = dataset2_TAACY)
dataset2_TAACY <- RunPCA(object = dataset2_TAACY)
dataset2_TAACY <- FindNeighbors(object = dataset2_TAACY, dims = 1:30)
dataset2_TAACY <- RunUMAP(object = dataset2_TAACY, reduction = "pca", dims = 1:30)
dataset2_TAACY <- FindClusters(object = dataset2_TAACY)

FeaturePlot(object = dataset2_TAACY, features = c("percent.mt"))
FeaturePlot(object = dataset2_TAACY, features = c("percent.ribo"))
FeaturePlot(object = dataset2_TAACY, features = c("nCount_RNA"))
FeaturePlot(object = dataset2_TAACY, features = c("nFeature_RNA"))
DimPlot(object = dataset2_TAACY, label = T)
DotPlot(dataset2_TAACY, features = markerGenes,cluster.idents = T) + RotatedAxis()
VlnPlot(dataset2_TAACY, features = "MRC1", pt.size = 0.25)

#, wt28

wt28.data <- Read10X_h5("/Users/lethithanhthuy/Desktop/20230410_Single\ fixed\ cell/Fixed\ cell\ TAA/filtered_feature_bc_matrix_WT28-TAA10wk-CYGB5wk.h5")
dim(wt28.data)
sampleName <-  'wt28'
dataset3_TAACY <- CreateSeuratObject(counts = wt28.data, project = sampleName, min.cells = 3, min.features = 200)
dataset3_TAACY
dataset3_TAACY[["percent.mt"]] <- PercentageFeatureSet(dataset3_TAACY, pattern = "^mt-")
dataset3_TAACY[["percent.ribo"]] <- PercentageFeatureSet(dataset3_TAACY, pattern = "^Rpl|^Rps|^Mrps|^Mrpl")
dataset3_TAACY<- subset(dataset3_TAACY, subset = percent.mt < 20 & nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 20000 )
VlnPlot(dataset3_TAACY, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4)
FeatureScatter(dataset3_TAACY, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(dataset3_TAACY, feature1 = "nFeature_RNA", feature2 = "percent.mt")
FeatureScatter(dataset3_TAACY, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dataset3_TAACY <- NormalizeData(object = dataset3_TAACY)
dataset3_TAACY <- FindVariableFeatures(object = dataset3_TAACY, nfeatures = 3000)
dataset3_TAACY <- ScaleData(object = dataset3_TAACY)
dataset3_TAACY <- RunPCA(object = dataset3_TAACY)
dataset3_TAACY <- FindNeighbors(object = dataset3_TAACY, dims = 1:30)
dataset3_TAACY <- RunUMAP(object = dataset3_TAACY, reduction = "pca", dims = 1:30)
dataset3_TAACY <- FindClusters(object = dataset3_TAACY)

FeaturePlot(object = dataset3_TAACY, features = c("percent.mt"))
FeaturePlot(object = dataset3_TAACY, features = c("percent.ribo"))
FeaturePlot(object = dataset3_TAACY, features = c("nCount_RNA"))
FeaturePlot(object = dataset3_TAACY, features = c("nFeature_RNA"))
DimPlot(object = dataset3_TAACY, label = T)
DotPlot(dataset3_TAACY, features = markerGenes,cluster.idents = T) + RotatedAxis()
VlnPlot(data.combined, features = "Col3a1", pt.size = 0.25)

#1, WT1815

wt1815.data <- Read10X_h5("/Users/lethithanhthuy/Desktop/20230410_Single\ fixed\ cell/Fixed\ cell\ recovery/filtered_feature_bc_matrix_WT1815.h5")
dim(wt1815.data)
sampleName <-  'wt1815'
dataset1_TAARE <- CreateSeuratObject(counts = wt1815.data, project = sampleName, min.cells = 3, min.features = 200)
dataset1_TAARE
dataset1_TAARE[["percent.mt"]] <- PercentageFeatureSet(dataset1_TAARE, pattern = "^mt-")
dataset1_TAARE[["percent.ribo"]] <- PercentageFeatureSet(dataset1_TAARE, pattern = "^Rpl|^Rps|^Mrps|^Mrpl")
dataset1_TAARE<- subset(dataset1_TAARE, subset = percent.mt < 20 & nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 20000 )
VlnPlot(dataset1_TAARE, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4)
FeatureScatter(dataset1_TAARE, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(dataset1_TAARE, feature1 = "nFeature_RNA", feature2 = "percent.mt")
FeatureScatter(dataset1_TAARE, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dataset1_TAARE <- NormalizeData(object = dataset1_TAARE)
dataset1_TAARE <- FindVariableFeatures(object = dataset1_TAARE, nfeatures = 3000)
dataset1_TAARE <- ScaleData(object = dataset1_TAARE)
dataset1_TAARE <- RunPCA(object = dataset1_TAARE)
dataset1_TAARE <- FindNeighbors(object = dataset1_TAARE, dims = 1:30)
dataset1_TAARE <- RunUMAP(object = dataset1_TAARE, reduction = "pca", dims = 1:30)
dataset1_TAARE <- FindClusters(object = dataset1_TAARE)

FeaturePlot(object = dataset1_TAARE, features = c("percent.mt"))
FeaturePlot(object = dataset1_TAARE, features = c("percent.ribo"))
FeaturePlot(object = dataset1_TAARE, features = c("nCount_RNA"))
FeaturePlot(object = dataset1_TAARE, features = c("nFeature_RNA"))
DimPlot(object = dataset1_TAARE, label = T)
DotPlot(dataset1_TAARE, features = markerGenes,cluster.idents = T) + RotatedAxis()
VlnPlot(dataset1_TAARE, features = "PECAM1", pt.size = 0.25)

#, wt1816

wt1816.data <- Read10X_h5("/Users/lethithanhthuy/Desktop/20230410_Single\ fixed\ cell/Fixed\ cell\ recovery/filtered_feature_bc_matrix_WT1816.h5")
dim(wt1816.data)
sampleName <-  'wt1816'
dataset2_TAARE <- CreateSeuratObject(counts = wt1816.data, project = sampleName, min.cells = 3, min.features = 200)
dataset2_TAARE
dataset2_TAARE[["percent.mt"]] <- PercentageFeatureSet(dataset2_TAARE, pattern = "^mt-")
dataset2_TAARE[["percent.ribo"]] <- PercentageFeatureSet(dataset2_TAARE, pattern = "^Rpl|^Rps|^Mrps|^Mrpl")
dataset2_TAARE<- subset(dataset2_TAARE, subset = percent.mt < 20 & nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 20000 )
VlnPlot(dataset2_TAARE, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4)
FeatureScatter(dataset2_TAARE, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(dataset2_TAARE, feature1 = "nFeature_RNA", feature2 = "percent.mt")
FeatureScatter(dataset2_TAARE, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dataset2_TAARE <- NormalizeData(object = dataset2_TAARE)
dataset2_TAARE <- FindVariableFeatures(object = dataset2_TAARE, nfeatures = 3000)
dataset2_TAARE <- ScaleData(object = dataset2_TAARE)
dataset2_TAARE <- RunPCA(object = dataset2_TAARE)
dataset2_TAARE <- FindNeighbors(object = dataset2_TAARE, dims = 1:30)
dataset2_TAARE <- RunUMAP(object = dataset2_TAARE, reduction = "pca", dims = 1:30)
dataset2_TAARE <- FindClusters(object = dataset2_TAARE)

FeaturePlot(object = dataset2_TAARE, features = c("percent.mt"))
FeaturePlot(object = dataset2_TAARE, features = c("percent.ribo"))
FeaturePlot(object = dataset2_TAARE, features = c("nCount_RNA"))
FeaturePlot(object = dataset2_TAARE, features = c("nFeature_RNA"))
DimPlot(object = dataset2_TAARE, label = T)
DotPlot(dataset2_TAARE, features = markerGenes,cluster.idents = T) + RotatedAxis()
VlnPlot(dataset2_TAARE, features = "MRC1", pt.size = 0.25)

#, wt1817

wt1817.data <- Read10X_h5("/Users/lethithanhthuy/Desktop/20230410_Single\ fixed\ cell/Fixed\ cell\ recovery/filtered_feature_bc_matrix_WT1817.h5")
dim(wt1817.data)
sampleName <-  'wt1817'
dataset3_TAARE <- CreateSeuratObject(counts = wt1817.data, project = sampleName, min.cells = 3, min.features = 200)
dataset3_TAARE
dataset3_TAARE[["percent.mt"]] <- PercentageFeatureSet(dataset3_TAARE, pattern = "^mt-")
dataset3_TAARE[["percent.ribo"]] <- PercentageFeatureSet(dataset3_TAARE, pattern = "^Rpl|^Rps|^Mrps|^Mrpl")
dataset3_TAARE<- subset(dataset3_TAARE, subset = percent.mt < 20 & nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 20000 )
VlnPlot(dataset3_TAARE, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4)
FeatureScatter(dataset3_TAARE, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(dataset3_TAARE, feature1 = "nFeature_RNA", feature2 = "percent.mt")
FeatureScatter(dataset3_TAARE, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dataset3_TAARE <- NormalizeData(object = dataset3_TAARE)
dataset3_TAARE <- FindVariableFeatures(object = dataset3_TAARE, nfeatures = 3000)
dataset3_TAARE <- ScaleData(object = dataset3_TAARE)
dataset3_TAARE <- RunPCA(object = dataset3_TAARE)
dataset3_TAARE <- FindNeighbors(object = dataset3_TAARE, dims = 1:30)
dataset3_TAARE <- RunUMAP(object = dataset3_TAARE, reduction = "pca", dims = 1:30)
dataset3_TAARE <- FindClusters(object = dataset3_TAARE)

FeaturePlot(object = dataset3_TAARE, features = c("percent.mt"))
FeaturePlot(object = dataset3_TAARE, features = c("percent.ribo"))
FeaturePlot(object = dataset3_TAARE, features = c("nCount_RNA"))
FeaturePlot(object = dataset3_TAARE, features = c("nFeature_RNA"))
DimPlot(object = dataset3_TAARE, label = T)
DotPlot(dataset3_TAARE, features = markerGenes,cluster.idents = T) + RotatedAxis()
VlnPlot(dataset3_TAARE, features = "LYVE1", pt.size = 0.25)

# Intergration

data.list <- list(dataset1,dataset2,dataset3, dataset1_TAA, dataset2_TAA, dataset3_TAA, dataset1_TAARE, dataset2_TAARE, dataset3_TAARE)
names(data.list) <- c('wt1198','wt1199','wt1200', 'wt11', 'wt13', 'wt14', 'wt1815', 'wt1816', 'wt1817')
names(data.list) 

features <- SelectIntegrationFeatures(object.list = data.list,nfeatures = 3000)

data.anchors <- FindIntegrationAnchors(object.list = data.list, reduction = "rpca", anchor.features = features) #
date()

data.combined <- IntegrateData(anchorset = data.anchors)

DefaultAssay(data.combined) <- "integrated"

data.combined <- ScaleData(data.combined, verbose = FALSE)
data.combined <- RunPCA(data.combined, npcs = 50, verbose = FALSE)
ElbowPlot(data.combined, ndims = 50)
ndim <- 22
data.combined <- FindNeighbors(data.combined, reduction = "pca", dims = 1:ndim)
data.combined <- FindClusters(data.combined)#, resolution = 0.5
data.combined <- RunUMAP(data.combined, reduction = "pca", dims = 1:ndim)
DimPlot(data.combined, label = T)
DimPlot(data.combined,group.by = "orig.ident", label = F,shuffle = T)
DimPlot(data.combined,group.by = "group", label = F,shuffle = T)
DimPlot(data.combined,split.by = "orig.ident", label = F,shuffle = T)


DefaultAssay(data.combined) <- "RNA"
DotPlot(data.combined, features = markerGenes,cluster.idents = T) + RotatedAxis()
DotPlot(data.combined, features = markerGenes) + RotatedAxis()
dim(data.combined@assays$RNA@data)
dim(data.combined@assays$integrated@data)
DimPlot(data.combined, label = T)
data.combined@meta.data

# Rename group

data.list_cont <- list(dataset1,dataset2,dataset3)
names(data.list_cont) <- c('wt1198','wt1199','wt1200')
names(data.list_cont)

data.list_taa <- list(dataset1_TAA,dataset2_TAA,dataset3_TAA)
names(data.list_taa) <- c('wt11','wt13','wt14')
names(data.list_taa)

data.list_taare <- list(dataset1_TAACY,dataset2_TAACY,dataset3_TAACY)
names(data.list_taare) <- c('wt18','wt27','wt28')
names(data.list_taare)

library(Seurat)

data.combined$group <- ifelse(data.combined$orig.ident %in% names(data.list_cont), "CONT",
                              ifelse(data.combined$orig.ident %in% names(data.list_taa), "TAA",
                                     ifelse(data.combined$orig.ident %in% names(data.list_taare), "TAACY",
                                            NA
                                     )
                              )
)
data.combined@meta.data


data.list_cont <- list(dataset1,dataset2,dataset3)
names(data.list_cont) <- c('wt1198','wt1199','wt1200')
names(data.list_cont)

data.list_taa <- list(dataset1_TAA,dataset2_TAA,dataset3_TAA)
names(data.list_taa) <- c('wt11','wt13','wt14')
names(data.list_taa)

data.list_taare <- list(dataset1_TAARE,dataset2_TAARE,dataset3_TAARE)
names(data.list_taare) <- c('wt1815','wt1816','wt1817')
names(data.list_taare)

library(Seurat)

data.combined$group <- ifelse(data.combined$orig.ident %in% names(data.list_cont), "CONT",
                              ifelse(data.combined$orig.ident %in% names(data.list_taa), "TAA",
                                     ifelse(data.combined$orig.ident %in% names(data.list_taare), "TAARE",
                                            NA
                                     )
                              )
)
data.combined@meta.data

# Name groups

new.cluster.ids <- c(
  "LSECs",
  "LSECs",
  "Hepatocytes",
  "Macrophages",
  "HSCs",
  "Macrophages",
  "BECs",
  "Hepatocytes",
  "Hepatocytes",
  "B cells",
  "Hepatocytes",
  "HSCs",
  "NKT cells",
  "Hepatocytes",
  "LSECs",
  "Macrophages",
  "LSECs",
  "HSCs",
  "LSECs",
  "Macrophages",
  "Mesothelial cells",
  "Neutrophils",
  "LSECs",
  "Hepatocytes",
  "Dendritic cells",
  "Macrophages", 
  "NKT cells",
  "Macrophages",
  "Macrophages",
  "Hepatocytes",
  "LSECs",
  "Macrophages")
names(new.cluster.ids) <- levels(data.combined)
data.combined <- RenameIdents(data.combined, new.cluster.ids)
DimPlot(data.combined, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
plot1 <- DimPlot(data.combined, reduction = "umap", label = TRUE)
plot1
data.combined[["Cell.types"]] <- Idents(object = data.combined)
DimPlot(data.combined, group.by =  "seurat_clusters", label = T, label.size = 5)
DimPlot(data.combined, group.by =  "group")
DimPlot(data.combined, split.by = "group") + theme(aspect.ratio = 1) + theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black"),legend.text = element_text(size = 12))
DimPlot(data.combined, label.size = 4) + theme(aspect.ratio = 1)
table(data.combined@meta.data$orig.ident)
data.combined@meta.data
DimPlot_scCustom(seurat_object = data.combined, figure_plot = TRUE)


# Find specific marker


cluster0.markers <- FindMarkers(data.combined, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 10)
cluster1.markers <- FindMarkers(data.combined, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 10)
cluster2.markers <- FindMarkers(data.combined, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 10)
cluster3.markers <- FindMarkers(data.combined, ident.1 = 3, min.pct = 0.25)
head(cluster3.markers, n = 10)
cluster4.markers <- FindMarkers(data.combined, ident.1 = 4, min.pct = 0.25)
head(cluster4.markers, n = 10)
cluster5.markers <- FindMarkers(data.combined, ident.1 = 5, min.pct = 0.25)
head(cluster5.markers, n = 10)
cluster6.markers <- FindMarkers(data.combined, ident.1 = 6, min.pct = 0.25)
head(cluster6.markers, n = 10)
cluster7.markers <- FindMarkers(data.combined, ident.1 = 7, min.pct = 0.25)
head(cluster7.markers, n = 10)
cluster8.markers <- FindMarkers(data.combined, ident.1 = 8, min.pct = 0.25)
head(cluster8.markers, n = 10)
cluster9.markers <- FindMarkers(data.combined, ident.1 = 9, min.pct = 0.25)
head(cluster9.markers, n = 10)
cluster10.markers <- FindMarkers(data.combined, ident.1 = 10, min.pct = 0.25)
head(cluster10.markers, n = 10)
cluster11.markers <- FindMarkers(data.combined, ident.1 = 11, min.pct = 0.25)
head(cluster11.markers, n = 10)
cluster12.markers <- FindMarkers(data.combined, ident.1 = 12, min.pct = 0.25)
head(cluster12.markers, n = 10)
cluster13.markers <- FindMarkers(data.combined, ident.1 = 13, min.pct = 0.25)
head(cluster13.markers, n = 10)
cluster14.markers <- FindMarkers(data.combined, ident.1 = 14, min.pct = 0.25)
head(cluster14.markers, n = 10)
cluster15.markers <- FindMarkers(data.combined, ident.1 = 15, min.pct = 0.25)
head(cluster15.markers, n = 10)
cluster16.markers <- FindMarkers(data.combined, ident.1 = 16, min.pct = 0.25)
head(cluster16.markers, n = 10)
cluster17.markers <- FindMarkers(data.combined, ident.1 = 17, min.pct = 0.25)
head(cluster17.markers, n = 10)
cluster18.markers <- FindMarkers(data.combined, ident.1 = 18, min.pct = 0.25)
head(cluster18.markers, n = 10)
cluster19.markers <- FindMarkers(data.combined, ident.1 = 19, min.pct = 0.25)
head(cluster19.markers, n = 10)
cluster20.markers <- FindMarkers(data.combined, ident.1 = 20, min.pct = 0.25)
head(cluster20.markers, n = 10)
cluster21.markers <- FindMarkers(data.combined, ident.1 = 21, min.pct = 0.25)
head(cluster21.markers, n = 10)
cluster22.markers <- FindMarkers(data.combined, ident.1 = 22, min.pct = 0.25)
head(cluster22.markers, n = 10)
cluster23.markers <- FindMarkers(data.combined, ident.1 = 23, min.pct = 0.25)
head(cluster23.markers, n = 10)
cluster24.markers <- FindMarkers(data.combined, ident.1 = 24, min.pct = 0.25)
head(cluster24.markers, n = 10)
cluster25.markers <- FindMarkers(data.combined, ident.1 = 25, min.pct = 0.25)
head(cluster25.markers, n = 10)
cluster25.markers <- FindMarkers(data.combined, ident.1 = 25, min.pct = 0.25)
head(cluster25.markers, n = 10)
cluster26.markers <- FindMarkers(data.combined, ident.1 = 26, min.pct = 0.25)
head(cluster26.markers, n = 10)
cluster27.markers <- FindMarkers(data.combined, ident.1 = 27, min.pct = 0.25)
head(cluster27.markers, n = 10)
cluster28.markers <- FindMarkers(data.combined, ident.1 = 28, min.pct = 0.25)
head(cluster28.markers, n = 10)
cluster29.markers <- FindMarkers(data.combined, ident.1 = 29, min.pct = 0.25)
head(cluster29.markers, n = 10)
cluster30.markers <- FindMarkers(data.combined, ident.1 = 30, min.pct = 0.25)
head(cluster30.markers, n = 10)
cluster31.markers <- FindMarkers(data.combined, ident.1 = 31, min.pct = 0.25)
head(cluster31.markers, n = 10)
cluster32.markers <- FindMarkers(data.combined, ident.1 = 32, min.pct = 0.25)
head(cluster32.markers, n = 10)
cluster33.markers <- FindMarkers(data.combined, ident.1 = 33, min.pct = 0.25)
head(cluster33.markers, n = 10)
cluster34.markers <- FindMarkers(data.combined, ident.1 = 34, min.pct = 0.25)
head(cluster34.markers, n = 10)
cluster35.markers <- FindMarkers(data.combined, ident.1 = 35, min.pct = 0.25)
head(cluster35.markers, n = 10)
cluster36.markers <- FindMarkers(data.combined, ident.1 = 36, min.pct = 0.25)
head(cluster36.markers, n = 10)
cluster37.markers <- FindMarkers(data.combined, ident.1 = 37, min.pct = 0.25)
head(cluster37.markers, n = 10)
cluster38.markers <- FindMarkers(data.combined, ident.1 = 38, min.pct = 0.25)
head(cluster38.markers, n = 10)

library(ggplot2)

HSC_subset <- subset(data.combined, idents = "HSCs")

VlnPlot(HSC_subset, "Emilin1", split.by = "group")

group_colors <- c("CONT" = "blue", "TAA" = "red", "TAACY" = "green", "TAARE" = "purple")

VlnPlot(HSC_subset, "Colec11", group.by = "group")

VlnPlot + scale_fill_manual(values = group_colors)

# Subset cluster

data.combined <- SetIdent(data.combined, value = "Cell.types")
req_subset <- subset(data.combined, idents = c("6"), invert=F)
req_subset
DimPlot(req_subset,label = T, group.by = "ident",pt.size = 1)
DimPlot(req_subset,label = T, group.by = "Cell.types")
DimPlot(req_subset,label = T, group.by = "Cell.clusters")
# DimPlot(req_subset,label = T, group.by = "Cell.types")

DefaultAssay(req_subset) <- "integrated"
req_subset <- NormalizeData(req_subset, normalization.method = "LogNormalize", scale.factor = 10000)
req_subset <- FindVariableFeatures(req_subset, selection.method = "vst", nfeatures = 2000)
# all.genes <- rownames(req_subset)
req_subset <- ScaleData(req_subset, features = VariableFeatures(object = req_subset), vars.to.regress = c("nCount_RNA", "percent.mt"))
req_subset <- RunPCA(req_subset, features = VariableFeatures(object = req_subset))
ElbowPlot(req_subset, ndims = 50)
ndim=10
req_subset <- FindNeighbors(req_subset, dims = 1:ndim)
req_subset <- FindClusters(req_subset)
req_subset <- RunUMAP(req_subset, reduction = "pca", dims = 1:ndim,min.dist = 0.9, spread = 1.745)
req_subset <- subset(req_subset, idents = c("6","7","8","9", "10"), invert =T)
DimPlot(req_subset,label = F, pt.size = 1.5)  + theme(axis.text.x = element_text(size = 20, color = "black"), axis.text.y = element_text(size = 20, color = "black"), legend.text = element_text(size = 18), axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20)) +  guides(color = guide_legend(override.aes = list(size=14), ncol=1) ) + ggtitle(NULL)
DimPlot(req_subset,label = F, group.by = "group") + theme(axis.text.x = element_text(size = 16, color = "black"), axis.text.y = element_text(size = 16, color = "black"), legend.text = element_text(size = 16), axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20)) +  guides(color = guide_legend(override.aes = list(size=6), ncol=1))+ ggtitle(NULL)
DimPlot(req_subset,label = F, group.by = "group") + theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black"), legend.text = element_text(size = 12))
DimPlot(req_subset, group.by = "group") + scale_fill_manual(values=colors_clusters)
FeaturePlot(req_subset, "Tnfrsf11b", label = F)
VlnPlot(req_subset, "Mmp7", idents = c(0,12,18))
VlnPlot(req_subset, "Vegfa", group.by = "group")
DimPlot(req_subset)  +  guides(color = guide_legend(override.aes = list(size=10), ncol=1) )

req_subset <- readRDS("/Volumes/LaCie/20251005_Data\ analysis\ /20250105_Data\ analysis\ /20250105_BECs.rds")
DimPlot(req_subset)


# Cluster marker

cluster0.markers <- FindMarkers(req_subset, ident.1 = 0)
head(cluster0.markers, n = 30)
cluster1.markers <- FindMarkers(req_subset, ident.1 = 1)
head(cluster1.markers, n = 30)
cluster2.markers <- FindMarkers(req_subset, ident.1 = 2)
head(cluster2.markers, n = 30)
cluster3.markers <- FindMarkers(req_subset, ident.1 = 3)
head(cluster3.markers, n = 30)
cluster4.markers <- FindMarkers(req_subset, ident.1 = 4)
head(cluster4.markers, n = 30)
cluster5.markers <- FindMarkers(req_subset, ident.1 = 5)
head(cluster5.markers, n = 30)
cluster6.markers <- FindMarkers(req_subset, ident.1 = 6)
head(cluster6.markers, n = 30)
cluster7.markers <- FindMarkers(req_subset, ident.1 = 7)
head(cluster7.markers, n = 30)
cluster8.markers <- FindMarkers(req_subset, ident.1 = 8)
head(cluster8.markers, n = 30)
cluster9.markers <- FindMarkers(req_subset, ident.1 = 9)
head(cluster9.markers, n = 30)



DefaultAssay(req_subset) <- "RNA"
VlnPlot(req_subset, "Spp1", group.by = "group")
VlnPlot(req_subset, "Col4a1", group.by = "group")
ChVlnPlot(req_subset, "Cxcl5")
FeaturePlot(req_subset, "Cdh5", label = T)


write.csv(cluster0.markers, "/Users/lethithanhthuy/Desktop/Cluster0.csv")
write.csv(cluster1.markers, "/Users/lethithanhthuy/Desktop/Cluster1.csv")
write.csv(cluster2.markers, "/Users/lethithanhthuy/Desktop/Cluster2.csv")
write.csv(cluster3.markers, "/Users/lethithanhthuy/Desktop/Cluster3.csv")
write.csv(cluster4.markers, "/Users/lethithanhthuy/Desktop/Cluster4.csv")
write.csv(cluster5.markers, "/Users/lethithanhthuy/Desktop/Cluster5.csv")
write.csv(cluster6.markers, "/Users/lethithanhthuy/Desktop/Cluster6.csv")
write.csv(cluster7.markers, "/Users/lethithanhthuy/Desktop/Cluster7.csv")
write.csv(cluster8.markers, "/Users/lethithanhthuy/Desktop/Cluster8.csv")
write.csv(cluster9.markers, "/Users/lethithanhthuy/Desktop/Cluster9.csv")




# Run Progeny


library(progeny)
library(ggplot2)
library(tidyr)
library(readr)
library(pheatmap)
library(tibble)
library(viper)
CellsClusters <- data.frame(Cell = names(Idents(req_subset)), 
                            CellType = as.character(Idents(req_subset)),
                            stringsAsFactors = FALSE)

taa <- progeny(req_subset, scale=FALSE, organism="Mouse", top=500, perm=1, 
               return_assay = TRUE)
taa <- Seurat::ScaleData(taa, assay = "progeny") 
progeny_scores_df <- 
  as.data.frame(t(GetAssayData(taa, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)
summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 
paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)
progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength, 
                      max(summarized_progeny_scores_df), 
                      length.out=floor(paletteLength/2)))
progeny_hmap = pheatmap(t(summarized_progeny_scores_df[,-1]),fontsize=14, 
                        fontsize_row = 10, 
                        color=myColor, breaks = progenyBreaks, 
                        main = "PROGENy (500)", angle_col = 45,
                        treeheight_col = 0,  border_color = NA)

# Run Dorothea

library(dorothea)
library(decoupleR)
library(ggplot2)
library(tidyr)
library(readr)
library(pheatmap)
library(tibble)
library(viper)
library(ggrepel)
net <- get_dorothea(organism='mouse', levels=c('A', 'B', 'C'))
net
mat <- as.matrix(req_subset@assays$RNA@data)
acts <- run_wmean(mat=mat, net=net, .source='source', .target='target',
                  .mor='mor', times = 100, minsize = 5)
acts
req_subset[['tfswmean']] <- acts %>%
  filter(statistic == 'norm_wmean') %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)
DefaultAssay(object = req_subset) <- "tfswmean"
req_subset <- ScaleData(req_subset)
req_subset@assays$tfswmean@data <- req_subset@assays$tfswmean@scale.data
p1 <- DimPlot(data.combined, reduction = "tsne", label = TRUE, pt.size = 0.5) + 
  NoLegend() + ggtitle('Cell types')
p2 <- (FeaturePlot(data.combined, features = c("Smad3")) & 
         scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) +
  ggtitle('Smad3 activity')
DefaultAssay(object = data.combined) <- "RNA"
p3 <- FeaturePlot(data.combined, features = c("Smad3")) + ggtitle('Smad3 expression')
DefaultAssay(object = data.combined) <- "tfswmean"
p1 | p2 | p3
n_tfs <- 30
df <- t(as.matrix(req_subset@assays$tfswmean@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(req_subset)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))
tfs <- df %>%
  group_by(source) %>%
  summarise(std = sd(mean)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)
top_acts_mat <- df %>%
  filter(source %in% tfs) %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 2, length.out=floor(palette_length/2)))
pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks)

# Run propotion of cell
req_subset@meta.data
group <- req_subset$group
clusters <- req_subset$seurat_clusters
cell_type_table <- table(clusters, group)
prop_table <- prop.table(cell_type_table, margin = 1)
print(prop_table)
write.csv(prop_table, "/Users/lethithanhthuy/Desktop/Prop_table.csv")
library(ggplot2)

prop_table <- as.data.frame(prop_table)
prop_table$clusters <- rownames(prop_table)

ggplot(prop_table, aes(x = clusters, y = Freq, fill = group)) +
  geom_col(position = "fill", width = 0.5) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) +
  xlab("Cluster") +
  ylab("Proportion") +
  theme_bw(base_size = 15) +
  theme(legend.title = element_blank())

group <- req_subset$seurat_clusters
clusters <- req_subset$seurat_clusters
cell_type_table <- table(clusters, group)
prop_table <- prop.table(cell_type_table, margin = 1)
print(prop_)

prop_table <- as.data.frame(prop_table)
prop_table$group <- as.character(prop_table$group)

library(ggplot2)

ggplot(prop_table, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Normal") +
  ylab("Proportion") +
  scale_fill_manual(values = colorer) +
  theme(legend.title = element_blank())

pt <- table(Idents(req_subset), req_subset$group)
pt
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)
library(ggplot2)

ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Normal") +
  ylab("Proportion") +
  scale_fill_manual(c("#E41A1C", "#377EB8", "#4DAF4A")) +
  theme(legend.title = element_blank())

# Run EnrichR

install.packages("enrichR")
library(enrichR)


DimPlot(req_subset,label = T, pt.size = 1)
DefaultAssay(req_subset) <- "RNA"
markers <- FindAllMarkers(req_subset)
markers <- markers[order(markers$cluster,markers$avg_log2FC,decreasing = T),]
head(markers,20)

clustSel <- "9"
tmppos <- which(markers$cluster==clustSel)
tmp <- markers[tmppos,]
tmp[1:20,]
tmppos <- which(tmp$avg_log2FC > 0.5)
tmppos <- tmppos[1:if(length(tmppos)>100) 100 else length(tmppos)] #keep a max of top 100 genes
genelist <- tmp$gene[tmppos]
dbs <- listEnrichrDbs()
dbs <- c("GO_Cellular_Component_2023")
enriched <- enrichr(genelist, dbs)
head(enriched)

write.csv(enriched, "/Users/lethithanhthuy/Desktop/Enrich_Cluster0.csv")
write.csv(enriched, "/Users/lethithanhthuy/Desktop/Enrich_Cluster1.csv")
write.csv(enriched, "/Users/lethithanhthuy/Desktop/Enrich_Cluster2.csv")
write.csv(enriched, "/Users/lethithanhthuy/Desktop/Enrich_Cluster3.csv")
write.csv(enriched, "/Users/lethithanhthuy/Desktop/Enrich_Cluster4.csv")
write.csv(enriched, "/Users/lethithanhthuy/Desktop/Enrich_Cluster5.csv")
write.csv(enriched, "/Users/lethithanhthuy/Desktop/Enrich_Cluster6.csv")
write.csv(enriched, "/Users/lethithanhthuy/Desktop/Enrich_Cluster7.csv")
write.csv(enriched, "/Users/lethithanhthuy/Desktop/Enrich_Cluster8.csv")
write.csv(enriched, "/Users/lethithanhthuy/Desktop/Enrich_Cluster9.csv")



clustSel <- "9"
tmppos <- which(markers$cluster == clustSel)
tmp <- markers[tmppos, ]
tmp <- tmp[order(tmp$avg_log2FC), ]
tmppos <- head(which(tmp$avg_log2FC < 0), 100)  # Select top 100 genes with negative log2 fold change
genelist <- tmp$gene[tmppos]

# Perform enrichment analysis
dbs <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2023")
enriched <- enrichr(genelist, dbs)


write.csv(enriched, "/Users/lethithanhthuy/Desktop/Enrich_Cluster0.csv")
write.csv(enriched, "/Users/lethithanhthuy/Desktop/Enrich_Cluster1.csv")
write.csv(enriched, "/Users/lethithanhthuy/Desktop/Enrich_Cluster2.csv")
write.csv(enriched, "/Users/lethithanhthuy/Desktop/Enrich_Cluster3.csv")
write.csv(enriched, "/Users/lethithanhthuy/Desktop/Enrich_Cluster4.csv")
write.csv(enriched, "/Users/lethithanhthuy/Desktop/Enrich_Cluster5.csv")
write.csv(enriched, "/Users/lethithanhthuy/Desktop/Enrich_Cluster6.csv")
write.csv(enriched, "/Users/lethithanhthuy/Desktop/Enrich_Cluster7.csv")
write.csv(enriched, "/Users/lethithanhthuy/Desktop/Enrich_Cluster8.csv")
write.csv(enriched, "/Users/lethithanhthuy/Desktop/Enrich_Cluster9.csv")


# Display the results
head(enriched)


write.csv(enriched, "/Users/lethithanhthuy/Desktop/Enrich_Cluster0.csv")
write.csv(enriched, "/Users/lethithanhthuy/Desktop/Enrich_Cluster3.csv")

write.csv(enriched, "/Users/lethithanhthuy/Desktop/Enrich_Cluster1.csv")
write.csv(enriched, "/Users/lethithanhthuy/Desktop/Enrich_Cluster4.csv")
write.csv(enriched, "/Users/lethithanhthuy/Desktop/Enrich_Cluster6.csv")
write.csv(enriched, "/Users/lethithanhthuy/Desktop/Enrich_Cluster11.csv")

write.csv(enriched, "/Users/lethithanhthuy/Desktop/Enrich_Cluster5.csv")
write.csv(enriched, "/Users/lethithanhthuy/Desktop/Enrich_Cluster7.csv")
Ti




# Running gene set enrichment analysis

library(fgsea)

hallmark <- gmtPathways("/Users/lethithanhthuy/Desktop/2023-04-13\ All\ scRNA\ seq\ analysis/GSEA/mh.all.v2023.1.Mm.symbols.gmt.txt")
reactome <- gmtPathways("/Users/lethithanhthuy/Desktop/2023-04-13\ All\ scRNA\ seq\ analysis/GSEA/m2.cp.reactome.v2023.1.Mm.symbols.gmt.txt")
go <- fgsea::gmtPathways("/Users/lethithanhthuy/Desktop/2023-04-13\ All\ scRNA\ seq\ analysis/GSEA/m5.go.bp.v2023.1.Mm.symbols.gmt.txt")
gene_sets <- c(go)

gene_sets[1]



seurat <- AddModuleScore(req_subset, features=gene_sets["GOBP_POSITIVE_REGULATION_OF_ANGIOGENESIS"],
                         name="POSITIVE_REGULATION_OF_ANGIOGENESIS")

hist(seurat$POSITIVE_REGULATION_OF_ANGIOGENESIS1, breaks=50)

FeaturePlot(seurat, features="POSITIVE_REGULATION_OF_ANGIOGENESIS1", cols=c('lightgrey', 'blue'), order=T)

VlnPlot(seurat, features="POSITIVE_REGULATION_OF_VASCULAR_DEVELOPMENT1", group.by = "group")

# Rename subcluster

new.cluster.ids <- c(
  "BECs_S0",
  "BECs_S1",
  "BECs_S2",
  "BECs_S3",
  "BECs_S4",
  "BECs_S5")
names(new.cluster.ids) <- levels(req_subset)
req_subset <- RenameIdents(req_subset, new.cluster.ids)
DimPlot(req_subset, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
plot1 <- DimPlot(req_subset, reduction = "umap", label = TRUE)
plot1
req_subset[["Celltypes"]] <- Idents(object = req_subset)
DimPlot(req_subset, group.by = "group")
table(req_subset@active.ident)
req_subset@meta.data

VlnPlot(req_subset, "Hgf", group.by = "group")
VlnPlot(req_subset1, "Cdh5", idents = c("LSEC_Q1","LSEC_Q2", "LSEC_A1", "LSEC_A4", "LSEC_R1"), sort = F)

metadata_column <- "Cell.types"

VlnPlot(req_subset, features = "Gpx1")

req_subset$Cell.types <- factor(req_subset$Cell.types, levels = c("HEP_PC_CONT", "HEP_PC_TAA", "HEP_PC_RE"))

VlnPlot(req_subset, features = "Lgals3", idents = c("HEP_PC_CONT", "HEP_PC_TAA", "HEP_PC_RE"))+ scale_x_discrete(limit=c("HEP_PC_CONT","HEP_PC_TAA","HEP_PC_RE"))

VlnPlot(req_subset, features = "Trp53inp1", cols = c("#FFA500", "#7B68EE", "#00BFFF"), idents = c("HEP_PP_CONT", "HEP_PP_TAA", "HEP_PP_RE"))+ scale_x_discrete(limit=c("HEP_PP_CONT","HEP_PP_TAA","HEP_PP_RE"))


# FeaturePlot

genes <- c(unlist(strsplit("Epcam",split = ", ")))
# markers required
ptsize <- 1.5
genes <- intersect(genes, rownames(req_subset@assays$RNA@data))
FeaturePlot(object = req_subset, features = c(genes[1]), max.cutoff = "q99", pt.size = ptsize, order = T) + scale_color_viridis_c() + theme(axis.text.x = element_text(size = 20, color = "black"), axis.text.y = element_text(size = 20, color = "black"), legend.text = element_text(size = 18), axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20), plot.title = element_text(size = 22))  

# Plot cell composition

library(Seurat)
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)
library(RColorBrewer)
library(ggsci)
library(gridExtra)
library(reshape)

colors_samples = c(brewer.pal(5, "Set1"), brewer.pal(8, "Dark2"), pal_igv("default")(51))
colors_clusters = c(pal_d3("category10")(10), pal_d3("category20b")(20), pal_igv("default")(51))

colors_clusters <- c("CONT"="red2","TAA"="blue", "TAACY"="forestgreen")


Group = req_subset@meta.data[,"group"]

Idents(req_subset) <- req_subset@meta.data$Celltypes

barplot = data.frame(BEC_Clusters = Idents(req_subset), Cluster = Group)
barplot$Sample <- req_subset@meta.data$Celltypes
bb=ggplot(barplot, aes(BEC_Clusters, fill=Group)) + geom_bar() + theme_bw() + ylab("Cell count")+ xlab(NULL)+ theme(axis.text.x = element_text(size = 20, color = "black", angle = 45, hjust = 1), axis.text.y = element_text(size = 20, color = "black"), legend.text = element_text(size = 20), axis.title.y = element_text(size = 20))
bb2=ggplot(barplot, aes(BEC_Clusters, fill=Group)) + geom_bar(position="fill")+ ylab("Cell composition")+ xlab(NULL) + theme(axis.text.x = element_text(size = 18, color = "black", angle = 45, hjust = 1), axis.text.y = element_text(size = 16, color = "black"), legend.text = element_text(size = 18), axis.title.y = element_text(size = 20))

bb
bb2

df = table(barplot$Cluster, barplot$Sample, barplot$LSEC_Clusters)
df
dfm = melt(df)
dfm = subset(dfm, value > 0)
colnames(dfm) <- c("Cluster", "Sample", "group", "value")
dfm$Cluster <- factor(dfm$Cluster, levels = levels(req_subset@meta.data$Group))

dfm$totals <- 0
samps <- as.character(unique(dfm$Sample))
for(i in 1:length(samps)){
  curr = subset(req_subset@meta.data, orig.ident == samps[i])
  totaler = as.numeric(curr$cell_total[1])
  dfm$totals[dfm$Sample == samps[i]] <- totaler
}

# Monocle

library(monocle3)
library(SeuratWrappers)

cds <- SeuratWrappers::as.cell_data_set(req_subset)

devtools::install_github('cole-trapnell-lab/monocle3', lib="/path/to/your/personal/R-libraries")

devtools::install_github("satijalab/seurat-wrappers")

cds <- cluster_cells(cds, cluster_method = "louvain")

p1 <- plot_cells(cds, show_trajectory_graph = FALSE, label_groups_by_cluster = T)
p1
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
p2
wrap_plots(p1, p2)
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = TRUE, label_leaves = FALSE, label_branch_points = FALSE)
plot_cells(cds, label_groups_by_cluster = F, label_leaves = T, label_branch_points = T,graph_label_size=7) + theme(axis.text.x = element_text(size = 20, color = "black"), axis.text.y = element_text(size = 20, color = "black"), axis.title.x = element_text(size = 20, color = "black"), axis.title.y = element_text(size = 20, color = "black"))
cds <- order_cells(cds)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5) + theme(axis.text.x = element_text(size = 12, color = "black"), axis.text.y = element_text(size = 12, color = "black"), legend.text = element_text(size = 12))

rowData(cds)$gene_short_name <- row.names(rowData(cds))

cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)

pt <- pseudotime(cds)
colData(cds)$pseudotime <- as.numeric(pt)

pseudotime_bin <- c()
pt_bin_number_of_cells <- list()

pt_bin_number_of_cells <- list()
for (t in pt) {
  bin <- ""
  if (t >= 0 & t < 2.5) {
    bin <- "0.0-2.5"
  } else if (t >= 2.5 & t < 5.0) {
    bin <- "2.5-5.0"
  } else if (t >= 5.0 & t < 7.5) {
    bin <- "5.0-7.5"
  } else if (t >= 7.5 & t < 10.0) {
    bin <- "7.5-10.0"
  } else if (t >= 10.0 & t < 12.5) {
    bin <- "10.0-12.5"
  } else if (t >= 12.5 & t < 15.0) {
    bin <- "12.5-15.0"
  } else if (t >= 15.0 & t < 17.5) {
    bin <- "15.0-17.5"
  } else if (t >= 17.5 & t < 20.0) {
    bin <- "17.5-20.0"
  } else if (t >= 20.0 & t < 22.5) {
    bin <- "20.0-22.5"
  } else if (t >= 22.5 & t < 25.0) {
    bin <- "22.5-25.0"
  } else if (t >= 25.0 & t < 27.5) {
    bin <- "25.0-27.5"
  } else if (t >= 27.5) {
    bin <- ">= 27.5"
  }
  
  pseudotime_bin <- c(pseudotime_bin, bin)
  
  if (bin %in% names(pt_bin_number_of_cells) == "FALSE") {
    pt_bin_number_of_cells[bin] <- 0
  } else {
    pt_bin_number_of_cells[bin] <- as.numeric(pt_bin_number_of_cells[bin]) + 1
  }
}

pt_bin_number_of_cells

colData(cds)$pseudotime_bin <- pseudotime_bin

colData(cds)$pseudotime_bin <- factor(x = pseudotime_bin, levels = c("0.0-2.5", "2.5-5.0", "5.0-7.5", "7.5-10.0", "10.0-12.5", "12.5-15.0", "15.0-17.5", "17.5-20.0", "20.0-22.5","22.5-25.0", "25.0-27.5", ">= 27.5"))

req_subset$pseudotime_bin <- colData(cds)$pseudotime_bin

head(colData(cds))

tail(colData(cds))

quantile_out <- quantile(cds_pr_test_res$q_value, probs = c(0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1), na.rm=TRUE)
quantile_out

q_val_cutoff <- 0.05

pr_deg_ids <- row.names(subset(cds_pr_test_res, q_value < q_val_cutoff))

resolution <- 1e-03

gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=resolution)

class(cds[pr_deg_ids,])


table(gene_module_df$module)


genes <- c("Saa1", "Apoa1", "Efemp1")

genes_cds <- cds[rowData(cds)$gene_short_name %in% genes,]

plot_genes_in_pseudotime(genes_cds,
                         color_cells_by="Cell.types",
                         min_expr=0.5,
                         cell_size = 0.5, 
                         panel_order = genes,
                         ncol=2) +
  scale_color_manual(values=c("#377eb8", "#ff7f00",  "#4daf4a", "#e41a1c")) +
  geom_hline(yintercept=0) +
  theme(legend.title = element_blank(), 
        legend.position = "bottom",
        strip.text.x = element_text(size = 10, colour = "black", face = "bold"),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8))


library(ggplot2)
library(ggsci)

plot_genes_in_pseudotime(genes_cds,
                         color_cells_by = "group",
                         min_expr = 0.5,
                         cell_size = 0.5, 
                         panel_order = genes,
                         ncol = 2) +
  scale_color_manual(values = c("#377eb8", "#ff7f00",  "#4daf4a", "#e41a1c","#FFA500", "#FF69B4", "#FFA07A", "#FFFF00", "#BDB76B",  "#FF00FF", "#9932CC", "#7B68EE", "#7CFC00" )) +
  geom_hline(yintercept = 0) +
  theme(legend.title = element_blank(), 
        legend.position = "bottom",
        strip.text.x = element_text(size = 10, colour = "black", face = "bold"),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8))

plot_genes_in_pseudotime(genes_cds,
                         color_cells_by = "group",
                         min_expr = 0.5,
                         cell_size = 0.5, 
                         panel_order = genes,
                         ncol = 2) +
  scale_color_manual(values = c("#377eb8", "#ff7f00",  "#4daf4a")) +
  geom_hline(yintercept = 0) +
  theme

plot_genes_in_pseudotime(cds, 
                         genes = c("Cd36", "Apoa1", "Efemp1"), 
                         color_cells_by = "group", 
                         min_expr = 0.5,
                         cell_size = 0.5,
                         panel_order = c("Cd36", "Apoa1", "Efemp1"),
                         ncol = 2)


# Run heatmap top marker


markers <- FindAllMarkers(req_subset, 
                          assay = "RNA",
                          min.pct = 0.25,
                          thresh.use = 0.25,
                          only.pos = TRUE)

markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

req_subset1$seurat_clusters <- factor(req_subset1$seurat_clusters, levels = c("HEP_PC_CONT", "HEP_PC_TAA", "HEP_PC_RE"))

DoHeatmap(req_subset1, features = top10$gene)


DoHeatmap(req_subset1, features = top15$gene,
          group.by = "seurat_clusters")

library(RColorBrewer)

palette <- scales::hue_pal()(3)

DoHeatmap(req_subset1, features = top10$gene,
          group.by = "seurat_clusters",
          col = palette,
          raster = TRUE)

my_levels <- c("HEP_PC_CONT", "HEP_PC_TAA", "HEP_PC_RE")

req_subset1@meta.data

req_subset1@active.ident<- factor(req_subset1@active.ident, label=my_levels)

markers <- FindAllMarkers(req_subset1,
                          min.pct = 0.25,
                          thresh.use = 0.25,
                          only.pos = TRUE)



markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top20

DoHeatmap(req_subset1, features = top20$gene) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12)) + scale_fill_gradientn(colors = c("green", "white", "red"))


colnames(markers)


# Single cell gene set enrichment analysis

library(devtools)
install_github("guokai8/scGSVA")

library(scGSVA)  

hsko<-buildAnnot(species="mouse",keytype="SYMBOL",anntype="GO")


res<-scgsva(req_subset1,hsko)

Heatmap(res,group_by="Cell.types")

findPathway(res,group = "Cell.types")

sigPathway(res, group = "Cell.types")

vlnPlot(res,features="i.kappab.kinase.nf.kappab.signaling", group_by="Cell.types")

ridgePlot(res,features="cytokine.mediated.signaling.pathway",group_by="Cell.types")

featurePlot(res,features="leucocyte.cell.cell.adhesion", group_by="group")

library(scGSVA)
library(GO.db)

hsko <- buildAnnot(species = "mouse", keytype = "SYMBOL", anntype = "GO")

str(hsko)


# Run monocle 2

library(monocle)

data <- as(as.matrix(req_subset@assays$RNA@data), 'sparseMatrix')

pd <- new('AnnotatedDataFrame', data = req_subset@meta.data)

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

cds <- newCellDataSet(data,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

var_genes <- req_subset[["RNA"]]@var.features
ordering_genes <- var_genes

cds <- setOrderingFilter(cds, ordering_genes)
print(dim(exprs(cds)))

cds <- reduceDimension(cds,norm_method="none", 
                       reduction_method="DDRTree",
                       max_components=4,
                       scaling=TRUE,
                       verbose=TRUE,
                       pseudo_expr=0)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

req_subset <- CellCycleScoring(req_subset, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

head(req_subset[[]])

req_subset <- RunPCA(req_subset, features = c(s.genes, g2m.genes))
DimPlot(req_subset)

# Run GSVA



req_subset1 <- subset(req_subset, idents = c("BECs_S7", "BECs_S0", "BECs_S2", "BECs_S3", "BECs_S4", "BECs_S5"))

req_subset1$Cell.types <- factor(req_subset1$Cell.types, levels = c("BECs_Q7", "BECs_A0", "BECs_A2", "BECs_A3", "BECs_R4", "BECs_R5"))

library(scGSVA) 
hsko<-buildAnnot(species="mouse",keytype="SYMBOL",anntype="GO")
res<-scgsva(req_subset1,hsko)

Heatmap(res, group_by="Cell.types")

vlnPlot(res,features="Collecting.duct.acid.secretion", group_by="Cell.types")



pathways <- "cellular.response.to.lipopolysaccharide" 
res@obj$Cell.types <- factor(res@obj$, levels = c("BECs_S7", "BECs_S0", "BECs_S2", "BECs_S3", "BECs_S4", "BECs_S5"))
p <- vlnPlot(res, features = pathways, group_by = "active.ident")
df <- p$data
colnames(df)[1] <- "pathways"
head(df)
tmppos <- which(df$group=="BECs_S7")
tmpx <- df$pathways[tmppos]
tmpy <- df$pathways[-tmppos]
tmp <- wilcox.test(x=tmpx, y=tmpy)
tmp$p.value

ggplot(data = df, aes(x=group, y=pathways)) +
  geom_violin(aes(fill=group),scale = "width") +
  geom_boxplot(width=0.1)  +
  theme_classic()+
  ylab(pathways)+ theme(axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1), axis.text.y = element_text(size = 14, color = "black"), legend.text = element_blank())


# Run SCENIC

library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(tidyr)
library(AUCell)
library(RcisTarget)
library(GENIE3)
library(SCENIC)
library(tidyverse)
library(patchwork)
library(arrow)
library(ComplexHeatmap)

dir <- "HCC_Stem_scRNAseq"

figure <- "Figure_6_and_S6/6D_motif_regulon15"
outs_subpath <- paste0(dir,"/outs/",figure)
data_subpath <- paste0(dir,"/data/",figure)
dir.create(outs_subpath, recursive=TRUE)
dir.create(data_subpath, recursive=TRUE)

other_figure <- "Figure_5_and_S5/5BEFG_S5AB_integration_and_annotation"
other_outs_subpath <- paste0(dir,"/outs/",other_figure)

SCENIC_dir <- paste0(outs_subpath,"/SCENIC")
dir.create(SCENIC_dir)

param <- "pc10_res0.1"

exprMat <- as.matrix(req_subset@assays$RNA@data)

db_num <- 1
db_num <- 2

if (db_num==1) {
  dbDir <- paste0(data_subpath,"/cisTarget_databases") # RcisTarget databases location
} else if (db_num==2) {
  dbDir <- paste0(data_subpath,"/cisTarget_databases_2") # RcisTarget databases location
}

org <- "mgi"
myDatasetTitle <- "HCC.final.SCT_SCENIC"
data(defaultDbNames)

if (db_num==1) {
  dbs <- defaultDbNames[[org]]
  minCountsPerGene_factor <- 3
} else if (db_num==2) {
  dbs <- c("mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather", "mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather")
  names(dbs) <- c("500bp", "10kb")
  minCountsPerGene_factor <- 1
}

scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10)


org<-"hngc"
dbDir<- "/Users/minhduc/HCC_Stem_scRNAseq/data/Figure_6_and_S6/6D_motif_regulon/cisTarget_databases/"
myDatasetTitle<-"SCENIC regulon analysis on mouse"
defaultDbNames$mgi[1]<-"mm9-500bp-upstream-7species.mc9nr.genes_vs_motifs.rankings.feather"
defaultDbNames$mgi[2]<-"mm9-tss-centered-10kb-7species.mc9nr.genes_vs_motifs.rankings.feather"
dbs <- defaultDbNames[[org]]
minCountsPerGene_factor <- 1
scenicOptions <- initializeScenic(org="hgnc", dbs=defaultDbNames[["hgnc"]], dbDir= dbDir, nCores=10)

org<-"mgi"
dbDir<- "/Users/lethithanhthuy/Desktop/data_subpath"
myDatasetTitle<-"SCENIC regulon analysis on mouse"
defaultDbNames$mgi[1]<-"mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
defaultDbNames$mgi[2]<-"mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
dbs <- defaultDbNames[[org]]
minCountsPerGene_factor <- 1
scenicOptions <- initializeScenic(org = org, dbDir = dbDir, dbs = dbs, datasetTitle = myDatasetTitle, nCores = 10)
data(list="motifAnnotations_mgi", package="RcisTarget")
motifAnnotations_mgi <- motifAnnotations_mgi

library(SCENIC)
org <- "mgi"
dbDir <- "/Users/lethithanhthuy/Desktop/data_subpath"
myDatasetTitle <- "SCENIC example on Mouse liver" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10) 

data(list="motifAnnotations_mgi_v9", package="RcisTarget")

motifAnnotations_mgi <- motifAnnotations_mgi_v9

genesKept <- geneFiltering(exprMat, 
                           scenicOptions, 
                           minCountsPerGene=minCountsPerGene_factor*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)

interestingGenes <- c("Smad3")
interestingGenes[which(!interestingGenes %in% genesKept)]

exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)

runCorrelation(exprMat_filtered, scenicOptions)


runGenie3(exprMat_filtered, scenicOptions)

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)

saveRDS(scenicOptions, file = "/Users/lethithanhthuy/Desktop/scenicOption.Rds")

library(BiocParallel)

library(parallel)

BPPARAM <- BiocParallel::bpparam()

BPPARAM$workers <- 1

scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)

library("doSNOW")

library("doParallel")

library("doMPI")
library("doMC")
library("R2HTML")

install.packages("doMPI")

install.packages("Rmpi")

scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat)

scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)

tsneAUC(scenicOptions, aucType = "AUC")

# Motif Entichment preview

motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")

tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Etv1"]

viewMotifs(tableSubset)

# regulon target info

regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")

tableSubset <- regulonTargetsInfo[TF="Nf1" & highContAnnot==TRUE]

# Cell.type specificregulator

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")

cellInfo <- data.frame(req_subset1@meta.data)

rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation = cellInfo[colnames(regulonAUC)], "Cell" )

rssPlot <- plotRSS(rss)

plotly:: ggplotly(rssPlot$plot)

# Heatmap creation

regulonAUCmatrix <- as.matrix(regulonAUC[,2:ncol(regulonAUC)]) %>% 
  
  rownames(regulonAUCmatrix ) <- colnames(regulonAUC)[2:ncol(regulonAUC)]

colnames(regulonAUCmatrix) <- rownames(req_subset1@meta.data)

regulonAUCmatrixSeurat <- CreateSeuratObject(counts = regulonAUCmatrix )

req_subset1@assays$RegulonAUC <- regulonAUCmatrixSeurat@assays$RNA

MSC_Regulon_ClusterID <- lapply(split(data.combined@meta.data,list(data.combined$Cell.types)), function(x)rownames(x))

MSC_RegulonClusterMean <- t(apply(regulonAUCmatrix,1,function(x){
  lapply(MSC_Regulon_ClusterID,function(ID){
    mean(x[ID]) 
  }) %>% unlist
}))

MSC_RegulonClusterMean <- MSC_RegulonClusterMean[apply(MSC_RegulonClusterMean,1,sum)>0,]

MSC_RegulonClusterMeanM <- MSC_RegulonClusterMean
MSC_RegulonClusterMeanM[,1:ncol(MSC_RegulonClusterMeanM)] <- t(apply(MSC_RegulonClusterMeanM,1,scale))
Sample_order <- seurat_clusters

MSC_RegulonClusterMeanM[MSC_RegulonClusterMeanM< -3] <- -3
MSC_RegulonClusterMeanM[MSC_RegulonClusterMeanM> 3] <- 3


regulons <- loadInt(scenicOptions, "aucell_regulons")
head(cbind(onlyNonDuplicatedExtended(names(regulons))))

regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="Stat6" & highConfAnnot==TRUE]
viewMotifs(tableSubset, options=list(pageLength=5))

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$Cell.types),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
write.csv(regulonActivity_byCellType_Scaled, "/Users/lethithanhthuy/Desktop/Cholangiocytes_regulonActivity_byCellType_Scaled.csv")
library(plotly)

hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity"), row_names_gp=grid::gpar(fontsize=6)))
regulonOrder <- rownames(regulonActivity_byCellType_Scaled)[row_order(hm)]
row_or
topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "Cell.types", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>1.5),]
viewTable(topRegulators)
library(grid)
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name = "Regulon activity"), row_names_gp = grid::gpar(fontsize = 6))


ComplexHeatmap::Heatmap(topRegulators, name="Regulon activity")
library(ggplot2)

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "Cell.types"])
rssPlot <- plotRSS(rss)*
  
  plotly::ggplotly(rssPlot$plot,) x.axis = 3

library(plotly)

Heatmap(selected_data, row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize= 10))

write.csv(regulonActivity_byCellType, "/Users/lethithanhthuy/Desktop/RegulonActivity_byCellType_TAA.csv")

saveRDS(scenicOptions, file = "/Users/lethithanhthuy/Desktop/scenicOption_TAA.Rds")

integrated.sub <- as.Seurat(cds, assay = NULL,project = "cell_data_set")
FeaturePlot(integrated.sub, "monocle3_pseudotime")

cellInfo <- data.frame(seuratCluster=Idents(req_subset))


regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$Cell.types),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

ComplexHeatmap::Heatmap(selected_data, name="Regulon activity", row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize= 10))

# Assuming you have regulonActivity_byCellType_Scaled data
# Check for NA/NaN/Inf values
any_na <- any(is.na(regulonActivity_byCellType_Scaled) | is.nan(regulonActivity_byCellType_Scaled) | is.infinite(regulonActivity_byCellType_Scaled))

if (any_na) {
  # Handle NA/NaN/Inf values here (e.g., impute or remove)
  # regulonActivity_byCellType_Scaled <- some_handling_function(regulonActivity_byCellType_Scaled)
}

# Z-score normalization
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType_Scaled), center = TRUE, scale = TRUE))

# Create the heatmap
hm <- ComplexHeatmap::Heatmap(selected_data, name = "Regulon activity")
hm


get_top_TFs <- function(cell_type_activity, N) {
  sorted_TFs <- sort(cell_type_activity, decreasing = TRUE)
  top_TFs <- names(sorted_TFs)[1:N]
  return(top_TFs)
}

top_N_TFs <- 5
top_TFs_list <- list()

cell_type <- c("LSEC_Q0", "LSEC_Q3", "LSEC_A2", "LSEC_A4", "LSEC_R1")

for (Cell.types in colnames(selected_data)) {
  top_TFs_list[[Cell.types]] <- get_top_TFs(selected_data[, Cell.types], top_N_TFs)
}

top_TFs_list

top_TFs <- unlist(top_TFs_list)

selected_data_top_TFs <- selected_data[, top_TFs]

regulonActivity_top_TFs <- selected_data[, list(top_TFs_list)]

valid_TFs <- intersect(unlist(top_TFs_list), colnames(regulonActivity_byCellType_Scaled))
top_TFs_list <- lapply(top_TFs_list, function(tf_list) intersect(tf_list, valid_TFs))

regulonActivity_top_TFs <- selected_data[, unlist(top_TFs_list)]
selected_data

heatmap(top_TFs, Rowv = NA, Colv = NA,
        col = colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(100),
        scale = "none", main = "Top 5 Transcription Factors by Cell Type",
        xlab = "Cell Type", ylab = "Transcription Factor")

head(regulonActivity_top_TFs)

valid_TFs <- intersect(unlist(top_TFs_list), colnames(regulonActivity_byCellType_Scaled))
top_TFs_list <- lapply(top_TFs_list, function(tf_list) intersect(tf_list, valid_TFs))
regulonActivity_top_TFs <- regulonActivity_byCellType_Scaled[, unlist(top_TFs_list)]



heatmap(regulonActivity_top_TFs, Rowv = NA, Colv = NA,
        col = colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(100),
        scale = "none", main = "Top 5 Transcription Factors by Cell Type",
        xlab = "Cell Type", ylab = "Transcription Factor")

is.na(regulonActivity_byCellType_Scaled) %>% table()

dim(regulonActivity_byCellType_Scaled)

regulonActivity_byCellType_Scaled <- t(regulonActivity_byCellType_Scaled[complete.cases(t(regulonActivity_byCellType_Scaled)), ])

regulonActivity_byCellType_Scaled <- t(regulonActivity_byCellType_Scaled)

heatmap(regulonActivity_byCellType_Scaled, features=rownames(regulonActivity_byCellType_Scaled)[1:3],
        center=TRUE, symmetric=TRUE)

library(dplyr)

regulonActivity_byCellType_Scaled %>% dplyr::select("LSEC_Q0", "LSEC_Q3", "LSEC_A2", "LSEC_A4", "LSEC_R1") 

dplyr::select(iris, Sepal.Width, Petal.Length, Species)

data_without_na <- regulonActivity_byCellType_Scaled[, complete.cases(regulonActivity_byCellType_Scaled)]

columns_to_keep <- c("LSEC_Q0", "LSEC_Q3", "LSEC_A2", "LSEC_A4", "LSEC_R1")

selected_data <- regulonActivity_byCellType_Scaled[, columns_to_keep]

for (Cell.types in colnames(selected_data)) {
  top_TFs_list[[Cell.types]] <- get_top_TFs(selected_data[, Cell.types], top_N = 5)
}

topRegulators <- reshape2::melt(selected_data)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)

ComplexHeatmap::Heatmap(topRegulators, name="Regulon activity")

ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity")

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "Cell.types"])
rssPlot <- plotRSS(rss)
library(plotly)
plotly_plot <- plotly::ggplotly(rssPlot$plot)

plotly_plot <- plotly_plot %>%
  plotly::layout(yaxis = list(tickfont = list(size = 6)))

plotly_plot

plotRSS_oneSet(rss, setName = "BECs_R5")

library(Seurat)
dr_coords <- Embeddings(req_subset1, reduction="tsne")

tfs <- c("Sox10","Irf1","Sox9", "Dlx5")
par(mfrow=c(2,2))
AUCell::AUCell_plotTSNE(dr_coords, cellsAUC=selectRegulons(regulonAUC, tfs), plots = "AUC")

motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Sox17"]
viewMotifs(tableSubset) 

Regulon <- readRDS("/Users/lethithanhthuy/Desktop/3.4_regulonAUC.Rds")
Regulon <- Regulon@assays@data@listData$AUC
dim(Regulon)

Regulon_low <- Regulon[grepl("extended", rownames(Regulon)),]
row.names(Regulon_low) <- str_split(row.names(Regulon_low), "_", simplify = T)[,1]

Regulon_low <- Regulon_low[rowSums(Regulon_low>0)>0.5*length(colnames(Regulon_low)),]
dim(Regulon_low)

req_subset1[["Regulon"]] <- CreateAssayObject(counts = Regulon_low) 

DefaultAssay(req_subset1) <- "Regulon"
req_subset$Cell.types <- factor(req_subset$Cell.types, levels = c("LSEC_S0", "LSEC_Q3", "LSEC_A2", "LSEC_A4", "LSEC_R1"))
req_subset1$group <- factor(req_subset1$group, levels = c("CONT", "TAA", "TAARE"))
dim(req_subset1)
TF <- c("Sox18")

for (TF in c("Hnf4a", "Cebpa", "Foxa3")){
  VlnPlot(req_subset1, features = TF, 
          group.by = "Cell.types", 
          pt.size=0,
          cols = c("#e41a1c", "#377eb8", "#4daf4a", "#ff7f00", "blue"))+
    geom_boxplot(width=0.1, fill="white", outlier.shape = NA, show.legend = FALSE)+
    labs(x="Cluster", y = "AUCell", title = paste0(TF, " motif regulon")) +
    scale_x_discrete(labels=c("3" = "LSEC_Q0", "0" = "LSEC_Q3",
                              "2" = "LSEC_A2", "1" = "LSEC_A4", "4" = "LSEC_R1") +
                       theme(axis.text.x = element_text(face="plain", color="black", 
                                                        size=10, angle=0, hjust = 0.5),
                             axis.text.y = element_text(face="plain", color="black", 
                                                        size=10, angle=0),
                             axis.title.x = element_text(face="bold", color="black", 
                                                         size=12, angle=0),
                             axis.title.y = element_text(face="bold", color="black", 
                                                         size=12, angle=90),
                             plot.title = element_text(color="black", size=14, face = "bold", hjust = 1), 
                             legend.position = "none")+
                       geom_signif(comparisons = list(c("LSEC_Q0", "LSEC_Q3"),
                                                      c("LSEC_A2", "LSEC_A4", "LSEC_R1")))
}

for (TF in c("Hnf4a")) {
  VlnPlot(req_subset1, features = TF, 
          group.by = "Cell.types", 
          pt.size = 0,
          cols = c("#e41a1c", "#377eb8", "#4daf4a", "#ff7f00", "blue")) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, show.legend = FALSE) +
    labs(x = "Cluster", y = "AUCell", title = paste0(TF, " motif regulon")) +
    scale_x_discrete(labels = c("3" = "LSEC_Q0", "0" = "LSEC_Q3",
                                "2" = "LSEC_A2", "1" = "LSEC_A4", "4" = "LSEC_R1")) +
    theme(axis.text.x = element_text(face = "plain", color = "black", 
                                     size = 10, angle = 0, hjust = 0.5),
          axis.text.y = element_text(face = "plain", color = "black", 
                                     size = 10, angle = 0),
          axis.title.x = element_text(face = "bold", color = "black", 
                                      size = 12, angle = 0),
          axis.title.y = element_text(face = "bold", color = "black", 
                                      size = 12, angle = 90),
          plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 1), 
          legend.position = "none") +
    geom_signif(comparisons = list(c("LSEC_Q0", "LSEC_Q3"),
                                   c("LSEC_A2", "LSEC_A4", "LSEC_R1")))
}

library(ggsignif)     

length(regulons)

sum(lengths(regulons)>=10)

viewTable(cbind(nGenes=lengths(regulons)), options=list(pageLength=10))

tableSubset <- motifEnrichment[TF=="Nf1"]
viewMotifs(tableSubset, colsToShow = c("logo", "NES", "TF" ,"Annotation"), options=list(pageLength=5))

binarizeAUC <- function(auc, thresholds)
{
  thresholds <- thresholds[intersect(names(thresholds), rownames(auc))]
  regulonsCells <- setNames(lapply(names(thresholds), 
                                   function(x) {
                                     trh <- thresholds[x]
                                     names(which(getAUC(auc)[x,]>trh))
                                   }),names(thresholds))
  
  regulonActivity <- reshape2::melt(regulonsCells)
  binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2]))
  class(binaryRegulonActivity) <- "matrix"  
  
  return(binaryRegulonActivity)
  
}

binaryRegulonActivity <- binarizeAUC(regulonAUC, regulonAucThresholds)
dim(binaryRegulonActivity)

nCells <- 1000
set.seed(123)
cellsSelected <- sample(colnames(regulonAUC), nCells) 
binAct_subset <- binaryRegulonActivity[, which(colnames(binaryRegulonActivity) %in% cellsSelected)]
dim(binAct_subset)

# Run Heatmap

"BECs_Q7","BECs_A0", "BECs_A2", "BECs_A3", "BECs_R4", "BECs_R5"

# Making heatmap

write.csv(regulonActivity_byCellType_Scaled, "/Users/lethithanhthuy/Desktop/Cholangiocytes_regulonActivity_byCellType_Scaled.csv")

heatmap <- read_csv("/Users/lethithanhthuy/Desktop/Cholangiocytes_regulonActivity_byCellType_Scaled.csv")

library(dplyr)

heatmap <- heatmap %>%
  rename(
    BECs_Q7 = BECs_S7,
    BECs_A0 = BECs_S0,
    BECs_A2 = BECs_S2,
    BECs_A3 = BECs_S3,
    BECs_R4 = BECs_S4,
    BECs_R5 = BECs_S5,
  )

selected_cell_types <- c("BECs_S7","BECs_S0", "BECs_S2", "BECs_S3", "BECs_S4", "BECs_S5")
selected_data <- heatmap[, c("...1", selected_cell_types)]

top_5_tfs <- lapply(selected_cell_types, function(cell_type) {
  sorted_data <- selected_data[order(selected_data[[cell_type]], decreasing = TRUE), ]
  top_5 <- head(sorted_data, 10)
  top_5
})

top_5_tfs_combined <- do.call(rbind, top_5_tfs)

top_5_tfs_combined$Cell_Type <- factor(rep(selected_cell_types, each = 10))

top_5_tfs_matrix <- as.matrix(top_5_tfs_combined[, selected_cell_types])
rownames(top_5_tfs_matrix) <- top_5_tfs_combined$'...1'

library(ComplexHeatmap)
library(circlize)

heatmap_object <- Heatmap(
  top_5_tfs_matrix,
  col = colorRamp2(c(min(top_5_tfs_matrix), 0, max(top_5_tfs_matrix)), c("green", "white", "red")),
  name = "Z-score",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  column_title = "",
  row_title = "Cholangiocytes regulons",
  row_names_gp = gpar(fontsize = 8),
  top_annotation = HeatmapAnnotation(df = data.frame(Cell_Type = selected_cell_types),
                                     show_legend = TRUE)
)

draw(heatmap_object)

# Violin PLot statistic

gene <-"Nfia"
p <- VlnPlot(req_subset1, features = gene, idents = c("BECs_Q7", "BECs_A0", "BECs_A2", "BECs_A3", "BECs_R4", "BECs_R5")) 
df <- p$data
colnames(df)[1] <- "gene"
head(df)
tmppos <- which(df$ident=="BECs_Q7")
tmpx <- df$gene[tmppos]
tmpy <- df$gene[-tmppos]
tmp <- wilcox.test(x=tmpx, y=tmpy)
tmp$p.value
ggplot(data = df, aes(x=ident, y=gene)) +
  geom_violin(aes(fill=ident),scale = "width") + 
  geom_boxplot(width=0.1)  +
  theme_classic()+
  ylab(gene)+
  ggtitle(paste0("Mm_",gene,"_Mann_whitney_p_",format(tmp$p.value,digits = 3))) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# AnalysisAalysis

plotter <- req_subset@meta.data
plotter$Mmp7 <- req_subset@assays$integrated@data["Mmp7",]
plotter$Cxcl16 <- req_subset@assays$integrated@data["Cxcl16",]


plotter$group <- plotter$group
plotter$pouch.status2[plotter$APOE < 0.25] <- "other"
plotter$pouch.status2[plotter$C1QC < 0.25] <- "other"

ggAPOE = ggplot(plotter, aes(Mmp7, Cxcl16, color=group)) + 
  geom_point(size=0.5) +
  theme_bw() + facet_wrap(~group, nrow=1) +
  geom_vline(xintercept =0.25) + geom_hline(yintercept = 0.25)
ggAPOE

VlnPlot(req_subset, "Mrs1")

# Stacked Violin Plot

DefaultAssay(data.combined)<- "RNA"
features <- c("Krt7","Krt19", "Mmp7", "Lcn2", "Scube3","Thbs1","Tgfb2", "Pdgfb", "Cxcl16", "Tnfrsf12a", "Icam1", "Vcam1")

VlnPlot(req_subset, features, stack = TRUE, sort = F, flip = F) +
  theme(legend.position = "none") + theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black"), legend.text = element_text(size = 16)) + ylab(NULL)

# Single cell

FeaturePlot(data.combined, "Vwf")

library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(assertthat)
library(tidyquant)
library(clusterProfiler)
library(viridis)
BiocManager::install("org.Hs.eg.db")
library(org.Mm.eg.db)

req_subset$Cell.types <- factor(req_subset$Cell.types, levels = c("BECs_S0", "BECs_S1", "BECs_S2", "BECs_S3", "BECs_S4", "BECs_S5"))

req_subset <- FindAllMarkers(req_subset)
genes <- bitr(unique(req_subset$gene), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = 'org.Mm.eg.db')
colnames(genes)[1] <- "gene"

req_subset <- req_subset %>% 
  left_join(genes[!duplicated(genes$gene),]) %>% 
  na.omit() %>% 
  group_by(cluster) %>% 
  filter(avg_log2FC>0)
diff_go_mf <- compareCluster(ENTREZID ~ cluster,
                             data=req_subset, 
                             fun = enrichGO,
                             OrgDb = 'org.Mm.eg.db',
                             ont = "BP")

clusterProfiler::dotplot(diff_go_mf)

diff_go_mf <- compareCluster(ENTREZID ~ cluster,
                             data=req_subset, 
                             fun = enrichGO,
                             OrgDb = 'org.Mm.eg.db',
                             ont = "BP")
top_pathways <- 3

clusterProfiler::dotplot(diff_go_mf, showCategory = top_pathways)

clusterProfiler::dotplot(diff_go_mf, showCategory = top_pathways) + theme(axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1), axis.text.y = element_text(size = 12, color = "black"), legend.text = element_text(size = 14))+ xlab(NULL)


# Enhanced Vocalno plot

c4 <- FindMarkers(req_subset, 
                  ident.1 = "TAA",
                  ident.2 = "TAARE", 
                  logfc.threshold = 0.01,
                  min.pct = 0.01,
                  group.by = "group")

library(EnhancedVolcano)
EnhancedVolcano(c4, lab = rownames(c4), 
                x="avg_log2FC", y="p_val", 
                FCcutoff = .1,
                pCutoff = .0001,
                #labSize = 6, 
                subtitle = element_blank(), title = "TAARE (left) vs. TAA (right)",
                caption = element_blank(),
                parseLabels = T,
                drawConnectors = T,
                max.overlaps = 20,
                colConnectors = "black",
                arrowheads = F,
                pointSize = 3,
                xlim = c(-1,1),
                ylim = c(0, 40)
) +
  theme(panel.grid.minor = element_blank()) 