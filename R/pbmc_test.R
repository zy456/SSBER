library(reticulate)
library(scran)
library(Seurat)
library(BiocNeighbors)
library(edgeR)
library(batchelor)
library(cowplot)
library(readxl)
library(mclust)
library(magrittr)
library(dplyr)
#get meta information
metadata <- readRDS("pbmc_metadata.rds")
# get data
pbmc<- readRDS("pbmc_expression_matrix.rds")
pbmc <- CreateSeuratObject(counts = pbmc, meta.data = metadata)

pbmc <- NormalizeData(pbmc, verbose = FALSE)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 6000, verbose = FALSE)
# clustering
pbmc <- ScaleData(pbmc, verbose = FALSE)
pbmc <- RunPCA(pbmc, npcs = 20, verbose = FALSE)
pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:20)
p1 <- DimPlot(pbmc, reduction = "umap", group.by = "tech",pt.size = 1)
p2 <- DimPlot(pbmc, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE,pt.size = 1,label.size = 8) 
plot_grid(p1, p2)
genename<-VariableFeatures(pbmc)

#harmony
matrix<-read.csv("./pbmc/integrated_data.csv")
#matrix<-t(matrix)
#meta_new<-rbind(data1,data2,data3)
colnames(matrix)<-row.names(metadata)
row.names(matrix)<-c("pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","pc11","pc12",
                     "pc13","pc14","pc15","pc16","pc17","pc18","pc19","pc20")


integrated_data <- CreateSeuratObject(counts = matrix,meta.data = metadata)

pca<-pbmc@reductions$pca
pca@cell.embeddings<-t(matrix)
pca@feature.loadings<-matrix(0)
integrated_data@reductions[["pca"]]<-pca
integrated_data <- RunUMAP(integrated_data, reduction = "pca", dims = 1:5)
p7 <- DimPlot(integrated_data, reduction = "umap", group.by = "tech",pt.size = 1)
p8 <- DimPlot(integrated_data, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE,pt.size = 1,label.size = 8) 
plot_grid(p1, p2)
integrated_data<- FindNeighbors(integrated_data, dims = 1:10)
integrated_data <- FindClusters(integrated_data, resolution = 0.09)
adjustedRandIndex(integrated_data@meta.data$celltype,integrated_data@meta.data$seurat_clusters)

subset_size <- 0.05 #subsample to 10% of the data
subset_id <- sample.int(n = length(integrated_data@meta.data$tech), size = floor(subset_size * length(integrated_data@meta.data$tech)), replace=FALSE)
batch.estimate <- kBET(t(integrated_data@assays$RNA@data)[subset_id,], integrated_data@meta.data$tech[subset_id])

lisi_res <- lisi::compute_lisi(integrated_data@reductions$pca@cell.embeddings, integrated_data@meta.data,c('tech', 'celltype'))


#seurat
pancreas.list <- SplitObject(pbmc, split.by = "tech")
all.genes <- rownames(pbmc)
for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
  pancreas.list[[i]] <- ScaleData(pancreas.list[[i]],features = all.genes)
  pancreas.list[[i]] <- RunPCA(pancreas.list[[i]], features = VariableFeatures(object = pancreas.list[[i]]))
  #pancreas.list[[i]] <- FindNeighbors(pancreas.list[[i]], dims = 1:10)
}
#clustering

pancreas.list[[1]] <- RunUMAP(pancreas.list[[1]], reduction = "pca", dims = 1:30)
DimPlot(pancreas.list[[1]], reduction = "umap")
pancreas.list[[2]] <- RunUMAP(pancreas.list[[2]], reduction = "pca", dims = 1:30)
DimPlot(pancreas.list[[2]], reduction = "umap")

#seurat find anchor
pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, dims = 1:30,k.anchor = 10,anchor.features=genename)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30,features.to.integrate	= genename)
#pancreas.list <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30,features.to.integrate	= genename)
pancreas.anchors<-0

DefaultAssay(pancreas.integrated) <- "integrated"
#DefaultAssay(pancreas.list) <- "integrated"
#pancreas.list<- ScaleData(pancreas.list, verbose = FALSE)

# Run the standard workflow for visualization and clustering
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- FindNeighbors(pancreas.integrated, 
                                     dims = 1:10)
pancreas.integrated <- FindClusters(pancreas.integrated, resolution = 0.09)
#pancreas.integrated <- RunTSNE(pbmc.atac,dims = 1:10,check_duplicates = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:20)
DimPlot(pancreas.integrated, reduction = "umap")
p5 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech",pt.size = 1)
p6 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype", label = TRUE, 
              repel = TRUE,pt.size = 1,label.size = 8) 
p1 + p2

# impute ARI
adjustedRandIndex(pancreas.integrated@meta.data$celltype,pancreas.integrated@meta.data$seurat_clusters)

subset_size <- 0.25 #subsample to 10% of the data
subset_id <- sample.int(n = length(pancreas.integrated@meta.data$tech), size = floor(subset_size * length(pancreas.integrated@meta.data$tech)), replace=FALSE)
batch.estimate <- kBET(t(pancreas.integrated@assays$RNA@data)[subset_id,], pancreas.integrated@meta.data$tech[subset_id])

lisi_res <- lisi::compute_lisi(pancreas.integrated@reductions$pca@cell.embeddings, pancreas.integrated@meta.data,c('tech', 'celltype'))

genename<-VariableFeatures(pbmc)
batch1_0<-SubsetData(pancreas.list[[1]],cell=pancreas.list[[1]]@meta.data$celltype =="B cell")
batch1_1<-SubsetData(pancreas.list[[1]],cell=pancreas.list[[1]]@meta.data$celltype =="CD4 T cell")
batch1_2<-SubsetData(pancreas.list[[1]],cell=pancreas.list[[1]]@meta.data$celltype =="CD8 T cell")
batch1_3<-SubsetData(pancreas.list[[1]],cell=pancreas.list[[1]]@meta.data$celltype =="Monocyte_CD14")
batch1_4<-SubsetData(pancreas.list[[1]],cell=pancreas.list[[1]]@meta.data$celltype =="Monocyte_FCGR3A")
batch1_5<-SubsetData(pancreas.list[[1]],cell=pancreas.list[[1]]@meta.data$celltype =="NK cell")
batch1_6<-SubsetData(pancreas.list[[1]],cell=pancreas.list[[1]]@meta.data$celltype =="Plasmacytoid dendritic cell")

batch2_0<-SubsetData(pancreas.list[[2]],cell=pancreas.list[[2]]@meta.data$celltype =="B cell")
batch2_1<-SubsetData(pancreas.list[[2]],cell=pancreas.list[[2]]@meta.data$celltype =="CD4 T cell")
batch2_2<-SubsetData(pancreas.list[[2]],cell=pancreas.list[[2]]@meta.data$celltype =="CD8 T cell")
batch2_3<-SubsetData(pancreas.list[[2]],cell=pancreas.list[[2]]@meta.data$celltype =="Monocyte_CD14")
batch2_4<-SubsetData(pancreas.list[[2]],cell=pancreas.list[[2]]@meta.data$celltype =="Monocyte_FCGR3A")
batch2_5<-SubsetData(pancreas.list[[2]],cell=pancreas.list[[2]]@meta.data$celltype =="NK cell")
batch2_6<-SubsetData(pancreas.list[[2]],cell=pancreas.list[[2]]@meta.data$celltype =="Plasmacytoid dendritic cell")

pancreas.list<-0
#integrated
#CD4
CD4<-merge(batch1_1,batch2_1)
batch1_1<-0
batch2_1<-0
CD4<-SplitObject(CD4, split.by = "tech")
for (i in 1:length(CD4)) {
  CD4[[i]] <- NormalizeData(CD4[[i]], verbose = FALSE)
  CD4[[i]] <- FindVariableFeatures(CD4[[i]], selection.method = "vst", 
                                      nfeatures = 2000, verbose = FALSE)
  CD4[[i]] <- ScaleData(CD4[[i]],features = all.genes)
  CD4[[i]] <- RunPCA(CD4[[i]], features = VariableFeatures(object = CD4[[i]]))
  #CD4[[i]] <- RunUMAP(CD4[[i]], reduction = "pca", dims = 1:50)
  #pancreas.list[[i]] <- FindNeighbors(pancreas.list[[i]], dims = 1:10)
}

CD4.anchor <- FindIntegrationAnchors(object.list = CD4, dims = 1:30,k.anchor = 10,k.filter=200,scale = TRUE,anchor.features=genename)
CD4_integrated <- IntegrateData(anchorset = CD4.anchor, dims = 1:50,features.to.integrate	= genename)
DefaultAssay(CD4_integrated) <- "integrated"
CD4<-0
# Run the standard workflow for visualization and clustering
CD4_integrated <- ScaleData(CD4_integrated, verbose = FALSE)
CD4_integrated <- RunPCA(CD4_integrated, npcs = 50, verbose = FALSE)
CD4_integrated <- RunUMAP(CD4_integrated, reduction = "pca", dims = 1:50)
p1 <- DimPlot(CD4_integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(CD4_integrated, reduction = "umap", group.by = "celltype", label = TRUE, 
              repel = TRUE) + NoLegend()
p1 + p2

#CD8
CD8<-merge(batch1_2,batch2_2)
batch1_2<-0
batch2_2<-0
CD8<-SplitObject(CD8, split.by = "tech")
for (i in 1:length(CD8)) {
  CD8[[i]] <- NormalizeData(CD8[[i]], verbose = FALSE)
  CD8[[i]] <- FindVariableFeatures(CD8[[i]], selection.method = "vst", 
                                   nfeatures = 2000, verbose = FALSE)
  CD8[[i]] <- ScaleData(CD8[[i]],features = all.genes)
  CD8[[i]] <- RunPCA(CD8[[i]], features = VariableFeatures(object = CD8[[i]]))
  #CD4[[i]] <- RunUMAP(CD4[[i]], reduction = "pca", dims = 1:50)
  #pancreas.list[[i]] <- FindNeighbors(pancreas.list[[i]], dims = 1:10)
}

CD8.anchor <- FindIntegrationAnchors(object.list = CD8, dims = 1:30,k.anchor = 10,k.filter=200,scale = TRUE,anchor.features=genename)
CD8_integrated <- IntegrateData(anchorset = CD8.anchor, dims = 1:50,features.to.integrate	= genename)
DefaultAssay(CD8_integrated) <- "integrated"
CD8<-0

#Monocyte_CD14
Monocyte_CD14 <-merge(batch1_3,batch2_3)
batch1_3<-0
batch2_3<-0
Monocyte_CD14<-SplitObject(Monocyte_CD14, split.by = "tech")
for (i in 1:length(Monocyte_CD14)) {
  Monocyte_CD14[[i]] <- NormalizeData(Monocyte_CD14[[i]], verbose = FALSE)
  Monocyte_CD14[[i]] <- FindVariableFeatures(Monocyte_CD14[[i]], selection.method = "vst", 
                                   nfeatures = 2000, verbose = FALSE)
  Monocyte_CD14[[i]] <- ScaleData(Monocyte_CD14[[i]],features = all.genes)
  Monocyte_CD14[[i]] <- RunPCA(Monocyte_CD14[[i]], features = VariableFeatures(object = Monocyte_CD14[[i]]))
  #CD4[[i]] <- RunUMAP(CD4[[i]], reduction = "pca", dims = 1:50)
  #pancreas.list[[i]] <- FindNeighbors(pancreas.list[[i]], dims = 1:10)
}

Monocyte_CD14.anchor <- FindIntegrationAnchors(object.list = Monocyte_CD14, dims = 1:30,k.anchor = 10,k.filter=200,scale = TRUE,anchor.features=genename)
Monocyte_CD14<-0
Monocyte_CD14_integrated <- IntegrateData(anchorset = Monocyte_CD14.anchor, dims = 1:50,features.to.integrate	= genename)
DefaultAssay(Monocyte_CD14_integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
Monocyte_CD14_integrated <- ScaleData(Monocyte_CD14_integrated, verbose = FALSE)
Monocyte_CD14_integrated <- RunPCA(Monocyte_CD14_integrated, npcs = 50, verbose = FALSE)
Monocyte_CD14_integrated <- RunUMAP(Monocyte_CD14_integrated, reduction = "pca", dims = 1:50)
p1 <- DimPlot(Monocyte_CD14_integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(Monocyte_CD14_integrated, reduction = "umap", group.by = "celltype", label = TRUE, 
              repel = TRUE) + NoLegend()
p1 + p2

#B
B<-merge(batch1_0,batch2_0)
batch1_0<-0
batch2_0<-0
B<-SplitObject(B, split.by = "tech")
for (i in 1:length(B)) {
  B[[i]] <- NormalizeData(B[[i]], verbose = FALSE)
  B[[i]] <- FindVariableFeatures(B[[i]], selection.method = "vst", 
                                   nfeatures = 2000, verbose = FALSE)
  B[[i]] <- ScaleData(B[[i]],features = all.genes)
  B[[i]] <- RunPCA(B[[i]], features = VariableFeatures(object = B[[i]]))
  #CD4[[i]] <- RunUMAP(CD4[[i]], reduction = "pca", dims = 1:50)
  #pancreas.list[[i]] <- FindNeighbors(pancreas.list[[i]], dims = 1:10)
}

B.anchor <- FindIntegrationAnchors(object.list = B, dims = 1:30,k.anchor = 10,k.filter=200,scale = TRUE,anchor.features=genename)
B_integrated <- IntegrateData(anchorset = B.anchor, dims = 1:50,features.to.integrate	= genename)
B<-0
DefaultAssay(B_integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
B_integrated <- ScaleData(B_integrated, verbose = FALSE)
B_integrated <- RunPCA(B_integrated, npcs = 50, verbose = FALSE)
B_integrated <- RunUMAP(B_integrated, reduction = "pca", dims = 1:50)
p1 <- DimPlot(B_integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(B_integrated, reduction = "umap", group.by = "celltype", label = TRUE, 
              repel = TRUE) + NoLegend()
p1 + p2


#NK
NK<-merge(batch1_5,batch2_5)
batch2_5<-0
batch1_5<-0
NK<-SplitObject(NK, split.by = "tech")
for (i in 1:length(NK)) {
  NK[[i]] <- NormalizeData(NK[[i]], verbose = FALSE)
  NK[[i]] <- FindVariableFeatures(NK[[i]], selection.method = "vst", 
                                 nfeatures = 2000, verbose = FALSE)
  NK[[i]] <- ScaleData(NK[[i]],features = all.genes)
  NK[[i]] <- RunPCA(NK[[i]], features = VariableFeatures(object = NK[[i]]))
  #CD4[[i]] <- RunUMAP(CD4[[i]], reduction = "pca", dims = 1:50)
  #pancreas.list[[i]] <- FindNeighbors(pancreas.list[[i]], dims = 1:10)
}

NK.anchor <- FindIntegrationAnchors(object.list = NK, dims = 1:30,k.anchor = 10,k.filter=200,scale = TRUE,anchor.features=genename)
NK_integrated <- IntegrateData(anchorset = NK.anchor, dims = 1:50,features.to.integrate	= genename)
DefaultAssay(NK_integrated) <- "integrated"
NK<-0
# Run the standard workflow for visualization and clustering
NK_integrated <- ScaleData(NK_integrated, verbose = FALSE)
NK_integrated <- RunPCA(NK_integrated, npcs = 50, verbose = FALSE)
NK_integrated <- RunUMAP(NK_integrated, reduction = "pca", dims = 1:50)
p1 <- DimPlot(NK_integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(NK_integrated, reduction = "umap", group.by = "celltype", label = TRUE, 
              repel = TRUE) + NoLegend()
p1 + p2

#Monocyte_fcgr3a
Monocyte_fcgr3a<-merge(batch1_4,batch2_4)
batch1_4<-0
batch2_4<-0
Monocyte_fcgr3a<-SplitObject(Monocyte_fcgr3a, split.by = "tech")
for (i in 1:length(Monocyte_fcgr3a)) {
  Monocyte_fcgr3a[[i]] <- NormalizeData(Monocyte_fcgr3a[[i]], verbose = FALSE)
  Monocyte_fcgr3a[[i]] <- FindVariableFeatures(Monocyte_fcgr3a[[i]], selection.method = "vst", 
                                 nfeatures = 2000, verbose = FALSE)
  Monocyte_fcgr3a[[i]] <- ScaleData(Monocyte_fcgr3a[[i]],features = all.genes)
  Monocyte_fcgr3a[[i]] <- RunPCA(Monocyte_fcgr3a[[i]], features = VariableFeatures(object = Monocyte_fcgr3a[[i]]))
  #CD4[[i]] <- RunUMAP(CD4[[i]], reduction = "pca", dims = 1:50)
  #pancreas.list[[i]] <- FindNeighbors(pancreas.list[[i]], dims = 1:10)
}

Monocyte_fcgr3a.anchor <- FindIntegrationAnchors(object.list = Monocyte_fcgr3a, dims = 1:30,k.anchor = 30,k.filter=200,scale = TRUE,anchor.features=genename)
Monocyte_fcgr3a<-0
Monocyte_fcgr3a_integrated <- IntegrateData(anchorset = Monocyte_fcgr3a.anchor, dims = 1:50,features.to.integrate	= genename)
DefaultAssay(Monocyte_fcgr3a_integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
Monocyte_fcgr3a_integrated <- ScaleData(Monocyte_fcgr3a_integrated, verbose = FALSE)
Monocyte_fcgr3a_integrated <- RunPCA(Monocyte_fcgr3a_integrated, npcs = 50, verbose = FALSE)
Monocyte_fcgr3a_integrated <- RunUMAP(Monocyte_fcgr3a_integrated, reduction = "pca", dims = 1:50)
p1 <- DimPlot(Monocyte_fcgr3a_integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(Monocyte_fcgr3a_integrated, reduction = "umap", group.by = "celltype", label = TRUE, 
              repel = TRUE) + NoLegend()
p1 + p2



#plas
plas<-merge(batch1_6,batch2_6)
batch1_6<-0
batch2_6<-0
plas<-SplitObject(plas, split.by = "tech")
for (i in 1:length(plas)) {
  plas[[i]] <- NormalizeData(plas[[i]], verbose = FALSE)
  plas[[i]] <- FindVariableFeatures(plas[[i]], selection.method = "vst", 
                                    nfeatures = 2000, verbose = FALSE)
  plas[[i]] <- ScaleData(plas[[i]],features = all.genes)
  plas[[i]] <- RunPCA(plas[[i]], features = VariableFeatures(object = plas[[i]]))
  #CD4[[i]] <- RunUMAP(CD4[[i]], reduction = "pca", dims = 1:50)
  #pancreas.list[[i]] <- FindNeighbors(pancreas.list[[i]], dims = 1:10)
}

plas.anchor <- FindIntegrationAnchors(object.list = plas, dims = 1:30,k.anchor = 20,k.filter=50,scale = TRUE,anchor.features=genename)
plas_integrated <- IntegrateData(anchorset = plas.anchor, dims = 1:50,features.to.integrate	= genename)
DefaultAssay(plas_integrated) <- "integrated"
plas<-0
# Run the standard workflow for visualization and clustering
plas_integrated <- ScaleData(plas_integrated, verbose = FALSE)
plas_integrated <- RunPCA(plas_integrated, npcs = 50, verbose = FALSE)
plas_integrated <- RunUMAP(plas_integrated, reduction = "pca", dims = 1:50)
p1 <- DimPlot(plas_integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(plas_integrated, reduction = "umap", group.by = "celltype", label = TRUE, 
              repel = TRUE) + NoLegend()
p1 + p2

#整合anchor
anchor<-rbind(CD4.anchor@anchors,Monocyte_CD14.anchor@anchors,Monocyte_fcgr3a.anchor@anchors,NK.anchor@anchors,plas.anchor@anchors,B.anchor@anchors)
anchor<-subset(anchor,anchor$score>0.7)
anchor_new<-pancreas.anchors
anchor_new@anchors<-anchor
anchor_new<-0
integrated <- IntegrateData(anchorset = anchor_new, dims = 1:30)
anchor_new<-0
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
#integrated <- FindNeighbors(integrated, dims = 1:10)

#integrated <- FindClusters(integrated, resolution = 0.14)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(integrated, reduction = "umap", group.by = "celltype", label = TRUE, 
              repel = TRUE) + NoLegend()
p1 + p2

meta_new<-rbind(CD4_integrated@meta.data,B_integrated@meta.data,Mega_integrated@meta.data,Monocyte_CD14_integrated@meta.data,
                Monocyte_fcgr3a_integrated@meta.data,NK_integrated@meta.data,plas_integrated@meta.data)

pancreas.list<-0
#impute vectors
cellname<-c(colnames(B_integrated),colnames(CD4_integrated),colnames(CD8_integrated),colnames(Monocyte_CD14_integrated),colnames(Monocyte_fcgr3a_integrated),colnames(NK_integrated),colnames(plas_integrated))
cellname_old<-setdiff(colnames(pbmc),cellname)
data_old<-SubsetData(pancreas.integrated,cell=cellname_old)
meta_old<-data_old@meta.data
data_old<-as.matrix(data_old@assays$integrated@data)
data_old<-subset(data_old,row.names(data_old)%in%genename)


CD4_integrated<-as.matrix(CD4_integrated@assays$integrated@data)
CD8_integrated<-as.matrix(CD8_integrated@assays$integrated@data)
CD4_integrated<-CD4_integrated[order(row.names(CD4_integrated),decreasing=F),]
B_integrated<-as.matrix(B_integrated@assays$integrated@data)
B_integrated<-B_integrated[order(row.names(B_integrated),decreasing=F),]

Monocyte_CD14_integrated<-as.matrix(Monocyte_CD14_integrated@assays$integrated@data)
Monocyte_CD14_integrated<-Monocyte_CD14_integrated[order(row.names(Monocyte_CD14_integrated),decreasing=F),]
Monocyte_fcgr3a_integrated<-as.matrix(Monocyte_fcgr3a_integrated@assays$integrated@data)
Monocyte_fcgr3a_integrated<-Monocyte_fcgr3a_integrated[order(row.names(Monocyte_fcgr3a_integrated),decreasing=F),]
NK_integrated<-as.matrix(NK_integrated@assays$integrated@data)
NK_integrated<-NK_integrated[order(row.names(NK_integrated),decreasing=F),]
plas_integrated<-as.matrix(plas_integrated@assays$integrated@data)
plas_integrated<-plas_integrated[order(row.names(plas_integrated),decreasing=F),]

data_new<-cbind(CD4_integrated,CD8_integrated,B_integrated,Monocyte_CD14_integrated,Monocyte_fcgr3a_integrated,NK_integrated,plas_integrated)

#pancreas.anchors<-0
pancreas.list<-0
data_new <- CreateSeuratObject(counts = data_new,meta.data = meta_new)
all.genes <- rownames(data_new)
data_new <- FindVariableFeatures(data_new, selection.method = "vst", nfeatures = 6000, verbose = FALSE)
data_new <- ScaleData(data_new, verbose = FALSE)
data_new <- RunPCA(data_new,npcs = 50, verbose = FALSE,features = VariableFeatures(object = data_new))
data_new <- FindNeighbors(data_new, dims = 1:10)

# clustering
data_new<- FindClusters(data_new, resolution = 0.055)
data_new <- RunUMAP(data_new, reduction = "pca", dims = 1:30)
DimPlot(data_new, reduction = "umap",label = TRUE)
p1 <- DimPlot(data_new, reduction = "umap", group.by = "tech")
p2 <- DimPlot(data_new, reduction = "umap", group.by = "celltype", label = TRUE, 
              repel = TRUE) + NoLegend()
p1 + p2

adjustedRandIndex(data_new@meta.data$celltype,data_new@meta.data$seurat_clusters)

