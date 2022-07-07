## Made by Christopher M. Horn, MS
## Kielian Lab data
## Analyzing craniotomy (blood) model scRNA-seq data w/Seurat
## 2020-03-23

## Load in packages
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SingleR)
library(SingleCellExperiment)
library(MAST)


## Cluster the cells w/Seurat

## Load in the data
blood.data <- Read10X(data.dir = '/Users/christopherhorn/Documents/Research/PhD/Kielian Lab/Data/scRNAseq/Craniotomy model- scRNA-seq Day 7 post-infection/R Analysis/Seurat/Blood/Data')

## Initialize Seurat object w/raw data
blood <- CreateSeuratObject(counts = blood.data, project = 'Blood_Cranio_scRNAseq', min.cells = 3, min.features = 200)

## Pre-processing and QC

## Find % of mitochondrional contamination
blood[['percent.mt']] <- PercentageFeatureSet(blood, pattern = '^mt-')

## Plot QC metrics
QC_metrics_plot <- VlnPlot(blood, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)

## Plot feature-feature relationships
count2mt <- FeatureScatter(blood, feature1 = 'nCount_RNA', feature2 = 'percent.mt')
count2feature <- FeatureScatter(blood, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')

## Select cells w/at least 200 features & less than 4500 features & percent.mt < 5
blood <- subset(blood, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 5)

## Normalize the data
blood <- NormalizeData(blood, normalization.method = 'LogNormalize', scale.factor = 10000)

## Find highly variable genes
blood <- FindVariableFeatures(blood, selection.method = 'vst', nfeatures = 2000)

## Find top 10 most variable
top10 <- head(VariableFeatures(blood), 10)

## Plot the vaible genes with and without labels
vary_no_label <- VariableFeaturePlot(blood)
vary_w_label <- LabelPoints(plot = vary_no_label, points = top10, repel = T)

## Scale the data & regress out mitochondrional contamination
all.genes <- rownames(blood)
blood <- ScaleData(blood, features = all.genes, vars.to.regress = 'percent.mt')

## Perform a linear dimensional reduction
blood <- RunPCA(blood, features = VariableFeatures(object = blood))
Blood_Cranio_PCA <- DimPlot(blood, reduction = 'pca')

## Find the dimensionality of the dataset

## Jackstraw method
blood <- JackStraw(blood, num.replicate = 1000, dims = 50)
blood <- ScoreJackStraw(blood, dims = 1:50)
js_plot <- JackStrawPlot(blood, dims = 1:50)

## Elbow plot method
e_plot <- ElbowPlot(blood)

## Clustering the cells
blood <- FindNeighbors(blood, dims = 1:20)
blood <- FindClusters(blood, resolution = 0.5)

## Non-linear dimensional reduction
## UMAP version
blood <- RunUMAP(blood, dims = 1:20)
Blood_Cranio_UMAP <- DimPlot(blood, reduction = 'umap', label = T) + labs(title = 'Craniotomy Model - Blood')
Blood_Cranio_UMAP_NoLegend <- DimPlot(blood, reduction = 'umap', label = T) + labs(title = 'Craniotomy Model - Blood') + NoLegend()

## TSNE version
## blood <- RunTSNE(blood, dims = 1:20, tsne.method = 'Flt-SNE')

## Save R data to skip computation time/effort
saveRDS(blood, file = '/Users/christopherhorn/Documents/Research/PhD/Kielian Lab/Data/scRNAseq/Craniotomy model- scRNA-seq Day 7 post-infection/R Analysis/Seurat/Blood/CranioBlood/Blood Cranio scRNA-seq Data.rds')

## Find differentially expressed markers
blood.markers <- FindAllMarkers(blood, min.pct = 0.25, logfc.threshold = 0.25, test.use = 'MAST')
blood.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

## Output differentially expressed markers in a .csv file
write.csv(blood.markers %>% group_by(cluster), file = '/Users/christopherhorn/Documents/Research/PhD/Kielian Lab/Data/scRNAseq/Craniotomy model- scRNA-seq Day 7 post-infection/R Analysis/Seurat/Blood/Output/Blood Cranio Differentially Expressed Markers.csv')

## Remove temp objects
rm(top10, all.genes)



## Annotate the cells w/SingleR

## Store the reference dataset
immgen.se <- ImmGenData()

## Convert Seurat object into a SingleCellExperiment object
blood_sce <- as.SingleCellExperiment(blood)

## Find the commonalities between the test and reference sets for efficiency
commonGenes <- intersect(rownames(blood_sce), rownames(immgen.se))
immgen.se <- immgen.se[commonGenes,]
blood_sce <- blood_sce[commonGenes,]

## Start the prediction

## w/Main labels
pred.blood <- SingleR(test = blood_sce, ref = immgen.se, labels = immgen.se$label.main)

## Plot annotation confidence
Blood_Cranio_AnnotateHeatMap <- plotScoreHeatmap(pred.blood, show.pruned = T)

## Add annotation data to the original Seurat object
blood[['SingleR.labels']] <- pred.blood$pruned.labels

## Store the oringal Ident information and switch out to SingleR labels
blood$store <- Idents(blood)
Idents(blood) <- blood$SingleR.labels

## Plot the annotations on UMAP

## w/Legend
Blood_Cranio_Annotate <- DimPlot(blood, reduction = 'umap', label = T) + labs(title = 'Craniotomy Model - Blood')

## w/o Legend
Blood_Cranio_Annotate_NoLegend <- DimPlot(blood, reduction = 'umap', label = T) + labs(title = 'Craniotomy Model - Blood') + NoLegend()

## Write a file containing all of the metadata
write.csv(blood@meta.data %>% select(-store) %>% arrange(seurat_clusters), file = '/Users/christopherhorn/Documents/Research/PhD/Kielian Lab/Data/scRNAseq/Craniotomy model- scRNA-seq Day 7 post-infection/R Analysis/Seurat/Blood/Output/Blood Cranio Annotations.csv')

## Plot Annotated UMAP alongside Clustered UMAP
Blood_Cranio_Annotate_vs_Cluster <- Blood_Cranio_UMAP_NoLegend + Blood_Cranio_Annotate_NoLegend

## Remove temp objects
rm(immgen.se, blood_sce, commonGenes, pred.blood)








## TEST AREA
markers.use <- subset(blood.markers, avg_logFC > 1 & p_val_adj < 0.05)$gene # Grabbing markers for heatmap
DoHeatmap(subset(blood, downsample = 100), features = markers.use, size = 4, angle = 90) + scale_fill_gradientn(colors = c('blue', 'white', 'red')) + NoLegend() + theme(text = element_text(size = 0)) # Making heatmap
new.cluster.ids <- c('Granulocytes', 'NK Cells I', 'B Cells I', 'T Cells I', 'T Cells II', '5', 'Monocytes I', 'NKT/T Cells', 'Monocytes II', '9', 'B Cells II', 'NK/NKT/T Cells', 'Monocytes III', 'NK Cells II', '14', 'Basophils', '16') # Preparing new cluster labels
names(new.cluster.ids) <- levels(blood)
blood <- RenameIdents(blood, new.cluster.ids)









