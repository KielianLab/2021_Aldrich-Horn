## Made by Christopher M. Horn, MS
## Kielian Lab data
## Analyzing craniotomy (galea) model scRNA-seq data w/Seurat
## 2020-03-24

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
galea.data <- Read10X(data.dir = '/Users/christopherhorn/Documents/Research/PhD/Kielian Lab/Data/scRNAseq/Craniotomy model- scRNA-seq Day 7 post-infection/R Analysis/Seurat/Galea/Data')

## Initialize Seurat object w/raw data
galea <- CreateSeuratObject(counts = galea.data, project = 'Galea_Cranio_scRNAseq', min.cells = 3, min.features = 200)

## Pre-processing and QC

## Find % of mitochondrional contamination
galea[['percent.mt']] <- PercentageFeatureSet(galea, pattern = '^mt-')

## Plot QC metrics
QC_metrics_plot <- VlnPlot(galea, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)

## Plot feature-feature relationships
count2mt <- FeatureScatter(galea, feature1 = 'nCount_RNA', feature2 = 'percent.mt')
count2feature <- FeatureScatter(galea, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')

## Select cells w/at least 200 features & less than 4500 features & percent.mt < 5
galea <- subset(galea, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 5)

## Normalize the data
galea <- NormalizeData(galea, normalization.method = 'LogNormalize', scale.factor = 10000)

## Find highly variable genes
galea <- FindVariableFeatures(galea, selection.method = 'vst', nfeatures = 2000)

## Find top 10 most variable
top10 <- head(VariableFeatures(galea), 10)

## Plot the vaible genes with and without labels
vary_no_label <- VariableFeaturePlot(galea)
vary_w_label <- LabelPoints(plot = vary_no_label, points = top10, repel = T)

## Scale the data & regress out mitochondrional contamination
all.genes <- rownames(galea)
galea <- ScaleData(galea, features = all.genes, vars.to.regress = 'percent.mt')

## Perform a linear dimensional reduction
galea <- RunPCA(galea, features = VariableFeatures(object = galea))
Galea_Cranio_PCA <- DimPlot(galea, reduction = 'pca')

## Find the dimensionality of the dataset

## Jackstraw method
galea <- JackStraw(galea, num.replicate = 1000, dims = 50)
galea <- ScoreJackStraw(galea, dims = 1:50)
js_plot <- JackStrawPlot(galea, dims = 1:50)

## Elbow plot method
e_plot <- ElbowPlot(galea)

## Clustering the cells
galea <- FindNeighbors(galea, dims = 1:30)
galea <- FindClusters(galea, resolution = 0.5)

## Non-linear dimensional reduction
## UMAP version
galea <- RunUMAP(galea, dims = 1:30)
Galea_Cranio_UMAP <- DimPlot(galea, reduction = 'umap', label = T) + labs(title = 'Craniotomy Model - Galea')
Galea_Cranio_UMAP_NoLegend <- DimPlot(galea, reduction = 'umap', label = T) + labs(title = 'Craniotomy Model - Galea') + NoLegend()

## TSNE version
## galea <- RunTSNE(galea, dims = 1:30, tsne.method = 'Flt-SNE')

## Save R data to skip computation time/effort
saveRDS(galea, file = '/Users/christopherhorn/Documents/Research/PhD/Kielian Lab/Data/scRNAseq/Craniotomy model- scRNA-seq Day 7 post-infection/R Analysis/Seurat/Galea/CranioGalea/Galea Cranio scRNA-seq Data.rds')

## Find differentially expressed markers
galea.markers <- FindAllMarkers(galea, min.pct = 0.25, logfc.threshold = 0.25, test.use = 'MAST')
galea.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

## Output differentially expressed markers in a .csv file
write.csv(galea.markers %>% group_by(cluster), file = '/Users/christopherhorn/Documents/Research/PhD/Kielian Lab/Data/scRNAseq/Craniotomy model- scRNA-seq Day 7 post-infection/R Analysis/Seurat/Galea/Output/Galea Cranio Differentially Expressed Markers.csv')

## Remove temp objects
rm(top10, all.genes)



## Annotate the cells w/SingleR

## Store the reference dataset
immgen.se <- ImmGenData()

## Convert Seurat object into a SingleCellExperiment object
galea_sce <- as.SingleCellExperiment(galea)

## Find the commonalities between the test and reference sets for efficiency
commonGenes <- intersect(rownames(galea_sce), rownames(immgen.se))
immgen.se <- immgen.se[commonGenes,]
galea_sce <- galea_sce[commonGenes,]

## Start the prediction

## w/Main labels
pred.galea <- SingleR(test = galea_sce, ref = immgen.se, labels = immgen.se$label.main)

## Plot annotation confidence
Galea_Cranio_AnnotateHeatMap <- plotScoreHeatmap(pred.galea, show.pruned = T)

## Add annotation data to the original Seurat object
galea[['SingleR.labels']] <- pred.galea$pruned.labels

## Store the oringal Ident information and switch out to SingleR labels
galea$store <- Idents(galea)
Idents(galea) <- galea$SingleR.labels

## Plot the annotations on UMAP

## w/Legend
Galea_Cranio_Annotate <- DimPlot(galea, reduction = 'umap', label = T) + labs(title = 'Craniotomy Model - Galea')

## w/o Legend
Galea_Cranio_Annotate_NoLegend <- DimPlot(galea, reduction = 'umap', label = T) + labs(title = 'Craniotomy Model - Galea') + NoLegend()

## Write a file containing all of the metadata
write.csv(galea@meta.data %>% select(-store) %>% arrange(seurat_clusters), file = '/Users/christopherhorn/Documents/Research/PhD/Kielian Lab/Data/scRNAseq/Craniotomy model- scRNA-seq Day 7 post-infection/R Analysis/Seurat/Galea/Output/Galea Cranio Annotations.csv')

## Plot Annotated UMAP alongside Clustered UMAP
Galea_Cranio_Annotate_vs_Cluster <- Galea_Cranio_UMAP_NoLegend + Galea_Cranio_Annotate_NoLegend

## Remove temp objects
rm(immgen.se, galea_sce, commonGenes, pred.galea)







## TEST AREA
markers.use <- subset(galea.markers, avg_logFC > 1 & p_val_adj < 0.05)$gene # Grabbing markers for heatmap
DoHeatmap(subset(galea, downsample = 100), features = markers.use, size = 4, angle = 90) + scale_fill_gradientn(colors = c('blue', 'white', 'red')) + NoLegend() + theme(text = element_text(size = 0)) # Making heatmap
new.cluster.ids <- c('Granulocytes I', 'Granulocytes II', 'Granulocytes III', 'Granulocytes IV', 'Granulocytes V', 'Monocytes/Macrophages', 'NK/NKT/T Cells', 'Granulocytes VI', 'Granulocytes VII', 'Dendritic Cells', 'ILC/NKT Cells') # Preparing new cluster labels
names(new.cluster.ids) <- levels(galea)
galea <- RenameIdents(galea, new.cluster.ids)







## IPA Canonical Pathways
data$Term <- factor(data$Term, levels = rev(unique(data$Term))) ## Preserves the order of the terms rather than auto-alphabetizing
ggplot(data, aes(x = Term, y = Annotation)) +
  geom_point(shape = 21, color = 'black', (aes(fill = Bin, size = Count))) +
  coord_flip() +
  theme_minimal() +
  scale_fill_viridis_b('p-value', labels = c('0.05', '0.01', '0.001'), option = 'C', direction = 1) +
  ggtitle(label = 'IPA Canonical Pathways in Galea') +
  theme(legend.position = 'right', axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, face = 'bold'), axis.text.y = element_text(face = 'bold'), legend.title = element_text(face = 'bold'), plot.title = element_text(size = rel(2)))




## Gene Dotplots
dot <- DotPlot(galea, features = features)

dot.data <- dot$data

ggplot(dot.data, aes(x = id, y = features.plot)) +
  geom_point(shape = 21, color = 'black', (aes(fill = avg.exp.scaled, size = pct.exp))) +
  theme_minimal() +
  scale_fill_gradientn('Scaled Expression', labels = c('-1', '0', '1', '2'), colors = c('blue', 'white', 'red'), breaks = c(-1, 0, 1, 2)) +
  theme(legend.position = 'bottom', panel.grid.major = element_blank(), axis.title = element_blank(), axis.text = element_text(face = 'bold', size = 14), legend.text = element_text(face = 'bold', size = 12), legend.title = element_text(face = 'bold', size = 12, vjust = 1), axis.text.x = element_text(angle = 90, hjust = 0)) +
  scale_x_discrete(position = 'top') +
  labs(size = 'Percent Expressed')









