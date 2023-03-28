## Made by Christopher M. Horn, MS
## Kielian Lab data
## Analyzing craniotomy (boneflap) model scRNA-seq data w/Seurat
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
boneflap.data <- Read10X(data.dir = '/Users/christopherhorn/Documents/Research/PhD/Kielian Lab/Data/scRNAseq/Craniotomy model- scRNA-seq Day 7 post-infection/R Analysis/Seurat/BoneFlap/Data')

## Initialize Seurat object w/raw data
boneflap <- CreateSeuratObject(counts = boneflap.data, project = 'BoneFlap_Cranio_scRNAseq', min.cells = 3, min.features = 200)

## Pre-processing and QC

## Find % of mitochondrional contamination
boneflap[['percent.mt']] <- PercentageFeatureSet(boneflap, pattern = '^mt-')

## Plot QC metrics
QC_metrics_plot <- VlnPlot(boneflap, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)

## Plot feature-feature relationships
count2mt <- FeatureScatter(boneflap, feature1 = 'nCount_RNA', feature2 = 'percent.mt')
count2feature <- FeatureScatter(boneflap, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')

## Select cells w/at least 200 features & less than 4500 features & percent.mt < 5
boneflap <- subset(boneflap, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

## Normalize the data
boneflap <- NormalizeData(boneflap, normalization.method = 'LogNormalize', scale.factor = 10000)

## Find highly variable genes
boneflap <- FindVariableFeatures(boneflap, selection.method = 'vst', nfeatures = 2000)

## Find top 10 most variable
top10 <- head(VariableFeatures(boneflap), 10)

## Plot the vaible genes with and without labels
vary_no_label <- VariableFeaturePlot(boneflap)
vary_w_label <- LabelPoints(plot = vary_no_label, points = top10, repel = T)

## Scale the data & regress out mitochondrional contamination
all.genes <- rownames(boneflap)
boneflap <- ScaleData(boneflap, features = all.genes, vars.to.regress = 'percent.mt')

## Perform a linear dimensional reduction
boneflap <- RunPCA(boneflap, features = VariableFeatures(object = boneflap))
BoneFlap_Cranio_PCA <- DimPlot(boneflap, reduction = 'pca')

## Find the dimensionality of the dataset

## Jackstraw method
boneflap <- JackStraw(boneflap, num.replicate = 1000, dims = 50)
boneflap <- ScoreJackStraw(boneflap, dims = 1:50)
js_plot <- JackStrawPlot(boneflap, dims = 1:50)

## Elbow plot method
e_plot <- ElbowPlot(boneflap)

## Clustering the cells
boneflap <- FindNeighbors(boneflap, dims = 1:20)
boneflap <- FindClusters(boneflap, resolution = 0.5)

## Non-linear dimensional reduction
## UMAP version
boneflap <- RunUMAP(boneflap, dims = 1:20)
BoneFlap_Cranio_UMAP <- DimPlot(boneflap, reduction = 'umap', label = T) + labs(title = 'Craniotomy Model - Bone Flap')
BoneFlap_Cranio_UMAP_NoLegend <- DimPlot(boneflap, reduction = 'umap', label = T) + labs(title = 'Craniotomy Model - Bone Flap') + NoLegend()

## TSNE version
## boneflap <- RunTSNE(boneflap, dims = 1:20, tsne.method = 'Flt-SNE')

## Save R data to skip computation time/effort
saveRDS(boneflap, file = '/Users/christopherhorn/Documents/Research/PhD/Kielian Lab/Data/scRNAseq/Craniotomy model- scRNA-seq Day 7 post-infection/R Analysis/Seurat/BoneFlap/CranioBoneFlap/BoneFlap Cranio scRNA-seq Data.rds')

## Find differentially expressed markers
boneflap.markers <- FindAllMarkers(boneflap, min.pct = 0.25, logfc.threshold = 0.25, test.use = 'MAST')
boneflap.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

## Output differentially expressed markers in a .csv file
write.csv(boneflap.markers %>% group_by(cluster), file = '/Users/christopherhorn/Documents/Research/PhD/Kielian Lab/Data/scRNAseq/Craniotomy model- scRNA-seq Day 7 post-infection/R Analysis/Seurat/BoneFlap/Output/BoneFlap Cranio Differentially Expressed Markers.csv')

## Remove temp objects
rm(top10, all.genes)



## Annotate the cells w/SingleR

## Store the reference dataset
immgen.se <- ImmGenData()

## Convert Seurat object into a SingleCellExperiment object
boneflap_sce <- as.SingleCellExperiment(boneflap)

## Find the commonalities between the test and reference sets for efficiency
commonGenes <- intersect(rownames(boneflap_sce), rownames(immgen.se))
immgen.se <- immgen.se[commonGenes,]
boneflap_sce <- boneflap_sce[commonGenes,]

## Start the prediction

## w/Main labels
pred.boneflap <- SingleR(test = boneflap_sce, ref = immgen.se, labels = immgen.se$label.main)

## w/Fine labels
## pred.boneflap <- SingleR(test = boneflap_sce, ref = immgen.se, labels = immgen.se$label.fine)

## Plot annotation confidence
BoneFlap_Cranio_AnnotateHeatMap <- plotScoreHeatmap(pred.boneflap, show.pruned = T)

## Add annotation data to the original Seurat object
boneflap[['SingleR.labels']] <- pred.boneflap$pruned.labels

## Store the oringal Ident information and switch out to SingleR labels
boneflap$store <- Idents(boneflap)
Idents(boneflap) <- boneflap$SingleR.labels

## Plot the annotations on UMAP

## w/Legend
BoneFlap_Cranio_Annotate <- DimPlot(boneflap, reduction = 'umap', label = T) + labs(title = 'Craniotomy Model - Bone Flap')

## w/o Legend
BoneFlap_Cranio_Annotate_NoLegend <- DimPlot(boneflap, reduction = 'umap', label = T) + labs(title = 'Craniotomy Model - Bone Flap') + NoLegend()

## Write a file containing all of the metadata
write.csv(boneflap@meta.data %>% select(-store) %>% arrange(seurat_clusters), file = '/Users/christopherhorn/Documents/Research/PhD/Kielian Lab/Data/scRNAseq/Craniotomy model- scRNA-seq Day 7 post-infection/R Analysis/Seurat/BoneFlap/Output/Bone Flap Cranio Annotations.csv')

## Plot Annotated UMAP alongside Clustered UMAP
BoneFlap_Cranio_Annotate_vs_Cluster <- BoneFlap_Cranio_UMAP_NoLegend + BoneFlap_Cranio_Annotate_NoLegend

## Remove temp objects
rm(immgen.se, boneflap_sce, commonGenes, pred.boneflap)





## TEST AREA
markers.use <- subset(boneflap.markers, avg_logFC > 1 & p_val_adj < 0.05)$gene # Grabbing markers for heatmap
DoHeatmap(subset(boneflap, downsample = 100), features = markers.use, size = 4, angle = 90) + scale_fill_gradientn(colors = c('blue', 'white', 'red')) + NoLegend() + theme(text = element_text(size = 0)) # Making heatmap
new.cluster.ids <- c('Granulocytes I', 'Granulocytes II', 'Granulocytes III', 'Granulocytes IV', 'Granulocytes V', 'Granulocytes VI', 'Granulocytes VII', 'Granulocytes VIII', '8', 'Granulocytes IX') # Preparing new cluster labels
names(new.cluster.ids) <- levels(boneflap)
boneflap <- RenameIdents(boneflap, new.cluster.ids)
boneflap.gran <- subset(boneflap, idents = c('Granulocytes I', 'Granulocytes II', 'Granulocytes III', 'Granulocytes IV', 'Granulocytes V', 'Granulocytes VI', 'Granulocytes VII', 'Granulocytes VIII', 'Granulocytes IX'))






## IPA Canonical Pathways
data$Term <- factor(data$Term, levels = rev(unique(data$Term))) ## Preserves the order of the terms rather than auto-alphabetizing
ggplot(data, aes(x = Term, y = Annotation)) +
  geom_point(shape = 21, color = 'black', (aes(fill = Bin, size = Count))) +
  coord_flip() +
  theme_minimal() +
  scale_fill_viridis_b('p-value', labels = c('0.05', '0.01', '0.001'), option = 'C', direction = 1) +
  ggtitle(label = 'IPA Canonical Pathways in Bone Flap') +
  theme(legend.position = 'right', axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, face = 'bold'), axis.text.y = element_text(face = 'bold'), legend.title = element_text(face = 'bold'), plot.title = element_text(size = rel(2)))





## Gene Dotplots
dot <- DotPlot(boneflap, features = features)

dot.data <- dot$data

ggplot(dot.data, aes(x = id, y = features.plot)) +
  geom_point(shape = 21, color = 'black', (aes(fill = avg.exp.scaled, size = pct.exp))) +
  theme_minimal() +
  scale_fill_gradientn('Scaled Expression', labels = c('-1', '0', '1'), colors = c('blue', 'white', 'red'), breaks = c(-1, 0, 1)) +
  theme(legend.position = 'bottom', panel.grid.major = element_blank(), axis.title = element_blank(), axis.text = element_text(face = 'bold', size = 14), legend.text = element_text(face = 'bold', size = 12), legend.title = element_text(face = 'bold', size = 12, vjust = 1), axis.text.x = element_text(angle = 90, hjust = 0)) +
  scale_x_discrete(position = 'top') +
  labs(size = 'Percent Expressed')









