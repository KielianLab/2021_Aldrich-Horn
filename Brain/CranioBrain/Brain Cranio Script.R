## Made by Christopher M. Horn, MS
## Kielian Lab data
## Analyzing craniotomy (brain) model scRNA-seq data w/Seurat
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
brain.data <- Read10X(data.dir = '/Users/christopherhorn/Documents/Research/PhD/Kielian Lab/Data/scRNAseq/Craniotomy model- scRNA-seq Day 7 post-infection/R Analysis/Seurat/Brain/Data')

## Initialize Seurat object w/raw data
brain <- CreateSeuratObject(counts = brain.data, project = 'Brain_Cranio_scRNAseq', min.cells = 3, min.features = 200)

## Pre-processing and QC

## Find % of mitochondrional contamination
brain[['percent.mt']] <- PercentageFeatureSet(brain, pattern = '^mt-')

## Plot QC metrics
QC_metrics_plot <- VlnPlot(brain, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)

## Plot feature-feature relationships
count2mt <- FeatureScatter(brain, feature1 = 'nCount_RNA', feature2 = 'percent.mt')
count2feature <- FeatureScatter(brain, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')

## Select cells w/at least 200 features & less than 4500 features & percent.mt < 5
brain <- subset(brain, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)

## Normalize the data
brain <- NormalizeData(brain, normalization.method = 'LogNormalize', scale.factor = 10000)

## Find highly variable genes
brain <- FindVariableFeatures(brain, selection.method = 'vst', nfeatures = 2000)

## Find top 10 most variable
top10 <- head(VariableFeatures(brain), 10)

## Plot the vaible genes with and without labels
vary_no_label <- VariableFeaturePlot(brain)
vary_w_label <- LabelPoints(plot = vary_no_label, points = top10, repel = T)

## Scale the data & regress out mitochondrional contamination
all.genes <- rownames(brain)
brain <- ScaleData(brain, features = all.genes, vars.to.regress = 'percent.mt')

## Perform a linear dimensional reduction
brain <- RunPCA(brain, features = VariableFeatures(object = brain))
Brain_Cranio_PCA <- DimPlot(brain, reduction = 'pca')

## Find the dimensionality of the dataset

## Jackstraw method
brain <- JackStraw(brain, num.replicate = 100, dims = 50)
brain <- ScoreJackStraw(brain, dims = 1:50)
js_plot <- JackStrawPlot(brain, dims = 1:50)

## Elbow plot method
e_plot <- ElbowPlot(brain)

## Clustering the cells
brain <- FindNeighbors(brain, dims = 1:30)
brain <- FindClusters(brain, resolution = 0.5)

## Non-linear dimensional reduction
## UMAP version
brain <- RunUMAP(brain, dims = 1:30)
Brain_Cranio_UMAP <- DimPlot(brain, reduction = 'umap', label = T) + labs(title = 'Craniotomy Model - Brain')
Brain_Cranio_UMAP_NoLegend <- DimPlot(brain, reduction = 'umap', label = T) + labs(title = 'Craniotomy Model - Brain') + NoLegend()

## TSNE version
## brain <- RunTSNE(brain, dims = 1:30, tsne.method = 'Flt-SNE')

## Save R data to skip computation time/effort
saveRDS(brain, file = '/Users/christopherhorn/Documents/Research/PhD/Kielian Lab/Data/scRNAseq/Craniotomy model- scRNA-seq Day 7 post-infection/R Analysis/Seurat/Brain/CranioBrain/Brain Cranio scRNA-seq Data.rds')

## Find differentially expressed markers
brain.markers <- FindAllMarkers(brain, min.pct = 0.25, logfc.threshold = 0.25, test.use = 'MAST')
brain.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

## Output differentially expressed markers in a .csv file
write.csv(brain.markers %>% group_by(cluster), file = '/Users/christopherhorn/Documents/Research/PhD/Kielian Lab/Data/scRNAseq/Craniotomy model- scRNA-seq Day 7 post-infection/R Analysis/Seurat/Brain/Output/Brain Cranio Differentially Expressed Markers.csv')

## Remove temp objects
rm(top10, all.genes)



## Annotate the cells w/SingleR

## Store the reference dataset
immgen.se <- ImmGenData()

## Convert Seurat object into a SingleCellExperiment object
brain_sce <- as.SingleCellExperiment(brain)

## Find the commonalities between the test and reference sets for efficiency
commonGenes <- intersect(rownames(brain_sce), rownames(immgen.se))
immgen.se <- immgen.se[commonGenes,]
brain_sce <- brain_sce[commonGenes,]

## Start the prediction

## w/Main labels
pred.brain <- SingleR(test = brain_sce, ref = immgen.se, labels = immgen.se$label.main)

## Plot annotation confidence
Brain_Cranio_AnnotateHeatMap <- plotScoreHeatmap(pred.brain, show.pruned = T)

## Add annotation data to the original Seurat object
brain[['SingleR.labels']] <- pred.brain$pruned.labels

## Store the oringal Ident information and switch out to SingleR labels
brain$store <- Idents(brain)
Idents(brain) <- brain$SingleR.labels

## Plot the annotations on UMAP

## w/Legend
Brain_Cranio_Annotate <- DimPlot(brain, reduction = 'umap', label = T) + labs(title = 'Craniotomy Model - Brain')

## w/o Legend
Brain_Cranio_Annotate_NoLegend <- DimPlot(brain, reduction = 'umap', label = T) + labs(title = 'Craniotomy Model - Brain') + NoLegend()

## Write a file containing all of the metadata
write.csv(brain@meta.data %>% select(-store) %>% arrange(seurat_clusters), file = '/Users/christopherhorn/Documents/Research/PhD/Kielian Lab/Data/scRNAseq/Craniotomy model- scRNA-seq Day 7 post-infection/R Analysis/Seurat/Brain/Output/Brain Cranio Annotations.csv')

## Remove temp objects
rm(immgen.se, brain_sce, commonGenes, pred.brain)



## Create 3D plots

plot.3d <- brain
## plot.3d <- brain.sub

# Re-run UMAPs that you have accurate calculations for all UMAP(s)
plot.3d <- RunUMAP(plot.3d, dims = 1:10, n.components = 3L)

# Extract tSNE information from Seurat Object
umap_1 <- plot.3d[["umap"]]@cell.embeddings[,1]
umap_2 <- plot.3d[["umap"]]@cell.embeddings[,2]
umap_3 <- plot.3d[["umap"]]@cell.embeddings[,3]

# Visualize what headings are called so that you can extract them to form a dataframe
Embeddings(object = plot.3d, reduction = "umap")

# Prepare a dataframe for cell plotting
plot.data <- FetchData(object = plot.3d, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "seurat_clusters", 'SingleR.labels'))

# Make a column of row name identities (these will be your cell/barcode names)
plot.data$label <- paste(rownames(plot.data))

# Plot your data
plot_ly(data = plot.data, 
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
        color = ~seurat_clusters, 
        type = "scatter3d", 
        mode = "markers", 
        marker = list(size = 5, width=2), # controls size of points
        text=~SingleR.labels,
        hoverinfo="text")

## 3D plot w/expression levels of a gene

goi <- "Nfkbia" ## Enter gene of interest here
plotting.data <- FetchData(object = plot.3d, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "Expression"=goi), slot = 'data')
plotting.data$seurat_clusters <- plot.3d$seurat_clusters
plotting.data$SingleR.labels <- plot.3d$SingleR.labels

# Say you want change the scale, so that every cell having an expression >1 will be one color
# Basically, you are re-adjusting the scale here, so that any cell having a certain expression will light up on your 3D plot

# First make another column in your dataframe, where all values above 1 are re-assigned a value of 1
# This information is stored in the 'Expression' column of your dataframe
# Cutoff <- 2
Cutoff <- quantile(plotting.data[,goi], probs = .95)
plotting.data$"Expr" <- ifelse(test = plotting.data[,goi] <Cutoff, yes = plotting.data[,goi], no = Cutoff)

# Add the label column, so that now the column has 'cellname-its expression value'
plotting.data$label <- paste(rownames(plotting.data)," - ", plotting.data[,goi], sep="")

# Plot your data
plot_ly(data = plotting.data,
        # name = goi,
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
        color = ~Expr,
        opacity = .5,
        colors = c('darkgrey', 'red'), 
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 5, width=2), 
        text=~SingleR.labels,
        hoverinfo="text"
) %>%layout(title=goi)





## TEST AREA
markers.use <- subset(brain.markers, avg_logFC > 1 & p_val_adj < 0.05)$gene # Grabbing markers for heatmap
DoHeatmap(subset(brain, downsample = 100), features = markers.use, size = 4, angle = 90) + scale_fill_gradientn(colors = c('blue', 'white', 'red')) + NoLegend() + theme(text = element_text(size = 0)) # Making heatmap
new.cluster.ids <- c('Microglia I', 'Monocytes/Macrophages', 'Microglia II', 'T/NKT Cells I', 'Granulocytes I', 'gd T Cells', 'Granulocytes II', 'NK/ILC Cells', 'Microglia III', 'Dendritic Cells I', 'Microglia IV', 'NKT/ILC Cells', 'T/NKT Cells II', '13', 'B Cells', 'Dendritic Cells II', 'Monocytes', 'Pro/B Cells', '18', 'Basophils/Eosinophils', 'Granulocytes III') # Preparing new cluster labels
names(new.cluster.ids) <- levels(brain)
brain <- RenameIdents(brain, new.cluster.ids)
microglia <- subset(brain, idents = c('Microglia I', 'Microglia II', 'Microglia III', 'Microglia IV'))
DoHeatmap(microglia, features = markers.use, size = 4, angle = 90) + scale_fill_gradientn(colors = c('blue', 'white', 'red')) + NoLegend() + theme(text = element_text(size = 0)) # Making heatmap
granulocytes <- subset(brain, idents = c('Granulocytes I', 'Granulocytes II', 'Granulocytes III'))
DoHeatmap(granulocytes, features = markers.use, size = 4, angle = 90) + scale_fill_gradientn(colors = c('blue', 'white', 'red')) + NoLegend() + theme(text = element_text(size = 0)) # Making heatmap




## IPA Canonical Pathways
data$Term <- factor(data$Term, levels = rev(unique(data$Term))) ## Preserves the order of the terms rather than auto-alphabetizing
ggplot(data, aes(x = Term, y = Annotation)) +
        geom_point(shape = 21, color = 'black', (aes(fill = Bin, size = Count))) +
        coord_flip() +
        theme_minimal() +
        scale_fill_viridis_b('p-value', labels = c('0.01', '0.001', '0.0001'), option = 'C', direction = 1) +
        ggtitle(label = 'IPA Canonical Pathways in Brain') +
        theme(legend.position = 'right', axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, face = 'bold'), axis.text.y = element_text(face = 'bold'), legend.title = element_text(face = 'bold'), plot.title = element_text(size = rel(2)))





## Gene Dotplots
dot <- DotPlot(microglia, features = features)

dot.data <- dot$data

ggplot(dot.data, aes(x = id, y = features.plot)) +
        geom_point(shape = 21, color = 'black', (aes(fill = avg.exp.scaled, size = pct.exp))) +
        theme_minimal() +
        scale_fill_gradientn('Scaled Expression', labels = c('-1', '0', '1'), colors = c('blue', 'white', 'red'), breaks = c(-1, 0, 1)) +
        theme(legend.position = 'bottom', panel.grid.major = element_blank(), axis.title = element_blank(), axis.text = element_text(face = 'bold', size = 14), legend.text = element_text(face = 'bold', size = 12), legend.title = element_text(face = 'bold', size = 12, vjust = 1), axis.text.x = element_text(angle = 90, hjust = 0)) +
        scale_x_discrete(position = 'top') +
        labs(size = 'Percent Expressed')




ggplot(micro.avg, aes(x = log10(`Microglia I`), y = log10(`Microglia II`))) +
        geom_point(shape = 21, color = 'black', fill = ifelse(abs(log10(micro.avg$`Microglia I`) - log10(micro.avg$`Microglia II`)) > 1 & (log10(micro.avg$`Microglia I`) > 0 | log10(micro.avg$`Microglia II`) > 0), 'red', 'white')) +
        geom_text_repel(aes(label = (ifelse(abs(log10(micro.avg$`Microglia I`) - log10(micro.avg$`Microglia II`)) > 1 & (log10(micro.avg$`Microglia I`) > 0 | log10(micro.avg$`Microglia II`) > 0), micro.avg$gene, '')))) +
        theme_classic() +
        labs(x = bquote('log(Microglia I Expression)'), y = 'log(Microglia II Expression)') +
        geom_abline(intercept = 0, slope = 1) +
        geom_abline(intercept = 1, slope = 1, linetype = 'dashed') +
        geom_abline(intercept = -1, slope = 1, linetype = 'dashed') +
        ggtitle('Microglia I vs Microglia II') +
        geom_vline(xintercept = 0, color = 'light grey', linetype = 'dashed') +
        geom_hline(yintercept = 0, color = 'light grey', linetype = 'dashed')




m3m4 <- micro.avg
m3m4 <- m3m4 %>% select('Microglia III', 'Microglia IV', gene)
m3m4 <- m3m4 %>% mutate(m3_log = log10(`Microglia III`))
m3m4 <- m3m4 %>% mutate(m4_log = log10(`Microglia IV`))
m3m4 <- m3m4 %>% mutate(log_diff = abs(m3_log - m4_log))
m3m4 <- m3m4 %>% filter(log_diff > 1)
m3m4 <- m3m4 %>% filter(m3_log > 0 | m4_log > 0)


## DYNO TRAJECTORY ANALYSIS

## Set seed
set.seed(12345)

## Create dyno object from existing Seurat object
object_counts <- Matrix::t(as(as.matrix(brain.gran@assays$RNA@counts), 'sparseMatrix'))
object_expression <- Matrix::t(as(as.matrix(brain.gran@assays$RNA@data), 'sparseMatrix'))
brain.gran_dyn <- wrap_expression(
  counts = object_counts, 
  expression = object_expression
)

rm(object_counts, object_expression)

## Add a dimensionality reduction
brain.gran_dimred <- dyndimred::dimred_umap(brain.gran_dyn$expression)

## Infer the trajectory
brain.gran_model <- infer_trajectory(brain.gran_dyn, ti_slingshot())

## Plot trajectory & pseudotime
brain.gran_milestone <- plot_dimred(brain.gran_model, label_milestones = T, dimred = brain.gran_dimred) + theme(legend.position = 'none')
brain.gran_traj <- plot_dimred(brain.gran_model, dimred = brain.gran_dimred, grouping = brain.gran@active.ident, color_density = 'grouping') + theme(legend.position = 'none')
brain.gran_pseudo <- plot_dimred(brain.gran_model, "pseudotime", pseudotime = calculate_pseudotime(brain.gran_model), dimred = brain.gran_dimred)

## Check trajectory and root if necessary
brain.gran_model <- add_root(brain.gran_model, root_milestone_id = "4")












