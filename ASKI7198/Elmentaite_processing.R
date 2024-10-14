
setwd("~/Anja/gutcellatlas/c colon immune atlas")
library(Seurat)
library(SeuratDisk)
library(ggplot2)

#  convert (once) and load the data
Convert("raw/Colon_cell_atlas.h5ad", "raw/raw.h5seurat")
data <- LoadH5Seurat("raw/raw.h5seurat")

# add metadata provided by the GCA
data <- AddMetaData(data, read.csv("raw/Colon_immune_metadata.csv"))

## QC checking ####

# cell filtering
png("plots/ncount_ngene.png", width=360, height=360);
FeatureScatter(data, 'n_genes', 'n_counts')+
  labs(subtitle = "min counts: 1114; min-max genes: 700-6463")+NoLegend()
dev.off()

min(data$n_counts);max(data$n_counts)
min(data$n_genes);max(data$n_genes)

#percen.mt
png("plots/percentmt.png", width=360, height=360);
VlnPlot(data, 'percent_mito', pt.size=0)+
  labs(subtitle = "Cutoff at 10% (removes 4% of data)")+NoLegend()+xlab("")+geom_hline(yintercept = 0.1)
dev.off()

table(data$percent_mito<0.1)
754/nrow(data)

## meta tables ####

meta <- data@meta.data
write.csv(meta,"raw/meta.csv")

table(data$donor)

#$Age
png("plots/MetaGender.png", width=360, height=200);
ggplot(meta, aes(gender, fill=gender)) +
  geom_bar(stat = 'count')+theme_classic()+RotatedAxis()+
  labs(title = "Gender")
dev.off()

#$region
png("plots/MetaRegion.png", width=360, height=200);
ggplot(meta, aes(region, fill=region)) +
  geom_bar(stat = 'count')+theme_classic()+RotatedAxis()+
  labs(title = "Region")
dev.off()

#$cell_type
png("plots/MetaCellType.png", width=360, height=200);
ggplot(meta, aes(cell_type, fill=Fraction)) +
  geom_bar(stat = 'count')+theme_classic()+RotatedAxis()+
  labs(title = "Annotations by the GCA")
dev.off()


### processing ####

#* 1 is all
#* 2 is myeloid
#* 3 is T cells
#* 4 is B cells

data <- subset(data, subset = percent_mito < 0.1)
data <- NormalizeData(data)

# 1. all

temp <- SplitObject(data, split.by='donor')
# temp <- temp[sapply(temp, function(x) ncol(x) > 80)]

temp <- lapply(X = temp, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = temp)
temp <- lapply(X = temp, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
temp <- FindIntegrationAnchors(object.list = temp, dims = 1:30, reduction='rpca')
temp <- IntegrateData(anchorset = temp, k.weight=100, #lower if low cell per donor 
                      dims = 1:30)

temp <- ScaleData(temp)
temp <- RunPCA(temp)
temp <- RunUMAP(temp, dims = 1:30)
UMAPPlot(temp, group.by='cell_type', label=T)
Idents(temp) <- "cell_type"
temp <- RunTSNE(temp, dims=1:30)
TSNEPlot(temp, label=T)

temp <- RenameIdents(temp, 
                     "CD8 T" = "T cell", 
                     "Tcm" = "T cell", 
                     "Treg" = "T cell", 
                     "Tfh" = "T cell", 
                     "Activated CD4 T" = "T cell", 
                     "Th1" = "T cell", 
                     "Th17" = "T cell", 
                     'gd T' = "T cell",
                     "NK" = "NK",
                     "ILC" = "ILC", 
                     "Mast" = "Mast", 
                     "B cell IgG Plasma" = "Plasma", 
                     "B cell IgA Plasma" = "Plasma", 
                     "cycling gd T" = "Cycling", 
                     "B cell cycling" = "Cycling", 
                     "B cell memory" = "B cell", 
                     "Follicular B cell" = "B cell", 
                     "Monocyte" = "MNP", 
                     "LYVE1 Macrophage" = "MNP", 
                     "Macrophage" = "MNP", 
                     "cycling DCs" = "MNP", 
                     "cDC1" = "MNP", 
                     "cDC2" = "MNP", 
                     "Lymphoid DC" = "MNP", 
                     "pDC" = "MNP"
)

TSNEPlot(temp, label=T, pt.size=2)

png("plots/All_TSNE.png", width=360, height=200);
TSNEPlot(temp, label=T, pt.size=2)
dev.off()
png("plots/All_UMAP.png", width=360, height=200);
UMAPPlot(temp, label=T, pt.size=2)
dev.off()
png("plots/All_UMAP_region.png", width=1200, height=350);
UMAPPlot(temp, label=F, pt.size=1.5, group.by="cell_type", split.by="region")
UMAPPlot(temp, label=F, pt.size=1.5, split.by="region")
dev.off()

DefaultAssay(temp) <- 'RNA'
temp <- ScaleData(temp)

saveRDS(temp, "processed/a_alldata.rds")
gc()
data <-temp

#* 2 is myeloid
temp <- subset(data, ident="MNP")

table(temp$donor)

temp <- SplitObject(temp, split.by='donor')

temp <- lapply(X = temp, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = temp)
temp <- lapply(X = temp, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE, npcs=30)
})
temp <- FindIntegrationAnchors(object.list = temp, dims = 1:30, reduction='rpca')
temp <- IntegrateData(anchorset = temp, k.weight=30, #lower if low cell per donor 
                      dims = 1:30)

temp <- ScaleData(temp)
temp <- RunPCA(temp)
temp <- RunUMAP(temp, dims = 1:30)
UMAPPlot(temp, group.by='cell_type', label=T)
Idents(temp) <- "cell_type"
temp <- RunTSNE(temp, dims=1:30)
TSNEPlot(temp, label=T)
TSNEPlot(temp, label=T, group.by="region")

DefaultAssay(temp) <- 'RNA'

png("plots/myeloid_TSNE.png", width=360, height=200);
TSNEPlot(temp, label=T, pt.size=2)
FeaturePlot(temp,"CD207", reduction='tsne')
dev.off()
png("plots/myeloid_UMAP.png", width=360, height=200);
UMAPPlot(temp, label=T, pt.size=2)
FeaturePlot(temp,"CD207", reduction='umap')
dev.off()

temp <- ScaleData(temp)

saveRDS(temp, "processed/b_myeloid.rds")
gc()


