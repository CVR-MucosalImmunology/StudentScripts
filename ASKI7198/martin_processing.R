library(Seurat)
library(GEOquery)
library(Matrix)
library(stringr)

setwd("/Users/thomasoneil/Library/CloudStorage/OneDrive-SharedLibraries-TheUniversityofSydney(Staff)/HIV Pathogenesis & Mucosal Immunology - Documents/Kirstie/Anja/martin/")
source("~/Desktop/.functions.r")

geo = GEOquery::getGEO("GSE134809")
pdat = pData(geo[[1]])

data = list()

i=1
if(pdat$`tissue:ch1`[i] == "Blood") {
  pid = str_split(pdat$title[i], ";")[[1]][1]
} else {
  pid = str_split(pdat$title[i], " ")[[1]][3]
}

mat = readMM(paste0("GSE134809_RAW/",pdat$geo_accession[i],"_",pid, "_matrix.mtx.gz"))
barc = read.csv(paste0("GSE134809_RAW/",pdat$geo_accession[i],"_",pid, "_barcodes.tsv.gz"), header=F)
feat = read.csv(paste0("GSE134809_RAW/",pdat$geo_accession[i],"_",pid, "_genes.tsv.gz"), header=F)
feat$V1 <- sapply(str_split(feat$V1, "\t"), function(x) x[2])
feat$V1 <- make.names(feat$V1, unique = T)

rownames(mat) <- feat$V1
colnames(mat) <- barc$V1

meta = data.frame(row.names = barc$V1, gsm = rep(pdat$geo_accession[i], nrow(barc)), pid = rep(pid, nrow(barc)), tissue = rep(pdat$`tissue:ch1`[i], nrow(barc)), status = rep(pdat$`status:ch1`[i], nrow(barc)), dataset = rep("Martin", nrow(barc)))
temp = CreateSeuratObject(mat, meta.data = meta)

data[[i]] = temp

for(i in 2:nrow(pdat)){
  if(pdat$`tissue:ch1`[i] == "Blood") {
    pid = str_split(pdat$title[i], ";")[[1]][1]
  } else {
    pid = str_split(pdat$title[i], " ")[[1]][3]
  }
  
  mat = readMM(paste0("GSE134809_RAW/",pdat$geo_accession[i],"_",pid, "_matrix.mtx.gz"))
  barc = read.csv(paste0("GSE134809_RAW/",pdat$geo_accession[i],"_",pid, "_barcodes.tsv.gz"), header=F)
  feat = read.csv(paste0("GSE134809_RAW/",pdat$geo_accession[i],"_",pid, "_genes.tsv.gz"), header=F)
  feat$V1 <- sapply(str_split(feat$V1, "\t"), function(x) x[2])
  feat$V1 <- make.names(feat$V1, unique = T)
  
  rownames(mat) <- feat$V1
  colnames(mat) <- barc$V1
  
  meta = data.frame(row.names = barc$V1, gsm = rep(pdat$geo_accession[i], nrow(barc)), pid = rep(pid, nrow(barc)), tissue = rep(pdat$`tissue:ch1`[i], nrow(barc)), status = rep(pdat$`status:ch1`[i], nrow(barc)), dataset = rep("Martin", nrow(barc)))
  temp = CreateSeuratObject(mat, meta.data = meta)
  temp$percent.mt <- PercentageFeatureSet(temp, pattern = "^MT[^a-zA-Z0-9]")
  table(temp$percent.mt >25)
  table(temp$nCount_RNA > 800)
  
  data[[i]] = temp
}
rm(i, pid, temp, meta, mat,barc, feat)
for (i in 1:nrow(pdat)) {
  temp <- data[[i]]
  temp$percent.mt <- PercentageFeatureSet(temp, pattern = "^MT[^a-zA-Z0-9]")
  temp <- subset(temp, subset = percent.mt<20 & nCount_RNA >800)
  data[[i]] = temp
}

for(i in 1:nrow(pdat)) {
  print(dim(data[[i]]))
}

names(data) <- pdat$geo_accession

pbmc <- data[pdat$`tissue:ch1`=='Blood']
tissue <- data[pdat$`tissue:ch1`!='Blood']

tissue <- lapply(X = tissue, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = tissue)

tissue <- lapply(X = tissue, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
  
})

tissue <- FindIntegrationAnchors(object.list = tissue, dims = 1:30, reduction='rpca')
gc()

tissue <- IntegrateData(anchorset = tissue, dims = 1:30)

saveRDS(data, "filtered_preint.rds")
saveRDS(pbmc, "filtered_preint_pbmc.rds")
saveRDS(tissue, "filtered_preint_tissue.rds")


tissue <- ScaleData(tissue, verbose = FALSE)
tissue <- RunPCA(tissue, npcs = 30, verbose = FALSE)
tissue <- RunUMAP(tissue, reduction = "pca", dims = 1:30)
tissue <- RunTSNE(tissue, reduction = "pca", dims = 1:30)
tissue <- FindNeighbors(tissue, dims=1:30)
tissue <- FindClusters(tissue, resolution = 0.05)

DefaultAssay(tissue) <- 'RNA'
tissue <- JoinLayers(tissue)
tissue <- NormalizeData(tissue)

DotPlot(tissue, features=c("PTPRC", "CD3E","CD8A", "FOXP3","NKG7", "KLRB1", "IL7R","JCHAIN", "MS4A1", "CD79A", "ITGAX", "CLEC10A","CD1C", "CLEC9A","XCR1","IL3RA", "HPGDS", "EPCAM", "COL1A1","LYVE1" ), cluster.idents = T)+RotatedAxis()

tissue <- RenameIdents(tissue, 
                       "0" = "T cells",
                       "1" = "Plasma",
                       "2" = "B cells",
                       "4" = "MNP",
                       "7" = "pDC",
                       
                       "6" = "ILCs",
                       "8" = "Mast",
                       
                       "3" = "Epi",
                       "5" = "Fb",
                       "9" = "Fb",
                       "10" = "Endothelial"
                       
                       )
DotPlot(tissue, features=c("PTPRC", "CD3E","CD8A", "FOXP3","NKG7", "JCHAIN", "MS4A1", "CD79A", "ITGAX", "CLEC10A","CD1C", "CLEC9A","XCR1","IL3RA", "KLRB1", "IL7R","HPGDS", "EPCAM", "COL1A1","LYVE1" ), cluster.idents = F)+RotatedAxis()
UMAPPlot(tissue, label=T, label.box=T)
tissue$broad = Idents(tissue)
saveRDS(tissue, "tissue_clusteredAndAnnotated_all.rds")


# recluster MNP
tissue <- subset(tissue, idents = c("MNP", 'pDC'))
tissue <- subset(tissue, subset = pid != "122" & pid!= "123")

tissue <- SplitObject(tissue, split.by = 'pid')

tissue <- lapply(X = tissue, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = tissue)

tissue <- lapply(X = tissue, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
  
})

tissue <- FindIntegrationAnchors(object.list = tissue, dims = 1:30, reduction='rpca')
gc()

tissue <- IntegrateData(anchorset = tissue, dims = 1:30, k.weight=70)

tissue <- ScaleData(tissue, verbose = FALSE)
tissue <- RunPCA(tissue, npcs = 30, verbose = FALSE)
tissue <- RunUMAP(tissue, reduction = "pca", dims = 1:30)
tissue <- RunTSNE(tissue, reduction = "pca", dims = 1:30)
tissue <- FindNeighbors(tissue, dims=1:30)
tissue <- FindClusters(tissue, resolution = 0.5)

DefaultAssay(tissue) <- 'RNA'
tissue <- JoinLayers(tissue)
tissue <- NormalizeData(tissue)

proportions(tissue, "ident", "status", "fill")

UMAPPlot(tissue,label=T, label.box=T)+DotPlot(tissue, features=c("XCR1","CLEC9A", "IL3RA", "AXL", "SIGLEC6", "MKI67", "GZMB"))+RotatedAxis() 
UMAPPlot(tissue,label=T, label.box=T)+DotPlot(tissue, features=c("F13A1", "CLEC10A", "CD207", "CD80"))+RotatedAxis() 
UMAPPlot(tissue,label=T, label.box=T)+DotPlot(tissue, features=c("F13A1", "CLEC10A", "CD207", "CD80", "JCHAIN"))+RotatedAxis()
tissue <- subset(tissue, idents = c(0,1,2,3,4,5,7,8,9,10,11,12)) #minus B and T cells

DefaultAssay(tissue) <- 'integrated'
tissue <- RunPCA(tissue, npcs = 30, verbose = FALSE)
tissue <- RunUMAP(tissue, reduction = "pca", dims = 1:30)
tissue <- FindNeighbors(tissue, dims=1:30)
tissue <- FindClusters(tissue, resolution = 1)

DefaultAssay(tissue) <- 'RNA'
UMAPPlot(tissue,label=T, label.box=T)+DotPlot(tissue, features=c("F13A1", "CLEC10A", "CD207", "CD80"))+RotatedAxis() 
UMAPPlot(tissue,label=T, label.box=T)+DotPlot(tissue, features=c("CLEC9A", "XCR1", "IL3RA"))+RotatedAxis()

proportions(tissue, "ident", "status", "fill")

marks = FindAllMarkers(tissue, only.pos=T)
marks <- marks[marks$p_val_adj<0.0001,]
topm <- marks %>%
  group_by(cluster) %>%
  top_n(n=50, wt=avg_log2FC)

UMAPPlot(tissue,label=T, label.box=T)+DotPlot(tissue, features=c("CD14", "F13A1", "LYVE1", "HPGDS"), idents=c(1,2,3,5,6,7,8,9,10,11,14))+RotatedAxis()

tissue <- RenameIdents(tissue, 
                       "0"="pDC",
                       "1"="Lang+ cDC2.1",#CD207+
                       "2"="Mf.1",
                       "3"="Mf.2",
                       "4"="pDC2",
                       "5"="Monocyte.1",
                       "6"="Mf3",
                       "7"="Monocyte.2",
                       "8"="Mig.DC",
                       "9"="Lang- cDC2",
                       "10"="Lang+ cDC2.2", #CD207+
                       "11"="Act.Mig.DC",
                       "12"="cDC1",
                       "13"="Cycling", 
                       "14" = "Doublet", 
                       "15"="pDC3")
saveRDS(tissue, "mnp_annotated.rds")

