setwd("C:/Users/erica/OneDrive - The University of Sydney (Staff)/PhD - HIV/RNA Analysis/Bujko")
source("https://raw.githubusercontent.com/DrThomasOneil/CVR-site/refs/heads/master/data/.functions.R")

library(Seurat)
library(ggplot2)
library(tidyverse)
library(dplyr)

# original file from TO

data <- readRDS("Mo_int.rds")

# my ongoing analysis

data <- readRDS("ericaanalysis_241217.rds")

# REMEMBER TO LOAD ENVIRONMENT WITH RDS AND R FILES

# to find clusters I increased resolution until I hit a max of clusters
# then I worked back down, choosing the last one to have the small populations
# present while still being identifiable

data <- FindClusters(data, resolution=0.6)

UMAPPlot(data, pt.size=2, label=T)

table(data$Init_Anno, Idents(data))

?Seurat::FindSubCluster
data <- FindSubCluster(data, "??", graph.name = "RNA_snn", resolution = 0.4, subcluster.name = "new_clust")
UMAPPlot(data, pt.size=2, label=T, group.by="new_clust")
data$old_clust <- Idents(data)
Idents(data) <- "new_clust"
UMAPPlot(data, pt.size=2, label=T)

?Seurat::FindSubCluster
data <- FindSubCluster(data, "Mf2>Mf3", graph.name = "RNA_snn", resolution = 0.1, subcluster.name = "new_clust2")
UMAPPlot(data, pt.size=2, label=T, group.by="new_clust2")
data$old_clust <- Idents(data)
Idents(data) <- "new_clust2"
UMAPPlot(data, pt.size=2, label=T)

?Seurat::FindSubCluster
data <- FindSubCluster(data, "Mf1.2", graph.name = "RNA_snn", resolution = 0.2, subcluster.name = "new_clust3")
UMAPPlot(data, pt.size=2, label=T, group.by="new_clust3")
data$old_clust <- Idents(data)
Idents(data) <- "new_clust3"
UMAPPlot(data, pt.size=2, label=T)

?Seurat::FindSubCluster
data <- FindSubCluster(data, "?Mono.1", graph.name = "RNA_snn", resolution = 0.2, subcluster.name = "new_clust4")
UMAPPlot(data, pt.size=2, label=T, group.by="new_clust4")
data$old_clust <- Idents(data)
Idents(data) <- "new_clust4"
UMAPPlot(data, pt.size=2, label=T)

data <- FindSubCluster(data, "?Mono.1", graph.name = "RNA_snn", resolution = 0.4, subcluster.name = "new_clust5")
UMAPPlot(data, pt.size=2, label=T, group.by="new_clust5")
data$old_clust <- Idents(data)
Idents(data) <- "new_clust5"
UMAPPlot(data, pt.size=2, label=T)

data <- FindSubCluster(data, "9", graph.name = "RNA_snn", resolution = 0.4, subcluster.name = "new_clust6")
UMAPPlot(data, pt.size=2, label=T, group.by="new_clust6")
data$old_clust <- Idents(data)
Idents(data) <- "new_clust6"
UMAPPlot(data, pt.size=2, label=T)

data <- FindSubCluster(data, "7", graph.name = "RNA_snn", resolution = 0.4, subcluster.name = "new_clust7")
UMAPPlot(data, pt.size=2, label=T, group.by="new_clust7")
data$old_clust <- Idents(data)
Idents(data) <- "new_clust7"
UMAPPlot(data, pt.size=2, label=T)

data <- FindSubCluster(data, "1", graph.name = "RNA_snn", resolution = 0.4, subcluster.name = "new_clust8")
UMAPPlot(data, pt.size=2, label=T, group.by="new_clust8")
data$old_clust <- Idents(data)
Idents(data) <- "new_clust8"
UMAPPlot(data, pt.size=2, label=T)


data <- FindSubCluster(data, "5", graph.name = "RNA_snn", resolution = 0.4, subcluster.name = "new_clust9")
UMAPPlot(data, pt.size=2, label=T, group.by="new_clust9")
data$old_clust <- Idents(data)
Idents(data) <- "new_clust9"
UMAPPlot(data, pt.size=2, label=T)


data <- FindSubCluster(data, "0", graph.name = "RNA_snn", resolution = 0.4, subcluster.name = "new_clust10")
UMAPPlot(data, pt.size=2, label=T, group.by="new_clust10")
data$old_clust <- Idents(data)
Idents(data) <- "new_clust10"
UMAPPlot(data, pt.size=2, label=T)


data <- FindSubCluster(data, "2", graph.name = "RNA_snn", resolution = 0.4, subcluster.name = "new_clust11")
UMAPPlot(data, pt.size=2, label=T, group.by="new_clust11")
data$old_clust <- Idents(data)
Idents(data) <- "new_clust11"
UMAPPlot(data, pt.size=2, label=T)


PCAPlot(data, pt.size=2, label=T)

# dotplot started as all the genes I could match to the INF MNP Flow Panel
# I then annotated the clusters based on that and confirmed it with differential genes
# from the Domanska et al. paper

DotPlot(data, 
        features= c("S100A12", "PROK2", "CLDN1", "POLN", "CCL13", "CLEC10A", "ACP5", "CD9", "CD81", "CXCL3", "IL1B", "FCGR3A", "S100A8", "S100A9", "C5AR1", "MKI67", "FCER1A", "CD1C", "CCR7", "ITGAX", "MRC1", "SAMHD1", "CD14", "HLA-DRB1", "CD163", "ITGAM", "F13A1", "LYVE1", "AXL", "SIGLEC1", "CD209", "CD4"), 
        cols=c("yellow", 'blue'))+ coord_flip()#+ RotatedAxis(45)


annotations = c(
  "10" = "Cycling DC",
  "5" = "DC2/3",
  "8" = "??",
  "3" = "?Mono.1",
  "6" = "?Mono.2",
  "13" = "CD16 Mono",
  "7" = "Mf1.1",
  "0" = "Mf1.2",
  "11" = "Mf2.1",
  "1" = "Mf1>Mf2",
  "2" = "Mf2.2",
  "4" = "Mf2>Mf3",
  "9" = "Mf3",
  "12" = "Nothing")

data <- RenameIdents(data, annotations)
data$Init_Anno <- Idents(data)


annotations = c(
  "Cycling DC" = "Cycling DC",
  "DC2/3" = "DC2/3",
  "Mf1>Mf2" = "MDDC",
  "CD16 Mono" = "CD16 Mono",
  "?Mono.1" = "Monocytes",
  "?Mono.2" = "Monocytes",
  "Mf1.1" = "MDM1",
  "Mf1.2" = "MDM2",
  "Mf2.1" = "MDM3",
  "??" = "Unknown",
  "Mf2.2" = "ED-Resident1",
  "Mf2>Mf3" = "ED-Resident2",
  "Mf3" = "M2 Resident",
  "Nothing" = "Junk")

data <- RenameIdents(data, annotations)
data$Init_Anno <- Idents(data)


annotations = c(
  "Cycling DC" = "Cycling DC",
  "DC2/3" = "DC2/3",
  "Mf1>Mf2" = "MDDC",
  "?Mono.2" = "Monocytes1",
  "?Mono.1" = "Monocytes2",
  "CD16 Mono" = "Blood Monocytes",
  "Mf1.1" = "MDM1",
  "Mf1.2" = "MDM2",
  "Mf2.1" = "MDM3",
  "??" = "Unknown",
  "Mf2.2" = "ED-Resident",
  "Mf2>Mf3" = "Resident",
  "Mf3" = "M2 Resident",
  "Nothing" = "Junk")

data <- RenameIdents(data, annotations)
data$Init_Anno <- Idents(data)

#mutate(data, ident = Init_Anno)

#data <- FindClusters(data, resolution=3)

UMAPPlot(data, pt.size=2, label=T)


DotPlot(data, 
        features= c("FCGR3A", "S100A8", "S100A9", "C5AR1", "MKI67", "FCER1A", "CD1C", "CCR7", "ITGAX", "MRC1", "CD14", "HLA-DRB1", "CD163", "ITGAM", "F13A1", "LYVE1", "AXL", "SIGLEC1", "CD209", "CD4"), 
        cols=c("yellow", 'blue'))+ coord_flip()#+ RotatedAxis(45)


DotPlot(data, 
        features= c("CD63", "ADAMDEC1", "DNASE1L3", "ITGAE", "CD274", "CD80", "CD86", "DHRS9", "CLEC10A", "ACP5", "CD9", "CD81", "CXCL3", "IL1B", "FCGR3A", "S100A8", "S100A9", "C5AR1", "MKI67", "FCER1A", "CD1C", "CCR7", "ITGAX", "MRC1", "SAMHD1", "CD14", "HLA-DRB1", "CD163", "ITGAM", "F13A1", "LYVE1", "AXL", "SIGLEC1", "CD209", "CD4"), 
        cols=c("yellow", 'blue'))+ coord_flip()#+ RotatedAxis(45)

?DotPlot


# ChatGPT: how to order clusters in a dotplot by the expression of one gene
# DOESN'T WORK!! - attempt to re-write now I kinda get it
# Assume 'gene_name' is your gene of interest
gene_name <- "HLA-DRB1"

#View(data)

mutate()

# Fetch and calculate average expression
# DOESN"T WORK YET:
average_expression <- FetchData(data, vars = gene_name) %>%
  data.frame() %>%
  cbind(idents = data@meta.data$Init_Anno) %>%
  group_by(Init_Anno) #>%
  summarise(average_expression = mean(!!sym(gene_name), na.rm = TRUE)) #%>%
  arrange(desc(average_expression))
  
?group_by
?cbind
?summarise
?ddply
  
# Reorder idents based on average expression
ordered_idents <- average_expression$idents

# Set new order in Seurat object
data <- SetIdent(data, value = factor(data$ident, levels = ordered_idents))

# Create the dot plot
DotPlot(data, features = gene_name) + 
  RotatedAxis()






# add below to above to include UMAP on same plot
# +UMAPPlot(data, pt.size=2, label=T)

# I checked all the feature plots for the genes of interest from the dotplot
# only kept those that had partitioned spatial expression

FeaturePlot(data, 
            c("ITGAM"), 
            pt.size=2, 
            label=F, order=T, min.cutoff = 1)

FeaturePlot(data, 
            c("FCGR3A",
              "S100A8",
              "F13A1",
              "CD1C",
              "ITGAM",
              "MKI67",
              "CCR7",
              "FCER1A",
              "AXL",
              "MRC1",
              "SIGLEC1",
              "CD209"), 
            pt.size=2, 
            label=F, order=T, min.cutoff = 1)

FeaturePlot(data, 
            c("CD9",
              "APC5",
              "DNASE1L3",
              "ADAMDEC1",
              "CD63",
              "ITGAX",
              "MRC1"), 
            pt.size=2, 
            label=F, order=T, min.cutoff = 1)



View(data@meta.data)


# not sure that the first two are for, but TO notes had it, so used

remotes::install_github('satijalab/seurat-wrappers')

data <- ScaleData(data)

# this version will take all clusters and compare to all other clusters. It'll pull out the most differential features for each subset. It can also take some time. 
all_markers <- FindAllMarkers(data, logfc.threshold = 0.5, only.pos= T) # we use only.pos to get the enriched markers, not the markers that are decreased on each subset.

# filter

all_markers_sig <- all_markers[all_markers$p_val_adj <0.001,]

top10 <- all_markers_sig %>% #this is a pipe
  group_by(cluster) %>% #by each cluster
  top_n(n = 10, wt = avg_log2FC)#take the top n (10) of the avg_log2FC 

#top10 is now a list of only the top 10 genes per cluster - 13 clusters = 130 genes
# and now we can do: 

#TO: to make sure there was more than the 2000 genes in the scale.data
data <- ScaleData(data, features=rownames(data))

DoHeatmap(data, features=top10$gene, size = 3) #or
DoHeatmap(AverageExpression(data, return.seurat = T), features=top10$gene, draw.lines=F)+NoLegend()

?DoHeatmap

top5 <- all_markers_sig %>% #this is a pipe
  group_by(cluster) %>% #by each cluster
  top_n(n = 5, wt = avg_log2FC)#take the top n (10) of the avg_log2FC 

#top10 is now a list of only the top 10 genes per cluster - 13 clusters = 130 genes
# and now we can do: 

#TO: to make sure there was more than the 2000 genes in the scale.data
data <- ScaleData(data, features=rownames(data))

DoHeatmap(data, features=top5$gene, size = 2) #or
DoHeatmap(AverageExpression(data, return.seurat = T), features=top5$gene, size = 3, draw.lines=F)


top3 <- all_markers_sig %>% #this is a pipe
  group_by(cluster) %>% #by each cluster
  top_n(n = 3, wt = avg_log2FC)#take the top n (10) of the avg_log2FC 

#top10 is now a list of only the top 10 genes per cluster - 13 clusters = 130 genes
# and now we can do: 

#TO: to make sure there was more than the 2000 genes in the scale.data
data <- ScaleData(data, features=rownames(data))

DoHeatmap(data, features=top3$gene, size = 3) #or
DoHeatmap(AverageExpression(data, return.seurat = T), features=top3$gene, size = 3, draw.lines=F)


#proportions of donor in each cluster

proportions <- data %>%
  group_by(seurat_cluster, donor) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

?

  
proportions <- function(data, ident.1, ident.2, position) {
  x<- FetchData(data,c(ident.1,ident.2))
  colnames(x) <- c('ident.2', 'ident.1')
  x%>% 
    group_by(ident.1) %>%
    mutate(prop=1/length(ident.2)) %>%
    ungroup() %>%
    group_by(ident.2,ident.1) %>%
    summarise(totprop=sum(prop)) %>%
    
    ggplot(aes(x=ident.2,fill=ident.1,y=totprop)) +
    geom_bar(position=position, stat='identity') + 
    theme(axis.text.x =element_text(angle = 45,hjust=1))+
    scale_y_continuous(name="Cluster Proportion")+ 
    theme_classic()
}
  
  
  
?proportions


# comparing the top10 differential genes for all the cluster that include Mf2
#logfc threshold - looking for bigger changes, delete for smaller changes or if you don't know

AllMf2_markers <- FindMarkers(data, logfc.threshold = 0.5, ident.1 = 'Mf2.1', ident.2='Mf2.2', ident.3 = 'Mf2>Mf3', ident.4 = 'Mf1>Mf2') 

# TO: I like to add a column here - if the pct is marginal, you might want to know that. 
# TO: e.g. 99% expressed vs 98% expressed 
# TO: however, amount of expression is just as important. 
AllMf2_markers$pct.diff <- AllMf2_markers$pct.1- AllMf2_markers$pct.2

AllMf2_markers$up <- ifelse(AllMf2_markers$avg_log2FC<0, "down", "up")
AllMf2_markers$gene <- rownames(AllMf2_markers)

top10 <- rbind(AllMf2_markers %>% filter(up == 'up') %>% top_n(n=20, wt=avg_log2FC),
               AllMf2_markers %>% filter(up == 'down') %>% top_n(n=-20, wt=avg_log2FC))

top10 <- rbind(AllMf2_markers %>% filter(up == 'up') %>% top_n(n=20, wt=avg_log2FC),
               AllMf2_markers %>% filter(up == 'down') %>% top_n(n=-20, wt=avg_log2FC))

# comparing the top10 differential genes for all the cluster that include Mf2

Mf3_markers <- FindMarkers(data, ident.1 = 'Mf3', ident.2='Mf2>Mf3') 
Mf3_markers[grep("ITGAM", rownames(Mf3_markers)),]
# TO: I like to add a column here - if the pct is marginal, you might want to know that. 
# TO: e.g. 99% expressed vs 98% expressed 
# TO: however, amount of expression is just as important. 
Mf3_markers$pct.diff <- Mf3_markers$pct.1- Mf3_markers$pct.2

Mf3_markers$up <- ifelse(Mf3_markers$avg_log2FC<0, "down", "up")
Mf3_markers$gene <- rownames(Mf3_markers)

top10 <- rbind(Mf3_markers %>% filter(up == 'up') %>% top_n(n=20, wt=avg_log2FC),
               Mf3_markers %>% filter(up == 'down') %>% top_n(n=-20, wt=avg_log2FC))



AllMono_markers <- FindMarkers(data, logfc.threshold = 0.5, ident.1 ='?Mono.1', ident.2 = '?Mono.2') 

# TO: I like to add a column here - if the pct is marginal, you might want to know that. 
# TO: e.g. 99% expressed vs 98% expressed 
# TO: however, amount of expression is just as important. 
AllMono_markers$pct.diff <- AllMono_markers$pct.1- AllMono_markers$pct.2

AllMono_markers$up <- ifelse(AllMono_markers$avg_log2FC<0, "down", "up")
AllMono_markers$gene <- rownames(AllMono_markers)

top10 <- rbind(AllMono_markers %>% filter(up == 'up') %>% top_n(n=20, wt=avg_log2FC),
               AllMono_markers %>% filter(up == 'down') %>% top_n(n=-20, wt=avg_log2FC))

DoHeatmap(data, features=top10$gene, size = 3) 


AllDC_markers <- FindMarkers(data, logfc.threshold = 0.5, ident.1 = 'Cycling DC', ident.2='DC2/3') 

# TO: I like to add a column here - if the pct is marginal, you might want to know that. 
# TO: e.g. 99% expressed vs 98% expressed 
# TO: however, amount of expression is just as important. 
AllDC_markers$pct.diff <- AllDC_markers$pct.1- AllDC_markers$pct.2

AllDC_markers$up <- ifelse(AllDC_markers$avg_log2FC<0, "down", "up")
AllDC_markers$gene <- rownames(AllDC_markers)

top10 <- rbind(AllDC_markers %>% filter(up == 'up') %>% top_n(n=20, wt=avg_log2FC),
               AllDC_markers %>% filter(up == 'down') %>% top_n(n=-20, wt=avg_log2FC))

DoHeatmap(data, features=top10$gene, size = 3) 


AllDCMf2_markers <- FindMarkers(data, logfc.threshold = 0.5, ident.1 = 'Mf1>Mf2', ident.2='DC2/3', ident.3='Mf2.1', ident.4='Mf1.2') 

# TO: I like to add a column here - if the pct is marginal, you might want to know that. 
# TO: e.g. 99% expressed vs 98% expressed 
# TO: however, amount of expression is just as important. 
AllDCMf2_markers$pct.diff <- AllDCMf2_markers$pct.1- AllDCMf2_markers$pct.2

AllDCMf2_markers$up <- ifelse(AllDCMf2_markers$avg_log2FC<0, "down", "up")
AllDCMf2_markers$gene <- rownames(AllDCMf2_markers)

top10 <- rbind(AllDCMf2_markers %>% filter(up == 'up') %>% top_n(n=20, wt=avg_log2FC),
               AllDCMf2_markers %>% filter(up == 'down') %>% top_n(n=-20, wt=avg_log2FC))

DoHeatmap(data, features=top10$gene, size = 3) 


AllMf2Mf1_markers <- FindMarkers(data, logfc.threshold = 0.5, ident.1 = 'Mf1>Mf2', ident.2='Mf2.2', ident.3='Mf1.1', ident.4='Mf1.2') 

# TO: I like to add a column here - if the pct is marginal, you might want to know that. 
# TO: e.g. 99% expressed vs 98% expressed 
# TO: however, amount of expression is just as important. 
AllMf2Mf1_markers$pct.diff <- AllMf2Mf1_markers$pct.1- AllMf2Mf1_markers$pct.2

AllMf2Mf1_markers$up <- ifelse(AllMf2Mf1_markers$avg_log2FC<0, "down", "up")
AllMf2Mf1_markers$gene <- rownames(AllMf2Mf1_markers)

top10 <- rbind(AllMf2Mf1_markers %>% filter(up == 'up') %>% top_n(n=20, wt=avg_log2FC),
               AllMf2Mf1_markers %>% filter(up == 'down') %>% top_n(n=-20, wt=avg_log2FC))

DoHeatmap(data, features=top10$gene, size = 3) 


AllMf3_markers <- FindMarkers(data, logfc.threshold = 0.5, ident.1 = 'Mf2>Mf3', ident.2='Mf2.2', ident.3='Mf3') 

# TO: I like to add a column here - if the pct is marginal, you might want to know that. 
# TO: e.g. 99% expressed vs 98% expressed 
# TO: however, amount of expression is just as important. 
AllMf3_markers$pct.diff <- AllMf3_markers$pct.1- AllMf3_markers$pct.2

AllMf3_markers$up <- ifelse(AllMf3_markers$avg_log2FC<0, "down", "up")
AllMf3_markers$gene <- rownames(AllMf3_markers)

top10 <- rbind(AllMf3_markers %>% filter(up == 'up') %>% top_n(n=20, wt=avg_log2FC),
               AllMf3_markers %>% filter(up == 'down') %>% top_n(n=-20, wt=avg_log2FC))

DoHeatmap(data, features=top10$gene, size = 3) 


AllMf3_markers <- FindMarkers(data, logfc.threshold = 0.5, ident.1 ='Mf2.2', ident.2='Mf3') 

# TO: I like to add a column here - if the pct is marginal, you might want to know that. 
# TO: e.g. 99% expressed vs 98% expressed 
# TO: however, amount of expression is just as important. 
AllMf3_markers$pct.diff <- AllMf3_markers$pct.1- AllMf3_markers$pct.2

AllMf3_markers$up <- ifelse(AllMf3_markers$avg_log2FC<0, "down", "up")
AllMf3_markers$gene <- rownames(AllMf3_markers)

top10 <- rbind(AllMf3_markers %>% filter(up == 'up') %>% top_n(n=20, wt=avg_log2FC),
               AllMf3_markers %>% filter(up == 'down') %>% top_n(n=-20, wt=avg_log2FC))

DoHeatmap(data, features=top10$gene, size = 3) 


AllMf3_markers <- FindMarkers(data, logfc.threshold = 0.5, ident.1 ='Mf2.2', ident.2='Mf2>Mf3') 

# TO: I like to add a column here - if the pct is marginal, you might want to know that. 
# TO: e.g. 99% expressed vs 98% expressed 
# TO: however, amount of expression is just as important. 
AllMf3_markers$pct.diff <- AllMf3_markers$pct.1- AllMf3_markers$pct.2

AllMf3_markers$up <- ifelse(AllMf3_markers$avg_log2FC<0, "down", "up")
AllMf3_markers$gene <- rownames(AllMf3_markers)

top10 <- rbind(AllMf3_markers %>% filter(up == 'up') %>% top_n(n=20, wt=avg_log2FC),
               AllMf3_markers %>% filter(up == 'down') %>% top_n(n=-20, wt=avg_log2FC))

DoHeatmap(data, features=top10$gene, size = 3) 


AllMf3_markers <- FindMarkers(data, logfc.threshold = 0.5, ident.1 ='Mf3', ident.2='Mf2>Mf3') 

# TO: I like to add a column here - if the pct is marginal, you might want to know that. 
# TO: e.g. 99% expressed vs 98% expressed 
# TO: however, amount of expression is just as important. 
AllMf3_markers$pct.diff <- AllMf3_markers$pct.1- AllMf3_markers$pct.2

AllMf3_markers$up <- ifelse(AllMf3_markers$avg_log2FC<0, "down", "up")
AllMf3_markers$gene <- rownames(AllMf3_markers)

top10 <- rbind(AllMf3_markers %>% filter(up == 'up') %>% top_n(n=20, wt=avg_log2FC),
               AllMf3_markers %>% filter(up == 'down') %>% top_n(n=-20, wt=avg_log2FC))

DoHeatmap(data, features=top10$gene, size = 3) 


AllMf3_markers <- FindMarkers(data, logfc.threshold = 0.5, ident.1 ='Mf3', ident.2='??') 

# TO: I like to add a column here - if the pct is marginal, you might want to know that. 
# TO: e.g. 99% expressed vs 98% expressed 
# TO: however, amount of expression is just as important. 
AllMf3_markers$pct.diff <- AllMf3_markers$pct.1- AllMf3_markers$pct.2

AllMf3_markers$up <- ifelse(AllMf3_markers$avg_log2FC<0, "down", "up")
AllMf3_markers$gene <- rownames(AllMf3_markers)

top10 <- rbind(AllMf3_markers %>% filter(up == 'up') %>% top_n(n=20, wt=avg_log2FC),
               AllMf3_markers %>% filter(up == 'down') %>% top_n(n=-20, wt=avg_log2FC))

DoHeatmap(data, features=top10$gene, size = 3) 


AllMf1_markers <- FindMarkers(data, logfc.threshold = 0.5, ident.1 ='Mf1.1', ident.2 ='Mf1.2', ident.3 ="Mf1>Mf2") 

# TO: I like to add a column here - if the pct is marginal, you might want to know that. 
# TO: e.g. 99% expressed vs 98% expressed 
# TO: however, amount of expression is just as important. 
AllMf1_markers$pct.diff <- AllMf1_markers$pct.1- AllMf1_markers$pct.2

AllMf1_markers$up <- ifelse(AllMf1_markers$avg_log2FC<0, "down", "up")
AllMf1_markers$gene <- rownames(AllMf1_markers)

top10 <- rbind(AllMf1_markers %>% filter(up == 'up') %>% top_n(n=20, wt=avg_log2FC),
               AllMf1_markers %>% filter(up == 'down') %>% top_n(n=-20, wt=avg_log2FC))

DoHeatmap(data, features=top10$gene, size = 3) 


AllMf1_markers <- FindMarkers(data, logfc.threshold = 0.5, ident.1 ='Mf1.2', ident.2 ='Mf1>Mf2') 

# TO: I like to add a column here - if the pct is marginal, you might want to know that. 
# TO: e.g. 99% expressed vs 98% expressed 
# TO: however, amount of expression is just as important. 
AllMf1_markers$pct.diff <- AllMf1_markers$pct.1- AllMf1_markers$pct.2

AllMf1_markers$up <- ifelse(AllMf1_markers$avg_log2FC<0, "down", "up")
AllMf1_markers$gene <- rownames(AllMf1_markers)

top10 <- rbind(AllMf1_markers %>% filter(up == 'up') %>% top_n(n=20, wt=avg_log2FC),
               AllMf1_markers %>% filter(up == 'down') %>% top_n(n=-20, wt=avg_log2FC))

DoHeatmap(data, features=top10$gene, size = 3) 


AllMf1_markers <- FindMarkers(data, logfc.threshold = 0.5, ident.1 ='MDM', ident.2 ='Unknown') 

# TO: I like to add a column here - if the pct is marginal, you might want to know that. 
# TO: e.g. 99% expressed vs 98% expressed 
# TO: however, amount of expression is just as important. 
AllMf1_markers$pct.diff <- AllMf1_markers$pct.1- AllMf1_markers$pct.2

AllMf1_markers$up <- ifelse(AllMf1_markers$avg_log2FC<0, "down", "up")
AllMf1_markers$gene <- rownames(AllMf1_markers)

top10 <- rbind(AllMf1_markers %>% filter(up == 'up') %>% top_n(n=20, wt=avg_log2FC),
               AllMf1_markers %>% filter(up == 'down') %>% top_n(n=-20, wt=avg_log2FC))

DoHeatmap(data, features=top10$gene, size = 3) 


colour_scale <- c("grey", "cyan", "black", "pink")

#data$(meta.data) so I can see what the labels on the data are

table(data$sample)
proportions(data, 'ident', 'sample', 'fill')
?proportions


# heatmap of top10 genes for all annotated clusters
# keep getting error that not all of the genes are listed in the scaled.data RNA assay,
# and at least of half of the top10 genes for all the clusters aren't there

?DoHeatmap

DoHeatmap(subset(data, idents=c("Mf2.1", "Mf2.2", "Mf1>Mf2", "Mf2>Mf3")), 
          features=top10$gene, 
          draw.lines=F,
          angle = 0,
          size = 4,
          group.color = colour_scale)+
          scale_fill_gradientn(colors = c("blue", "white", "red"))


DotPlot(subset(data, idents=c("Mf3_0", "Mf3_1")), 
            features=top10$gene)+coord_flip()


DotPlot(data, 
        features=top10$gene)+coord_flip()



# set colours
install.packages("RColorBrewer")

colour_scale <- brewer.pal(n = 3, name = "Set1")

custom_gradient <- colorRampPalette(c("blue", "white", "red"))(100)

#hex colours
custom_colours <- c("#FF0000", "#00FF00", "#0000FF")



# this method, you can compare two subsets specifically. 
cDC1cDC2_markers <- FindMarkers(data, logfc.threshold = 0.5, ident.1 = 'cDC1', ident.2='cDC2') 
# I like to add a column here - if the pct is marginal, you might want to know that. 
# e.g. 99% expressed vs 98% expressed 
# however, amount of expression is just as important. 
cDC1cDC2_markers$pct.diff <- cDC1cDC2_markers$pct.1- cDC1cDC2_markers$pct.2

cDC1cDC2_markers$up <- ifelse(cDC1cDC2_markers$avg_log2FC<0, "down", "up")
cDC1cDC2_markers$gene <- rownames(cDC1cDC2_markers)

top10 <- rbind(cDC1cDC2_markers %>% filter(up == 'up') %>% top_n(n=20, wt=avg_log2FC),
               cDC1cDC2_markers %>% filter(up == 'down') %>% top_n(n=-20, wt=avg_log2FC))


DoHeatmap(subset(data, idents=c("cDC1", "cDC2")), 
          features=top10$gene, 
          draw.lines=F)+NoLegend()


##### Measuring Expression #####

FeatureScatter(subset(data, idents=c("Mf1.1", "Mf1.2")), "HLA-DRB1","ITGAM",pt.size=2)
FeatureScatter(subset(data, idents=c('Cycling DC', "Mf3_1")), "CCR7","LYVE1",pt.size=2)
FeatureScatter(subset(data, idents=c("Mf3", "Mf1.2")), "ITGAM","CCL13",pt.size=2)


FeatureScatter(subset(data, idents=c("Mf3", "Mf2>Mf3", "Mf2.2", "??")), "CD14","F13A1",pt.size=2)
FeatureScatter(subset(data, idents=c("Mf3", "Mf2>Mf3", "Mf2.2", "??")), "ITGAM","FCGR3A",pt.size=2)
FeatureScatter(subset(data), "HLA-DRB1","S100A8",pt.size=2)
FeatureScatter(subset(data, idents=c("Cycling DC", "DC2/3")), "F13A1","S100A8",pt.size=2)


# here we'll just use a pre-defined set of genes

# tip: use grep to find gene names - e.g. sometimes HLA-DRB1 is stored as HLA.DRB1
grep("^CCR", # some pattern to look for - the ^ at the start specifies that the pattern 'starts with'
     rownames(data), # where to look
     value=T) #return the actual names





# Attempting to add pseudotime analysis to UMAP

install.packages("BiocManager")
BiocManager::install("monocle")

browseVignettes("monocle")

library(monocle)



# instructions at https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                        'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                        'SummarizedExperiment', 'batchelor', 'HDF5Array',
                        'terra', 'ggrastr'))


install.packages(c("BiocGenerics", "DelayedArray", "DelayedMatrixStats", "limma", 
                   "lme4", "S4Vectors", "SingleCellExperiment", "SummarizedExperiment", 
                   "batchelor", "HDF5Array", "ggrastr"), 
                 dependencies = TRUE, 
                 force = TRUE)

install.packages("devtools")

# my_check <- check
# rm(check)

# pkgbuild::check_build_tools(debug = TRUE)

devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle3)


# expression_matrix <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_expression.rds"))
expression_matrix <- data[["RNA"]]$counts

# cell_metadata <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_colData.rds"))
cell_metadata <- data@meta.data

#gene_annotation <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/packer_embryo_rowData.rds"))
gene_annotation <- data.frame(gene_short_name = rownames(data), row.names = rownames(data))

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 50)
#cds <- preprocess_method()

cds <- align_cds(cds, alignment_group = "sample")

cds <- reduce_dimension(cds, preprocess_method = 'Aligned', reduction_method = 'UMAP')
plot_cells(cds, cell_size = 0.4, label_cell_groups=T,  color_cells_by = "Init_Anno", group_label_size = 4)
cds <- cluster_cells(cds, reduction_method = 'UMAP')
plot_cells(cds, color_cells_by = "partition")
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "cluster",
           label_principal_points=T, label_branch_points=T,  group_label_size = 6)
cds <- order_cells(cds, root_pr_nodes='Y_76', reduction_method = 'UMAP') #Y_76 is the Monocytes
cds <- order_cells(cds, root_pr_nodes='Y_133', reduction_method = 'UMAP') #Y_133 is the cDC2


?plot_cells


## NEW - run this (takes a bit of time) to get the top genes
pr_graph_test <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)

# q value is statistic - so grabbing most signif
top_genes <- subset(pr_graph_test, q_value < 0.001)

# Morans I range between -1 and 1. If closer to 1, you have a very strong 'tightness' or 'smoothness' correlating with pseudotime. I.e. changes along pseudotime are very consistent.
top_genes2 <- top_genes[rev(order(top_genes$morans_I)),]

# plot top 20 genes. 
cds_subset <- cds[rownames(head(top_genes2, 20)), ]
plot_genes_in_pseudotime(cds_subset, color_cells_by = "pseudotime", ncol=4)

?plot_genes_in_pseudotime




hist(cds$pseudotime, breaks = 30)
?hist

cds@principal_graph_aux$UMAP$pseudotime


data$pseudo <- c(cds@principal_graph_aux$UMAP$pseudotime)

FeaturePlot(subset(data, subset = pseudo != Inf), "pseudo")


{
gene="TIMP1"
FetchData(data, c("pseudo", gene, 'ident')) %>% 
  mutate(kill = ifelse(is.infinite(pseudo), 0,1))%>%
  filter(kill !=0)->ddf
ddf$feature = ddf[,gene]
ggplot(ddf,aes_string("pseudo", "feature"))+geom_point(aes(color=ident))+geom_smooth(color='black')+theme_pubr()+ggtitle(gene)
}

{
  gene="FCGR3A"
  FetchData(data, c("pseudo", gene, 'ident')) %>% 
    mutate(kill = ifelse(is.infinite(pseudo), 0,1))%>%
    filter(kill !=0)->ddf
  ddf$feature = ddf[,gene]
  ggplot(ddf,aes_string("pseudo", "feature"))+geom_point(aes(color=ident))+geom_smooth(color='black')+theme_pubr()+ggtitle(gene)
}

{
  gene="LYVE1"
  FetchData(data, c("pseudo", gene, 'ident')) %>% 
    mutate(kill = ifelse(is.infinite(pseudo), 0,1))%>%
    filter(kill !=0)->ddf
  ddf$feature = ddf[,gene]
  ggplot(ddf,aes_string("pseudo", "feature"))+geom_point(aes(color=ident))+geom_smooth(color='black')+theme_pubr()+ggtitle(gene)
}

{
  gene="SIGLEC1"
  FetchData(data, c("pseudo", gene, 'ident')) %>% 
    mutate(kill = ifelse(is.infinite(pseudo), 0,1))%>%
    filter(kill !=0)->ddf
  ddf$feature = ddf[,gene]
  ggplot(ddf,aes_string("pseudo", "feature"))+geom_point(aes(color=ident))+geom_smooth(color='black')+theme_pubr()+ggtitle(gene)
}

{
  gene="HLA-DRB1"
  FetchData(data, c("pseudo", gene, 'ident')) %>% 
    mutate(kill = ifelse(is.infinite(pseudo), 0,1))%>%
    filter(kill !=0)->ddf
  ddf$feature = ddf[,gene]
  ggplot(ddf,aes_string("pseudo", "feature"))+geom_point(aes(color=ident))+geom_smooth(color='black')+theme_pubr()+ggtitle(gene)
}

check(data, "XCR")

pseudot <- function(obj = data, genes = c("LYVE1")) {
  if(length(genes) ==1) {
    ddf <- FetchData(data, c("pseudo", gene, 'ident')) %>% 
      mutate(kill = ifelse(is.infinite(pseudo), 0,1))%>%
      filter(kill !=0)
    ddf$feature = ddf[,gene]
    ggplot(ddf,aes_string("pseudo", "feature"))+geom_point(aes(color=ident))+geom_smooth(color='black')+theme_pubr()+ggtitle(gene)
    
  } else if(length(genes) >1) {
    ddf <- FetchData(data, c("pseudo", gene, 'ident')) %>% 
      mutate(kill = ifelse(is.infinite(pseudo), 0,1))%>%
      filter(kill !=0) %>%
      pivot_longer(cols=genes, names_to = "feature", values_to = "value")
    ddf$feature = ddf[,gene]
    ggplot(ddf,aes_string("pseudo", "value"))+geom_point(aes(color=ident))+geom_smooth(color='black')+theme_pubr()+ggtitle(paste("Genes:", paste0(genes,collapse = ", ")))+facet_wrap(~feature)
    
  }
    
}
pseudot(data, c("ITGAM"))

tmp <- subset(data, cells = Cells(data)[!is.infinite(data$pseudo)])       
pseudotime <- tmp@meta.data$pseudo
coords <- tmp@reductions$umap@cell.embeddings[, 1:2]  # or wherever your spatial data is
neighbors <- spdep::knearneigh(coords, k = 10)  # Adjust k as needed
nb <- spdep::knn2nb(neighbors)

spdep::knearneigh

# Create the spatial weights list
listw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
moran_test <- spdep::moran.test(pseudotime, listw)
print(moran_test)

# Example to visualize the pseudotime values
data$pseudotime <- pseudotime  # Store back to the Seurat object
FeaturePlot(data, features = "pseudotime")


install.packages("destiny")
library(destiny)

BiocManager::install("destiny")

# Create a matrix using scaled data
destiny_data <- GetAssayData(data, slot = "data")
is_matrix <- is.matrix(destiny_data)
if (!is_matrix) {
  destiny_data <- as.matrix(destiny_data)
}



# Run Destiny

diffusion_result <- DiffusionMap(destiny_data)

# Check the output

summary(diffusion_result)

# Visualize diffusion map
plot(diffusion_result)



saveRDS(data, "C:/Users/erica/OneDrive - The University of Sydney (Staff)/PhD - HIV/RNA Analysis/Bujko/ForTom.rds")



ls()
sapply(ls(), function(x) class(get(x)))

if (is.null(cds)) {
  stop("CellDataSet 'cds' is NULL. Ensure it's correctly initialized.")
}

rm(list = ls()[sapply(ls(), function(x) is.null(get(x)))])

save.image(file = "enviropack3.RData")
