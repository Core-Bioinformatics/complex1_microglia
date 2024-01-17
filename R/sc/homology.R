library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(scales)
library(reshape2)
library(tidyverse)

figs_dir = '../figs/paper/'

so_absinta = readRDS("../processed/absinta.rds")
cremato_mini = readRDS('../processed/subcluster.rds')

absinta_microclia_clusters = c('6', '9', '20')
absinta_microglia = subset(so_absinta, seurat_clusters %in% absinta_microclia_clusters) 
a = as.data.frame(absinta_microglia@reductions$umap@cell.embeddings)
cells = rownames(a[a$UMAP_1 < 3 ,])
absinta_microglia = subset(absinta_microglia, cells = cells)
cells = rownames(a[a$UMAP_2 < -4,])
absinta_microglia = subset(absinta_microglia, cells = cells)

absinta_microglia <- FindNeighbors(absinta_microglia, reduction = "harmony")
absinta_microglia <- FindClusters(absinta_microglia, resolution = 1.4)

absinta_microglia@meta.data$in_mams = absinta_microglia@meta.data$seurat_clusters %in% c('3')
absinta_subcluster = subset(absinta_microglia, in_mams == T)

markers = read.csv('../processed/markers_4.csv')
cluster4_genes <- markers$gene
genes_to_consider = lapply(cluster4_genes, toupper)

rows2 = rownames(absinta_subcluster@assays$SCT@data) 
genes_to_index_absinta = intersect(genes_to_consider, rows2)

x = cremato_mini@assays$RNA@data[as.character(lapply(genes_to_index_absinta, str_to_title)), ]
df_cremato = as.data.frame(log(matrixStats::rowMeans2(as.matrix(x)) + 1))
colnames(df_cremato) = 'Cremato'
rownames(df_cremato) = lapply(rownames(df_cremato), toupper)

x = absinta_subcluster@assays$SCT@data[as.character(genes_to_index_absinta), ]
df_absinta = as.data.frame(log(matrixStats::rowMeans2(as.matrix(x)) + 1))
colnames(df_absinta) = 'Absinta'

cremato_absinta = merge(df_cremato, df_absinta, by = 0)

g = ggplot(cremato_absinta, aes(x=Cremato, y=Absinta)) + geom_point() + geom_smooth(method='lm')

ggsave(paste(figs_dir, 'schirmer_absinta/cremato_absinta_scatter.pdf', sep=""), plot=g)
