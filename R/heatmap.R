library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(scales)
library(reshape2)
library(tidyverse)

figs_dir = '../figs/paper/'

so = readRDS('../processed/merged/so.rds')

g = DimPlot(so, group.by = 'state', cols=c('#F59546', '#48C133','#4087D8'),
            pt.size=0.05, shuffle=T)

ggsave(paste(figs_dir, "fig1/1b.pdf", sep=''), plot=g)

markers = FindAllMarkers(so, logfc.threshold = 0.25, min.pct = 0.1, only.pos = TRUE, test.use = 'roc')

heatmap_genes = (markers %>% group_by(cluster) %>% top_n(5, avg_log2FC))$gene

obj = so
features = c(c('Siglech', 'Cx3cr1'), heatmap_genes)
nclusters = max(as.integer(obj@meta.data$seurat_clusters))
avgexpr = AverageExpression(obj, features = features, group.by='seurat_clusters', assay='RNA')

mat = apply(avgexpr$RNA, 1, rescale)

g = ggplot(melt(mat), aes(Var1, Var2, fill=value)) +
  geom_tile(color='white', lwd=0.1, linetype=1)  +
  scale_x_continuous(limits = c(-1,nclusters), breaks = round(seq(0,nclusters))) + 
  theme(axis.text.x=element_text(size=6), axis.text.y=element_text(size=6, hjust = 0))  + coord_equal() +
  xlab('Cluster') + ylab('Gene') + scale_y_discrete(expand=c(0, 0)) +  scale_y_discrete(limits=rev) +
  scale_fill_gradientn(colors=rev(brewer.pal(11, 'RdYlBu')))

ggsave(paste(figs_dir, 'new_heatmap.pdf', sep=''), plot=g)
