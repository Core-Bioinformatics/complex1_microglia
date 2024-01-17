library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(scales)
library(reshape2)
library(tidyverse)

figs_dir = '../figs/paper/'

complex1 <- c('Ndufs1', 'Ndufs2', 'Ndufs3', 'Ndufs4', 'Ndufs5', 'Ndufs6', 'Ndufs7', 'Ndufs8',
              'Ndufv1', 'Ndufv2', 'Ndufv3', 'Ndufab1', 'Ndufa1', 'Ndufa2', 'Ndufa3', 'Ndufa5',
              'Ndufa6', 'Ndufa7', 'Ndufa8', 'Ndufa9', 'Ndufa10', 'Ndufa11', 'Ndufa12', 'Ndufa13',
              'Ndufb1', 'Ndufb2', 'Ndufb3', 'Ndufb4', 'Ndufb5', 'Ndufb6', 'Ndufb7', 'Ndufb8', 'Ndufb9',
              'Ndufb10', 'Ndufb11', 'Ndufc1', 'Ndufc2')
complex2 <- c('Sdha', 'Sdhb', 'Sdhc', 'Sdhd')

plot_genes<- function(list_of_genes, title, so, expr_threshold=0, relaxation=0) {
    expr_matrix <- as.matrix(so@assays$RNA@data[rownames(so@assays$RNA@data) %in% list_of_genes,])
    new_matrix = matrixStats::colMeans2(expr_matrix)
    a = so@reductions$umap@cell.embeddings
    a_2 = cbind(a, list(new_matrix)[[1]])
    colnames(a_2) = c('UMAP_1', 'UMAP_2', 'Mean gene set expression')
    df = as.data.frame(a_2)
    x = ggplot(df %>% arrange(`Mean gene set expression`), aes(x=UMAP_1, y=UMAP_2)) + geom_point(aes(color=`Mean gene set expression`), size=2) +
        scale_color_gradient2(low='#291043', mid='#1f8f8d', high='#fde726',  midpoint = 0.5) +
        ggtitle(title) + theme_bw()
    print(x)
    x
}

to_plot = list('Complex I'=complex1, 'Complex II'=complex2)

for (i in seq_along(to_plot)){
    l = to_plot[[i]]
    name = names(to_plot)[i]
    g = plot_genes(l, name, so_cluster) 
    ggsave(paste(figs_dir, "umap", "_", i, ".pdf", sep=""), plot=g, width=10)
}
