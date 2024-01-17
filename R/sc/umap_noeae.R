library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(scales)
library(reshape2)
library(tidyverse)

figs_dir = '../figs/paper/'

june2 = readRDS('/servers/iss-corescratch/lp488/Luca_June2023/processed/new_objects/june2_so.rds')

temp_so = june2
ko_eae = rownames(temp_so@meta.data[as.array(temp_so@meta.data$eae == 'EAE' & temp_so@meta.data$ctrl == 'no'), ])
wt_eae = rownames(temp_so@meta.data[as.array(temp_so@meta.data$eae == 'EAE' & temp_so@meta.data$ctrl == 'yes'), ])

names = c()

for (rowname in rownames(gg)) {
    if (rowname %in% ko_eae) {
        names = append(names, 'KO')
    } else if (rowname %in% wt_eae) {
        names = append(names, 'WT')
    } else {
        names = append(names, 'No EAE')
    }
}

gg['color'] = names
desired_order = c('No EAE', 'KO', 'WT')
gg$color <- factor(gg$color, levels=desired_order)
gg <- gg[order(gg$color),]
g = ggplot(gg, aes(x=UMAP_1, y=UMAP_2, color=color)) + geom_point(alpha=0.5, size=0.6) + 
    scale_color_manual(values = c("WT" = "#1e81b0", "KO" = "#d62762", 'No EAE' = "#aaaaaa"))

ggsave(paste(figs_dir, "june2_umap.pdf", sep=''), plot=g)
