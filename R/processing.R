library(ggplot2)
library(Seurat)
library(dplyr)
library(harmony)
library(data.table)
library(ShinyCell)

metadata = data.frame(
  'ID'= c('SIGAF10', 'SIGAB10', 'SIGAC10', 'SIGAB1', 'SIGAG12', 'SIGAD10', 'SIGAC1', 'SIGAD3',
          'SIGAA3', 'SIGAE12', 'SIGAE3', 'SIGAB3', 'SIGAF12'),
  'rep'= c('2','1','2','1','3','2','1','3', '2', '1', '3','2','1'),
  'batch' = c('10', '10', '10', '1', '12', '10', '1', '3', '3', '12', '3', '3', '12'),
  'stage' = c('peak', 'peak', 'peak', 'peak', 'chronic', 'peak', 'peak',
              'chronic', 'chronic', 'chronic', 'chronic', 'chronic', 'chronic'),
  'condition' = c('eae', 'eae', 'eae', 'eae', 'ctrl', 'ctrl', 'ctrl', 'eae', 'eae', 'eae', 'eae', 'eae', 'eae'),
  'stage_control' = c('peak', 'peak', 'peak', 'peak', 'ctrl', 'ctrl', 'ctrl',
              'chronic', 'chronic', 'chronic', 'chronic', 'chronic', 'chronic'),
  'fl' = c('resident', 'resident', 'infiltrating', 'infiltrating', 'control', 'control', 'control', 'resident',
           'resident', 'resident', 'infiltrating', 'infiltrating', 'infiltrating'),
  'fl_cre' = c('CreRFP', 'CreRFP', 'CreYFP', 'CreYFP', 'CreRFP', 'CreRFP', 'CreRFP', 'CreRFP',
          'CreRFP', 'CreRFP', 'CreYFP', 'CreYFP', 'CreYFP'),
  'state' = c('A-EAE', 'A-EAE', 'A-EAE', 'A-EAE', 'Ctrl', 'Ctrl', 'Ctrl', 'C-EAE', 'C-EAE', 'C-EAE', 'C-EAE', 'C-EAE', 'C-EAE')
)

full_names = vector()
for (ix in seq(nrow(metadata))) {
  row = metadata[ix,]
  full_names = append(full_names, paste(row$stage, '_', row$condition,j '_', row$fl, '_rep', row$rep, sep=''))
}
metadata$full_name = full_names

for (row_ix in seq(nrow(metadata))) {
  row = metadata[row_ix, ]
  if (row$ID == "SIGAG10") {
    next
  }
  path = row$ID
  sample.so = Read10X_h5(paste('../cellranger',
                               path,
                               'outs/filtered_feature_bc_matrix.h5',
                               sep='/'))
  
  sample.so = CreateSeuratObject(sample.so, project = row$full_name)
  sample.so$ID = row$ID
  sample.so$condition = row$condition
  sample.so$rep = row$rep
  sample.so$full_name = row$full_name
  sample.so$batch = row$batch
  sample.so$fl = row$fl
  sample.so$fl_cre = row$fl_cre
  sample.so$stage = row$stage
  sample.so$state = row$state
  sample.so$stage_control = row$stage_control
  
  if (row_ix == 1) {
    so = sample.so
  } else {
    so = merge(so, y = c(sample.so), add.cell.ids = c('', row$full_name), project = "Luca")
  }
}

mt_genes = grep('^mt-', row.names(so), value=TRUE)
rp_genes = grep('^Rp[sl]', row.names(so), value=TRUE)

so$percent_mt = PercentageFeatureSet(so, features=mt_genes)
so$percent_rp = PercentageFeatureSet(so, features=rp_genes)

so_ours = subset(so,
                   nFeature_RNA > 400 &
                     nFeature_RNA < 8000 &
                     nCount_RNA > 800 &
                     nCount_RNA < 35000 &
                     percent_rp > 3 &
                     percent_rp < 30 &
                     percent_mt < 5
)

nfeatures = 5000
ndims = 50

set.seed(0)

so_ours <- NormalizeData(so_ours)
so_ours <- FindVariableFeatures(so_ours, nfeatures = nfeatures)
so_ours <- ScaleData(so_ours, features = so_ours@assays$RNA@var.features)
so_ours <- RunPCA(so_ours, npcs = ndims, verbose = FALSE, approx = FALSE,
                  features = so_ours@assays$RNA@var.features[1:nfeatures])
so_ours <- RunHarmony(so_ours, group.by.vars = 'batch')
so_ours <- RunUMAP(so_ours, reduction = "harmony", dims = 1:ndims)

so_ours <- FindNeighbors(so_ours, reduction = "harmony", dims = 1:ndims)
so_ours <- FindClusters(so_ours, resolution=0.3)

# Cluster renaming
so_ours$seurat_clusters[so_ours$seurat_clusters == 11] = 2
so_ours$seurat_clusters[so_ours$seurat_clusters == 4] = 11
so_ours$seurat_clusters[so_ours$seurat_clusters == 5] = 4
so_ours$seurat_clusters[so_ours$seurat_clusters == 11] = 5
so_ours$seurat_clusters[so_ours$seurat_clusters == 12] = 11
so_ours$seurat_clusters[so_ours$seurat_clusters == 13] = 12
so_ours = SetIdent(so_ours, value='seurat_clusters')


markers = FindAllMarkers(so_ours, logfc.threshold = 0.25, min.pct = 0.1,
                             only.pos = TRUE, test.use = 'roc')


gprofiler_list = list()

for (cluster_str in unique(markers$cluster)) {
    idx = as.numeric(cluster_str) + 1
    gprofiler_results = gprofiler2::gost((markers %>% filter(cluster == cluster_str))$gene,
                                     organism='mmusculus',
                                     sources=c('GO:BP'),
                                     correction_method = 'fdr',
                                     evcodes=T)
    gprofiler_list[[idx]] = gprofiler_results$result
}

for (cluster_str in unique(markers$cluster)) {
    print(cluster_str)
    write.csv(markers %>% filter(cluster == cluster_str), paste('../processed//clusters/markers_' , cluster_str, '.csv', sep=''))
    idx = as.numeric(cluster_str) + 1
    if (typeof(gprofiler_list[[idx]]) == 'NULL') {
        next
    }
    fwrite(gprofiler_list[[idx]], paste('../processed/clusters/go_' , cluster_str, '.csv', sep=''))
}

saveRDS(so_ours, file.path('../processed/so.rds'))

scConf = createConfig(so_ours)
makeShinyApp(so_ours, scConf, gene.mapping = TRUE,
             shiny.title = "Clusters", shiny.dir = "../processed/shinyapp/")

