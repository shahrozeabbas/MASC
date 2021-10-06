args <- commandArgs(trailingOnly=TRUE)

#set wd and load packages
home <- '/data/ShernData/'
setwd('/data/ShernData/PerturbSeqData/')

pkgs <- c('data.table', 'Seurat', 'dplyr', 'ggplot2', 'future', 'reticulate', 'parallel')
invisible(sapply(pkgs, require, character.only=TRUE))

output <- 'default_output'

h <- 30
w <- 45
options(ggrepel.max.overlaps=Inf)
plan('multicore', workers=availableCores() / 2)
options(future.globals.maxSize=100 * 1000 * 1024^2)

invisible(lapply(X=list.files(path=paste0(home, 'scripts_mixscape'),
                                pattern='.R', full.names=TRUE), FUN=source))

use_python(system('which python', intern=TRUE), required=TRUE)
source_python(paste0(home, 'scripts_NF1/scrublet_py.py'))


swisnf_data <- list.files('SWISNF', full.names=TRUE)

swisnf <- mclapply(swisnf_data, function(lane) {
    setwd(lane)

    data <- Read10X_h5('filtered_feature_bc_matrix.h5', use.names=TRUE, unique.features=TRUE)
    object <- CreateSeuratObject(counts=data[['Gene Expression']], project='SWISNF', min.cells=0, min.features=0)
    object <- object %>% add_metadata(protospacers='protospacer_calls_per_cell.csv')

}, mc.cores=length(swisnf_data))

swisnf <- merge(x=swisnf[[1]], y=swisnf[2:3], add.cell.ids=c('1', '2', '3'), project='SWISNF')
swisnf <- swisnf %>% AddMetaData(metadata=rep('SWISNF', ncol(swisnf)), col.name='exP')


prc1_data7 <- list.files('PRC1', full.names=TRUE)[1:2]

PRC1.data <- mclapply(prc1_data7, function(lane) {
    setwd(lane)

    data <- Read10X_h5('filtered_feature_bc_matrix.h5', use.names=TRUE, unique.features=TRUE)
    object <- CreateSeuratObject(counts=data[['Gene Expression']], project='PRC1', min.cells=0, min.features=0)
    object <- object %>% add_metadata(protospacers='protospacer_calls_per_cell.csv')

}, mc.cores=length(prc1_data7))

PRC1.data <- merge(x=PRC1.data[[1]], y=PRC1.data[2], add.cell.ids=c('4', '5'), project='PRC1_Day7')
PRC1.data <- PRC1.data %>% AddMetaData(metadata=rep('PRC1_Day7', ncol(PRC1.data)), col.name='exP')


prc1_data10 <- list.files('PRC1', full.names=TRUE)[3:4]

PRC1.data2 <- mclapply(prc1_data10, function(lane) {
    setwd(lane)

    data <- Read10X_h5('filtered_feature_bc_matrix.h5', use.names=TRUE, unique.features=TRUE)
    object <- CreateSeuratObject(counts=data[['Gene Expression']], project='PRC1', min.cells=0, min.features=0)
    object <- object %>% add_metadata(protospacers='protospacer_calls_per_cell.csv')

}, mc.cores=length(prc1_data10))

PRC1.data2 <- merge(x=PRC1.data2[[1]], y=PRC1.data2[2], add.cell.ids=c('6', '7'), project='PRC1_Day10')
PRC1.data2 <- PRC1.data2 %>% AddMetaData(metadata=rep('PRC1_Day10', ncol(PRC1.data2)), col.name='exP')


##Integrate Seurat Objects
integration_list <- list(swisnf, PRC1.data, PRC1.data2)

#normalize objects and find variable featyres
integration_list <- lapply(X=integration_list, FUN=function(x) {
  x <- x %>%
        filter_cells() %>%
        subset(subset=predicted_doublets == 'singlet') %>%

        NormalizeData() %>%
        FindVariableFeatures()
})

#select features for downstream integration and run PCA on each object in the list
# k <- 20
k <- as.numeric(args[1])
noise <- c('nCount_RNA', 'nFeature_RNA', 'percent.mt')
features <- SelectIntegrationFeatures(object.list=integration_list)

integration_list <- lapply(X=integration_list, FUN=function(x) {
  x <- x %>%
        ScaleData(features=features, vars.to.regress=noise, verbose=FALSE) %>%
        RunPCA(features=features, verbose=FALSE)
})

#Find integration anchors
anchors <- FindIntegrationAnchors(integration_list, anchor.features=features, reduction='rpca', k.anchor=k, scale=FALSE)

combined <- IntegrateData(anchorset=anchors)
DefaultAssay(combined) <- 'integrated'

combined <- combined %>%
  ScaleData(verbose=FALSE) %>%
  RunPCA(npcs=30, verbose=FALSE) %>%
  RunUMAP(dims=seq(30))


invisible(lapply(c('exP', 'genes'), function(x) {
    p <- combined %>% DimPlot(reduction='umap', group.by=x, pt.size=4)
    ggsave(plot=p, height=h, width=w, filename=paste0(output, '/plots/', DefaultAssay(combined), '_', x, '_umap.pdf'))
}))


#Re-visualise metrics after subset
# v <- VlnPlot(combined, features=noise, ncol=3)
v <- VlnPlot(combined, features=c(noise, 'doublet_scores'), ncol=4, group.by='orig.ident')
ggsave(plot=v, height=h, width=w, filename=paste0(output, '/plots/noise_vln.pdf'))

plot1 <- FeatureScatter(combined, feature1='nCount_RNA', feature2='percent.mt')
plot2 <- FeatureScatter(combined, feature1='nCount_RNA', feature2='nFeature_RNA')
plot3 <- FeatureScatter(combined, feature1='nCount_RNA', feature2='doublet_scores')

scatter <- plot1 + plot2 + plot3
ggsave(plot=scatter, height=h, width=w, filename=paste0(output, '/plots/', DefaultAssay(combined), '_feature_scatter.pdf'))


##Run perturbation signature calculation function
combined <- combined %>%
                CalcPerturbSig(assay='RNA', slot='data', gd.class ='NT', nt.cell.class='NT',
                    reduction='pca', ndims=30, num.neighbors=20, split.by='exP', new.assay.name='PRTB')

# Prepare PRTB assay for dimensionality reduction:
# Normalize data, find variable features and center data.
DefaultAssay(combined) <- 'RNA'
combined <- combined %>% NormalizeData() %>% FindVariableFeatures(nfeatures=5000)

DefaultAssay(object=combined) <- 'PRTB'
# Use variable features from RNA assay.
VariableFeatures(object=combined) <- VariableFeatures(object=combined[['RNA']])

combined <- combined %>%

                ScaleData(do.scale=FALSE, do.center=TRUE) %>%
                RunPCA(reduction.key='prtbpca', reduction.name='prtbpca') %>%
                RunUMAP(dims=1:40, reduction='prtbpca', reduction.key='prtbumap', reduction.name='prtbumap') %>%

                RunMixscape(assay='PRTB', slot='scale.data', labels='NT', nt.class.name='NT',
                            new.class.name='mixscape_class', min.de.genes=5, iter.num=10, de.assay='RNA')

saveRDS(combined, file=paste0(output, '/objects/combined_integrated_mixscape.rds'))


Idents(combined) <- 'guide_RNA'
df <- prop.table(table(combined$mixscape_class.global, combined$guide_RNA), 2)

df2 <- reshape2::melt(df)
df2$Var2 <- as.character(df2$Var2)
test <- df2[which(df2$Var1 == 'KO'),]
test <- test[order(test$value, decreasing=T),]
new.levels <- test$Var2
df2$Var2 <- factor(df2$Var2, levels=new.levels )
df2$Var1 <- factor(df2$Var1, levels=c('NT', 'NP', 'KO'))
df2$gene <- sapply(as.character(df2$Var2), function(x) strsplit(x, split='_')[[1]][1])
df2$guide_number <- sapply(as.character(df2$Var2),
                           function(x) strsplit(x, split='_')[[1]][2])
df3 <- df2[-c(which(df2$gene == 'NT')),]

p1 <- ggplot(df3, aes(x=guide_number, y=value*100, fill= Var1)) +
  geom_bar(stat= 'identity') +
  theme_classic()+
  scale_fill_manual(values=c('grey49', 'grey79','coral1')) +
  ylab('% of cells') +
  xlab('sgRNA')

p1 <- p1 + theme(axis.text.x=element_text(size=18, hjust=1),
           axis.text.y=element_text(size=18),
           axis.title=element_text(size=16),
           strip.text=element_text(size=16, face='bold')) +
  facet_wrap(vars(gene),ncol=5, scales='free') +
  labs(fill='mixscape class') + theme(legend.title=element_text(size=14),
                                       legend.text=element_text(size=12))

ggsave(plot=p1, height=w, width=h, filename=paste0(output, '/plots/ko_barplot.pdf'))

#sub summary characteristics
Idents(combined) <- 'genes'
All_sub <- count(combined@meta.data, vars=genes)
fwrite(All_sub, paste0(output, '/summaries/intergrated_counts.csv'))

# Remove non-perturbed cells and run LDA to reduce the dimensionality of the data.
Idents(combined) <- 'mixscape_class.global'
sub <- subset(combined, idents=c('KO', 'NT'))

sub <- sub %>%
            MixscapeLDA(assay='RNA', pc.assay='PRTB', labels='NT', nt.label='NT', npcs=10, logfc.threshold=0.25, verbose=FALSE) %>%
            RunUMAP(dims=1:79, reduction='lda', reduction.key='ldaumap', reduction.name='ldaumap')

saveRDS(sub, file=paste0(output, '/objects/sub_ko_cells.rds'))

# Visualize UMAP clustering results.
p <- DimPlot(sub, reduction='ldaumap', label=TRUE, repel=TRUE, label.size=7, group.by='mixscape_class', pt.size=2) + NoLegend()
ggsave(plot=p, height=h, width=w, filename=paste0(output, '/plots/mixscape_class_group_umap.pdf'))

p <- DimPlot(sub, reduction='ldaumap', split.by='mixscape_class', ncol=8, label=TRUE)
ggsave(plot=p, height=h, width=w, filename=paste0(output, '/plots/mixscape_class_split_umap.pdf'))


#sub summary characteristics
Idents(sub) <- 'mixscape_class'
All_sub <- count(sub@meta.data, vars=mixscape_class)
fwrite(All_sub, paste0(output, '/summaries/integrated_sub_counts.csv'))


genes <- unique(sub@meta.data$NT)
genes <- genes[genes != 'NT']

# out <- vector('list')
for (gene in genes) {

  markers <- FindMarkers(sub, ident.1='NT', ident.2=gene, group.by='NT')

  setDT(markers, keep.rownames=TRUE)
  m <- copy(markers[p_val_adj < 0.05])

  fwrite(m, file=paste0(output, '/de_genes/signif_', gene, '.csv'), row.names=TRUE)
  fwrite(markers, file=paste0(output, '/de_genes/All_', gene, '.csv'), row.names=TRUE)

}
