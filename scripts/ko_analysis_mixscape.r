pkgs <- c('data.table', 'Seurat', 'dplyr', 'ggplot2', 'future', 'reticulate', 'wesanderson')
invisible(sapply(pkgs, require, character.only=TRUE))

home <- '/gpfs/gsfs12/users/abbass2/'

custom_theme <- theme(
  plot.title=element_text(size=16, hjust=0.5), 
  legend.key.size=unit(0.7, 'cm'), 
  legend.text=element_text(size=14))

guide_colors <- c('chartreuse3', 'cornflowerblue', 'darkorchid2', 'lightcoral', 'turquoise1', 
                  'yellow3', 'deeppink2', 'seagreen1', 'orange2', 'royalblue4', 
                  'maroon', 'blue', 'black')

cores <- availableCores() / 2
plan('multicore', workers=cores)
setDTthreads(threads=1)

use_python(system('which python', intern=TRUE), required=TRUE)
source_python(paste0(home, 'scrublet_py.py'))

invisible(lapply(X=list.files(path=paste0(home, 'scripts_mixscape'), 
                                pattern='.R', full.names=TRUE), FUN=source)) 

data_dir <- paste0(home, 'SNF1CrispriData/SCAF2095/')

setwd(data_dir)
data <- CreateSeuratObject(counts=Read10X_h5(list.files(pattern='h5'))[['Gene Expression']],
                           project='pilot_dox', min.cells=0, min.features=0)

d <- 60
ngbs <- 100
nhvg <- 5000
assay <- 'RNA'
k <- round(0.5 * sqrt(ncol(data)))
guides <- fread('barcode_guide_nUMI_target.csv')
options(future.globals.maxSize=ngbs * 1000 * 1024^2)
noise <- c('nFeature_RNA', 'nCount_RNA', 'percent.mt', 'n_guides')


var_method <- 'mvp'
use_top_var_genes <- TRUE

data <- data %>% 
    
    add_metadata(metadata=guides, ncores=cores, scrub_doublets=TRUE) %>% 
    filter_cells() %>% NormalizeData() %>% 
    
    { if (use_top_var_genes) FindVariableFeatures(., selection.method=var_method, nfeatures=nhvg) else . } 

genes <- if (use_top_var_genes) VariableFeatures(data) else rownames(data)

data <- data %>% 

    ScaleData(vars.to.regress=noise, features=genes) %>% 
    
    RunPCA(features=genes, npcs=d) %>% 
    
    CalcPerturbSig(assay=assay, slot='data', gd.class='targets', nt.cell.class='NT', 
                    reduction='pca', num.neighbors=k, ndims=d, new.assay.name='PRTB', 
                    features=genes)


DefaultAssay(data) <- 'PRTB'

data <- data %>% 
    
    ScaleData(do.scale=FALSE, do.center=TRUE, features=genes) %>%
    
    RunPCA(reduction.key='prtbpca', reduction.name='prtbpca', features=genes, npcs=d) %>%
    
    RunUMAP(dims=seq(d), reduction='prtbpca', n.neighbors=k,
            reduction.key='prtbumap', reduction.name='prtbumap') %>%
    
    RunMixscape(assay='PRTB', slot='scale.data', labels='targets', 
                nt.class.name='NT', new.class.name='mixscape_class', min.de.genes=2, 
                logfc.threshold=0.1, iter.num=10, de.assay=assay, verbose=FALSE, 
                fine.mode=FALSE, fine.mode.labels='guide_RNA', min.cells=1)


Idents(data) <- 'mixscape_class.global'
sub <- subset(data, idents=c('KO', 'NT'))
Project(sub) <- paste0('ko_', Project(data))

sub <- sub %>% 
    MixscapeLDA(assay=assay, pc.assay='PRTB', labels='targets', nt.label='NT', 
                npcs=10, logfc.threshold=0.01, verbose=FALSE) %>% 
    
    RunUMAP(dims=seq(8), reduction='lda', reduction.key='ldaumap', reduction.name='ldaumap')


out_path <- paste0(data_dir, '/metadata/')

Idents(sub) <- 'mixscape_class'
sub$mixscape_class <- as.factor(sub$mixscape_class)

if (length(levels(sub$mixscape_class)) < length(guide_colors)) {
  g_c <- sample(x=guide_colors, size=length(levels(sub$mixscape_class)))
} 

names(g_c) <- levels(sub$mixscape_class)

ko_umap <- sub %>% 
      DimPlot(reduction='ldaumap', label=TRUE, repel=TRUE, label.size=5, pt.size=0.5) +
          ylab('LDAUMAP 2') + xlab('LDAUMAP 1') + custom_theme + scale_color_manual(values=g_c, drop=FALSE)  

ggsave(filename=paste0(out_path, 'final_', Project(sub), '_umap.pdf'), plot=ko_umap, width=14, height=10)

data %>% plots_outputs(out_dir=out_path, group='targets', plot=TRUE, markers=FALSE)
data %>% SaveH5Seurat(filename=paste0(out_path, 'final_', Project(data), '.h5Seurat'), overwrite=TRUE)

# sub %>% plots_outputs(group='targets', plot=TRUE, markers=FALSE)
# sub %>% SaveH5Seurat(filename=paste0(out_path, 'final_', Project(sub), '.h5Seurat'), overwrite=TRUE)

message('Number of Cells: ', ncol(data))
message('Number of KO/NT Cells: ', ncol(sub))



# ##end of mixscape v3.1

# #in mixscape v3.2
# # Calculate percentage of KO cells for all target gene classes.
# df <- prop.table(table(data$mixscape_class.global, data$NT),2)

# df2 <- reshape2::melt(df)
# df2$Var2 <- as.character(df2$Var2)
# test <- df2[which(df2$Var1 == 'KO'),]
# test <- test[order(test$value, decreasing=T),]
# new.levels <- test$Var2
# df2$Var2 <- factor(df2$Var2, levels=new.levels )
# df2$Var1 <- factor(df2$Var1, levels=c('NT', 'NP', 'KO'))
# df2$gene <- sapply(as.character(df2$Var2), function(x) strsplit(x, split='g')[[1]][1])
# df2$guide_number <- sapply(as.character(df2$Var2), 
#                            function(x) strsplit(x, split='g')[[1]][2])
# df3 <- df2[-c(which(df2$gene == 'NT')),]

# p1 <- ggplot(df3, aes(x=guide_number, y=value*100, fill= Var1)) +
#   geom_bar(stat= 'identity') +
#   theme_classic()+
#   scale_fill_manual(values=c('red', 'grey79','blue')) + 
#   ylab('% of cells') +
#   xlab('sgRNA')

# p1 + theme(axis.text.x=element_text(size=18, hjust=1), 
#            axis.text.y=element_text(size=18), 
#            axis.title=element_text(size=16), 
#            strip.text=element_text(size=16, face='bold')) + 
#   facet_wrap(vars(gene),ncol=5, scales='free') +
#   labs(fill='mixscape class') +theme(legend.title=element_text(size=14),
#                                        legend.text=element_text(size=12))


# # Explore the perturbation scores of cells.
# PlotPerturbScore(object=data, 
#                  target.gene.ident='TMSB10',
#                  group.by='mixscape_class', 
#                  col='coral2') + labs(fill='mixscape class')

# # Inspect the posterior probability values in NP and KO cells.
# VlnPlot(data, 'mixscape_class_p_ko', idents=c('NT', 'TMSB10 KO', 'TMSB10 NP')) +
#   theme(axis.text.x=element_text(angle=0, hjust=0.5)) + 
#   NoLegend() +
#   ggtitle('mixscape posterior probabilities')



# # Run DE analysis and visualize results on a heatmap ordering cells by their posterior 
# # probability values.
# Idents(object=data) <- 'NT'
# MixscapeHeatmap(object=data, 
#                 ident.1='NT', 
#                 ident.2='ACTL6A', 
#                 balanced=F, 
#                 assay=assay, 
#                 logfc.threshold=0,
#                 max.genes=20, angle=0, 
#                 group.by='mixscape_class', 
#                 max.cells.group=300, 
#                 size=3.5) + NoLegend()

