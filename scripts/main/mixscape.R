working_dir <- '/data/ShernData/SWISNF/'; setwd(working_dir)

source('scripts/main/load_packages.r')

plan('multicore', workers=snakemake@threads)
options(future.globals.maxSize=50 * 1000 * 1024^2)


object <- readRDS(snakemake@input[['seurat_object']])

features <- rownames(object[['integrated']])

object <- object %>% 
            
            CalcPerturbSig(assay='RNA', slot='data', gd.class ='target', features=features,
            nt.cell.class='CONTROL', reduction='pca', ndims=30, num.neighbors=20, new.assay.name='PRTB')

DefaultAssay(object) <- 'PRTB'
VariableFeatures(object) <- features

object <- object %>% 
            
            ScaleData(do.scale=FALSE, do.center=TRUE) %>% 
            RunPCA(reduction.key='prtbpca', reduction.name='prtbpca') %>% 
            RunUMAP(dims=1:30, reduction='prtbpca', reduction.key='prtbumap', reduction.name='prtbumap') %>%
            
            RunMixscape(assay='PRTB', slot='scale.data', labels='target', nt.class.name='CONTROL', 
                        new.class.name='mixscape_class', de.assay='RNA', min.de.genes=snakemake@params[['min_de_genes']])


saveRDS(object, snakemake@output[['seurat_object']])


Idents(object) <- 'mixscape_class.global'
sub <- object %>% subset(idents=c('KO', 'CONTROL'))

dims <- length(table(sub$mixscape_class)) - 1

sub <- sub %>% 
            MixscapeLDA(assay='RNA', pc.assay='PRTB', labels='target', nt.label='CONTROL', 
                        npcs=10, logfc.threshold=snakemake@params[['log_FC']]) %>%

            RunUMAP(dims=1:dims, reduction='lda', reduction.key='ldaumap', reduction.name='ldaumap')


saveRDS(sub, snakemake@output[['mixscape_sub_object']])