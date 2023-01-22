working_dir <- '/data/ShernData/SWISNF/'; setwd(working_dir)

source('scripts/main/load_packages.r')

plan('multicore', workers=snakemake@threads)
options(future.globals.maxSize=50 * 1000 * 1024^2)


object.list <- readRDS(snakemake@input[['seurat_object']])

noise <- c('nCount_RNA', 'S.Score', 'G2M.Score', 'nFeature_RNA', 'percent.mt', 'percent.rb')

object.list <- sapply(object.list, CleanVarGenes, nHVG=3000, simplify=FALSE)

features <- object.list %>% SelectIntegrationFeatures(nfeatures=5000)

object.list <- sapply(object.list, function(object) {
    object %>% 
      ScaleData(features=features, vars.to.regress=noise, verbose=FALSE) %>% 
      RunPCA(features=features, verbose=FALSE)
})

anchors <- object.list %>% FindIntegrationAnchors(anchor.features=features, reduction='rpca', k.anchor=7, scale=FALSE)
object <- IntegrateData(anchorset=anchors, new.assay.name='integrated', dims=1:30)

DefaultAssay(object) <- 'integrated'

object <- object %>% 
            ScaleData(verbose=FALSE, vars.to.regress='doublet_scores') %>% 
            RunPCA(verbose=FALSE)

Project(object) <- 'SWISNF'

saveRDS(object, snakemake@output[['seurat_object']])

