working_dir <- '/data/ShernData/SWISNF/'; setwd(working_dir)

source('scripts/main/load_packages.r')

plan('multicore', workers=snakemake@threads)
options(future.globals.maxSize=50 * 1000 * 1024^2)


object <- readRDS(snakemake@input[['seurat_object']])

object.list <- object %>% 
    FilterCells() %>% NormalizeData() %>% SplitObject(split='lane') %>% PlotQC(project.name='SWISNF')


s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

object.list <- sapply(object.list, function(object) {
    genes <- rownames(object)

    object %>% 
        ScaleData(features=genes, verbose=FALSE) %>% 
        CellCycleScoring(s.features=s.genes, g2m.features=g2m.genes)

}, simplify=FALSE)


saveRDS(object.list, snakemake@output[['seurat_object']])