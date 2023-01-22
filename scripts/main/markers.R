working_dir <- '/data/ShernData/SWISNF/'; setwd(working_dir)

source('scripts/main/load_packages.r')

plan('multicore', workers=snakemake@threads)
options(future.globals.maxSize=50 * 1000 * 1024^2)


object <- readRDS(snakemake@input[['mixscape_object_object']])

genes <- unique(object@meta.data$target)
genes <- genes[genes != 'CONTROL']


de.genes <- rbindlist(lapply(genes, function(gene) {
    
    markers <- object %>% FindMarkers(ident.1='CONTROL', ident.2=gene, group.by='target')

    setDT(markers, keep.rownames='marker')
    
    markers[, target := gene]

    centers <- markers[, kmeans(x=avg_log2FC, centers=2)$centers]
    
    markers[avg_log2FC > min(centers) & avg_log2FC < max(centers), grade := 'low']
    markers[!(avg_log2FC > min(centers) & avg_log2FC < max(centers)), grade := 'mid']
    markers[avg_log2FC < quantile(avg_log2FC, 0.01) | avg_log2FC > quantile(avg_log2FC, 0.99), grade := 'high']

    markers
  
}))

fwrite(x=de.genes, file=snakemake@output[['markers']])