working_dir <- '/data/ShernData/SWISNF/'; setwd(working_dir)

source('scripts/main/load_packages.r')

samples <- fread(snakemake@input[['samples']])[[1]]

filtered.path <- sapply(samples, function(lane) paste0(data_path, lane, '/filtered_feature_bc_matrix.h5'), simplify=FALSE) 

filtered.input.list <- sapply(filtered.path, Read10X_h5, simplify=FALSE)
filtered.input.list <- sapply(filtered.input.list, function(input.list) input.list[['Gene Expression']], simplify=FALSE)

protospacer_calls <- sapply(samples, function(lane) {
    
    pc <- fread(paste0(data_path, lane, '/protospacer_calls_per_cell.csv'))
    pc <- cSplit(pc, c('feature_call', 'num_umis'), direction='long', sep='|')
    pc[pc[, .I[which.max(num_umis)], cell_barcode]$V1]

}, simplify=FALSE)



object.list <- sapply(samples, function(lane) {

    object <- CreateSeuratObject(filtered.input.list[[lane]], min.cells=0, min.features=0, project=lane)

    object[['percent.mt']] <- PercentageFeatureSet(object, pattern='^MT-')
    object[['percent.rb']] <- PercentageFeatureSet(object, pattern='^RP[SL]')

    pc <- protospacer_calls[[lane]]

    pc[, target := tstrsplit(feature_call, '_')[1]]
    pc[, target := fifelse(target %chin% 'Human', 'CONTROL', target)]

    t <- pc[, target]
    g <- pc[, feature_call]
    s <- pc[, fifelse(target %chin% 'CONTROL', 'non-targeting', 'perturbed')]

    names(t) <- names(s) <- names(g) <- pc[, cell_barcode]
    
    doublet_rate <- (ncol(object) / 1000) * 0.008

    object %>% 
        AddMetaData(metadata=factor(t), col.name='target') %>% 
        AddMetaData(metadata=factor(g), col.name='guide') %>% 
        AddMetaData(metadata=factor(s), col.name='status') %>%
        AddMetaData(metadata=factor(lane), col.name='lane') %>% 

        scrublet(n_prin_comps=30, expected_doublet_rate=doublet_rate)
    
}, simplify=FALSE)

object <- merge(x=object.list[[1]], y=object.list[-1])

saveRDS(object, snakemake@output[['seurat_object']])