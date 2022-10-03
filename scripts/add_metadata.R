add_metadata <- function(object, protospacers=NULL, scrub_doublets=TRUE) {


    guides <- fread(protospacers)[num_features == 1]
    guides[, target := tstrsplit(feature_call, '_')[1]]

    # add targets
    t <- guides[, target]
    
    # add NT
    nt <- guides[, fifelse(target %chin% 'NEG', 'NT', target)]
    
    # add guides
    g <- guides[, fifelse(feature_call %like% 'NEG', 'NT', feature_call)]

    # add status
    s <- guides[, fifelse(target %chin% 'NEG', 'non-targeting', 'perturbed')]
    
    names(g) <- names(t) <- names(s) <- names(nt) <- guides[, cell_barcode]
    
    
    object <- object %>% 
                subset(cells=guides[, cell_barcode]) %>% 

                AddMetaData(metadata=factor(t), col.name='genes') %>% 
                AddMetaData(metadata=factor(g), col.name='guide_RNA') %>% 
                AddMetaData(metadata=factor(s), col.name='Status') %>% 
                AddMetaData(metadata=factor(nt), col.name='NT')

    object[['percent.mt']] <- PercentageFeatureSet(object, pattern='^MT-', assay='RNA')
    
    
    doublet_rate <- (ncol(object) / 1000) * 0.008

    if (scrub_doublets) object <- object %>% scrublet(expected_doublet_rate=doublet_rate)


    return(object)
}
