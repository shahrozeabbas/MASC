add_metadata <- function(object, metadata=NULL, ncores=1, scrub_doublets=FALSE, find_doublets=FALSE) {


    if (is.null(metadata)) {
        root <- '/gpfs/gsfs12/users/abbass2/BegaPerturbData/Seq1_GEXnFBC_wNTG/SCAF1819_Dox'
        guides <- fread(paste0(root, 'barcode_guide_nUMI_target.csv'))
    } else {
        guides <- metadata
    }

    mug <- copy(guides[guides[, .I[which.max(nUMI)], by=cell_barcode]$V1])

    multiple_guides <- guides[, .N, keyby=cell_barcode][N > 1, cell_barcode]

    mcd <- copy(guides[cell_barcode %chin% multiple_guides, 
                .(msd=sqrt(mean((nUMI - max(nUMI))^2)) / max(nUMI)), keyby=cell_barcode][
                    , .(cell_barcode, mcd=1 - msd)])

    mcd[, cluster := fifelse(kmeans(x=mcd, centers=2, 
                                nstart=round(length(mcd) * 0.05))$cluster == 1, 'singlet', 'multiplet')]


    mug <- mcd[mug, on='cell_barcode'][
    , `:=` (mcd=fifelse(is.na(mcd), 0.0, mcd), 
            cluster=fifelse(is.na(cluster), 'singlet', cluster))]

    # add meta object
    g <- mug[, guide]
    names(g) <- mug[, cell_barcode]

    t <- mug[, fifelse(target %chin% 'Human', 'NT', target)]
    names(t) <- mug[, cell_barcode]

    d <- mug[, mcd]
    names(d) <- mug[, cell_barcode]

    z <- guides[, .N, keyby=cell_barcode][, N]
    names(z) <- guides[, .N, keyby=cell_barcode][, cell_barcode]


    object <- object %>% 
                subset(cells=mug[, cell_barcode]) %>% 
                
                AddMetaData(metadata=d, col.name='mcd_score') %>%
                AddMetaData(metadata=factor(t), col.name='targets') %>% 
                AddMetaData(metadata=factor(g), col.name='guide_RNA') %>%
                AddMetaData(metadata=factor(z), col.name='n_guides')


    object[['percent.mt']] <- object %>% PercentageFeatureSet(pattern='^MT-')

    
    doublet_rate <- (ncol(object) / 1000) * 0.008

    if (scrub_doublets) object <- object %>% scrublet(expected_doublet_rate=doublet_rate)
    
    if (find_doublets) object <- object %>% doublet_finder(expected_doublet_rate=doublet_rate, cores=ncores)


    return(object)
}