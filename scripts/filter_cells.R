filter_cells <- function(object=NULL, lower=0.02, upper=0.98, filter.by=NULL) {
  
  require(data.table); require(Seurat)

  metadata <- copy(object@meta.data)
  
  setDT(metadata, keep.rownames=TRUE)
  new_names <- c('cells', 'nUMI', 'nGene')
  old_names <- c('rn', 'nFeature_RNA', 'nCount_RNA')
  setnames(metadata, old=old_names, new=new_names)
  
  cells2keep <- metadata[fintersect(fintersect(metadata[, 
                                    .(cellNo=.I[between(nUMI, quantile(nUMI, lower), quantile(nUMI, upper))]), 
                                    by=filter.by], 
                           metadata[, 
                                    .(cellNo=.I[between(nGene, quantile(nGene, lower), quantile(nGene, upper))]), 
                                    by=filter.by]), metadata[, 
                                                            .(cellNo=.I[between(percent.mt, quantile(percent.mt, lower), quantile(percent.mt, upper))]), 
                                                            by=filter.by])[, cellNo]][, cells]
  
  
  
  return(object %>% 
    subset(cells=cells2keep))

}



