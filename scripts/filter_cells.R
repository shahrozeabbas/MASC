filter_cells <- function(object=NULL, lower=0.02, upper=0.98, filter.by=NULL) {
  
  require(data.table); require(Seurat)

  metadata <- copy(object@meta.data)
  
  setDT(metadata, keep.rownames=TRUE)
  new_names <- c('cells', 'nGene', 'nUMI')
  old_names <- c('rn', 'nFeature_RNA', 'nCount_RNA')
  setnames(metadata, old=old_names, new=new_names)
  
  metadata <- metadata[fintersect(fintersect(metadata[, 
                                    .(cellNo=.I[between(nUMI, quantile(nUMI, lower), quantile(nUMI, upper))]), 
                                    by=filter.by], 
                           metadata[, 
                                    .(cellNo=.I[between(nGene, quantile(nGene, lower), quantile(nGene, upper))]), 
                                    by=filter.by]), metadata[, 
                                                            .(cellNo=.I[between(percent.mt, quantile(percent.mt, lower), quantile(percent.mt, upper))]), 
                                                            by=filter.by])[, cellNo]]
  
  
  if (Project(object) %like% 'PRC1') {
    cells2keep <- metadata[nGene > 3000 & nGene < 7000][nUMI > 10000 & nUMI < 45000][percent.mt > 2, cells]
  } else {
    cells2keep <- metadata[, cells]
  }
  
  return(object %>% 
    subset(cells=cells2keep))

}



