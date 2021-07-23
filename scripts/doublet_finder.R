require(DoubletFinder)

doublet_finder <- function(object, npcs=50, sct=FALSE, GT=FALSE, expected_doublet_rate=0.06, cores=1, k=NULL) {

  if (is.null(k)) k <- round(0.5 * sqrt(ncol(object)))
  noise <- c('nFeature_RNA', 'nCount_RNA', 'percent.mt')

  object <- object %>% 
              NormalizeData() %>% FindVariableFeatures() %>% 
              ScaleData(vars.to.regress=noise) %>% RunPCA(npcs=npcs) %>% RunUMAP(dims=seq(npcs)) %>%
              FindNeighbors(dims=seq(npcs), k.param=k) %>% FindClusters(algorithm=4) 

  pilot_sweep <- paramSweep_v3(object, PCs=seq(npcs), sct=sct, num.cores=cores)
  
  sweep_stats <- summarizeSweep(pilot_sweep, GT=GT)
  pK_pilot <- find.pK(sweep_stats); setDT(pK_pilot)
  
  nExp_poi <- round(expected_doublet_rate * ncol(object))
  pK <- pK_pilot[pK_pilot[, .I[which.max(BCmetric)]], as.numeric(levels(pK))[pK]]
  
  object <- object %>% doubletFinder_v3(PCs=seq(npcs), pN=0.25, 
                                 pK=pK, nExp=nExp_poi, reuse.pANN=FALSE, sct=sct)
  
  setnames(object@meta.data, old=which(colnames(object@meta.data) %like% 'DF'), new='DF_calls')
  setnames(object@meta.data, old=which(colnames(object@meta.data) %like% 'pANN'), new='DF_score')

  return(object)
}
  