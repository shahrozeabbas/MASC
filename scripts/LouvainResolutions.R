


rs <- seq(0.5, 2, 0.1)

nClusters <- function(object = NULL, resolutions = NULL) {

  n <- vector('list', length = length(resolutions))
  
  for (r in seq_along(resolutions)) {
    
    object <- FindClusters(object = object, resolution = resolutions[r])
    n[r] <- length(levels(object$seurat_clusters))
    
  }
  
  n <- unlist(n)

  return(n)
}



pcs <- seq(5, 50)
pc_list <- vector('list', length = length(pcs))
for (pc in 1:length(pcs)) pc_list[pc] <- js[, sum(PC[1:pcs[pc]] > pcs[pc])]
max(unlist(pc_list))



