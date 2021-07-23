plots_outputs <- function(object, out_dir=getwd(), w=14, h=10, group='seurat_clusters', plot=FALSE, output=TRUE, markers=TRUE) {
  
  # plot <- FALSE; output <- FALSE
  require(data.table); require(Seurat); require(ggplot2); require(dplyr); require(SeuratDisk)

  to_disk <- '/gpfs/gsfs12/users/ShernData/'
  root <- paste0(out_dir, 'final_', Project(object))

  if (markers) {
    system(paste('sbatch', paste0('--job-name=', Project(object), '_markers'), paste0(to_disk, 'markers.sh'), root, group))
  }
  
  
  if (plot) output <- TRUE

  if (output) {

    md <- copy(object@meta.data)
    setDT(md, keep.rownames=TRUE)[, orig.ident := NULL]
    setnames(md, old=c('rn', 'nCount_RNA', 'nFeature_RNA'), 
            new=c('cells', 'nUMI', 'nGene'))
    
    umap <- copy(as.data.frame(object[['prtbumap']]@cell.embeddings))
    setDT(umap, keep.rownames='cells')
    
    md <- md[umap, on='cells']; rm(umap)
    fwrite(md, paste0(root, '_metadata.csv'))


    avg <- as.data.frame(object %>% AverageExpression(group.by=group, assay='RNA'))
    setDT(avg, keep.rownames='genes')

    fwrite(x=avg, file=paste0(root, '_cluster_average_expression.csv'))

    genes <- avg[, genes]
    avg[, genes := NULL]
    setDF(avg, rownames=genes)
    
    adj <- WGCNA::adjacency(datExpr=avg, type='signed', power=1)
    hm_colors <- colorRampPalette(wesanderson::wes_palette('Zissou1', type='discrete'))
    pheatmap::pheatmap(mat=adj, show_rownames=TRUE, show_colnames=TRUE, 
            filename=paste0(root, '_targets_average_correlation.pdf'), 
            color=hm_colors(100))


    for (var in c('mixscape_class.global', 'guide_RNA', 'n_guides', 'predicted_doublets')) {
        
        o <- paste0(root, '_', group, '_', var, '_composition.csv')
        
        composition <- md[, {
                            totwt=.N
                            .SD[, .(percent=.N / totwt), by=var]
                                }, by=group][md[, .(nCells=.N), keyby=c(group, var)], 
                          on=c(group, var)]
                    
        fwrite(x=composition, file=o) 
        system(paste('Rscript /gpfs/gsfs12/users/abbass2/scripts_mixscape/composition_pies.r', o, var, root))
    }
    
    
    # md <- fread(paste0(root, '_metadata.csv'))
    if (plot) {

      g <- ggplot(data=md, aes(x=prtbumap_1, y=prtbumap_2)) + theme_classic()

      pal <- wesanderson::wes_palette('Zissou1', type='discrete')[c(1, 3, 5)]
      continuous_stats <- c('nUMI', 'nGene', 'percent.mt', 'doublet_scores', 'mixscape_class_p_ko', 'mcd_score')
      
      invisible(lapply(X=continuous_stats, FUN=function(md_col) {
        
        middle <- md[, mean(c(max(get(md_col)), min(get(md_col))))]
        ggsave(plot=g +
                geom_point(aes(color=get(md_col)), alpha=0.7, size=0.5) + 
                scale_color_gradient2(low=pal[1], mid=pal[2], high=pal[3], midpoint=middle) +
                labs(color=md_col),
              filename=paste0(out_dir, '/', md_col, '_umap.pdf'), width=w, height=h)
        
      }))
        
      
      pal <- colorRampPalette(wesanderson::wes_palette('Moonrise3', type='discrete'))
      discrete_stats <- c('mixscape_class.global', 'predicted_doublets', 'targets', 'guide_RNA', 'mixscape_class', 'n_guides')
      
      invisible(lapply(X=discrete_stats, FUN=function(md_col) {

        n_colors <- nrow(md[, .N, keyby=get(md_col)])  
        ggsave(plot=g + 
                geom_point(aes(color=get(md_col)), size=0.5) + 
                scale_color_manual(values=pal(n_colors)) +
                labs(color=md_col) + 
                guides(color=guide_legend(override.aes=list(size=5))), 
              filename=paste0(out_dir, '/', md_col, '_umap.pdf'), width=w, height=h)
        
      }))
    }
  }
  
}


