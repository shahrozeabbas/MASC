working_dir <- '/data/ShernData/SWISNF/'; setwd(working_dir)

source('scripts/main/load_packages.r')

options(ggrepel.max.overlaps=Inf)

object <- readRDS(snakemake@input[['seurat_object']])


Idents(object) <- 'guide'
df <- prop.table(table(object$mixscape_class.global, object$guide), 2)

df2 <- reshape2::melt(df)
df2$Var2 <- as.character(df2$Var2)
test <- df2[which(df2$Var1 == 'KO'),]
test <- test[order(test$value, decreasing=TRUE),]
new.levels <- test$Var2
df2$Var2 <- factor(df2$Var2, levels=new.levels)
df2$Var1 <- factor(df2$Var1, levels=c('CONTROL', 'NP', 'KO'))
df2$gene <- sapply(as.character(df2$Var2), function(x) strsplit(x, split='_')[[1]][1])
df2$guide_number <- sapply(as.character(df2$Var2), function(x) strsplit(x, split='_')[[1]][2])
df3 <- df2[-c(which(df2$gene == 'NEG')),]
df3 <- df3[-c(which(df3$gene == 'Human')),]

p1 <- ggplot(df3, aes(x=guide_number, y=value*100, fill=Var1)) +
  geom_bar(stat='identity') +
  theme_classic()+
  scale_fill_manual(values=c('grey49', 'grey79','coral1')) +
  ylab('% of cells') + xlab('sgRNA')

p1 <- p1 + theme(axis.text.x=element_text(size=8, hjust=0.5),
           axis.text.y=element_text(size=8),
           axis.title=element_text(size=16),
           strip.text=element_text(size=16, face='bold')) +
            facet_wrap(vars(gene),ncol=5, scales='free') +
            labs(fill='mixscape class') + theme(legend.title=element_text(size=14), legend.text=element_text(size=12))

ggsave(plot=p1, height=22, width=15, filename=snakemake@output[['ko_barplot']])


sub <- readRDS(snakemake@input[['mixscape_sub_object']])

g <- sub %>% DimPlot(reduction='ldaumap', label=TRUE, repel=TRUE, label.size=5, group.by='mixscape_class', pt.size=0.7) + 
      NoLegend() + labs(x='LDAUMAP 1', y='LDAUMAP 2', title='SWISNF Mixscape')  + 
      theme(axis.text=element_text(size=18), axis.title=element_text(size=20), plot.title=element_text(size=30))

ggsave(plot=g, height=9, width=15, filename=snakemake@output[['mixscape_umap']])


markers <- fread(snakemake@input[['markers']])

genes <- markers[target != 'CONTROL', .N, keyby=target][, target]

for (gene in genes) {

    g <- markers[target %chin% gene] %>% ggplot(aes(x=avg_log2FC, y=-log10(p_val), col=grade)) + 
            geom_point(alpha=0.5) + theme_bw() + 
    
            scale_x_continuous(breaks=scales::pretty_breaks(10)) +
            scale_y_continuous(breaks=scales::pretty_breaks(10)) + 
    
            theme(legend.title.align=0.5) + 
            ggtitle(paste('CONTROL vs', gene, 'Markers'))

    ggsave(plot=g, height=8, width=12, filename=paste0('plots/volcano/', gene, '_volcano_plot.pdf'))

}
