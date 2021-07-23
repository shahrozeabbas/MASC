args <- commandArgs(trailingOnly=TRUE)
require(ggplot2); require(dplyr)

composition <- data.table::fread(args[1]); group <- args[2]; root <- args[3]



blank_theme <- theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    axis.text.x=element_blank(),
    plot.title=element_text(size=14, face='bold')
  )


n_colors <- nrow(composition[, .N, keyby = group])
pal <- colorRampPalette(wesanderson::wes_palette('Darjeeling1', type = 'discrete')[c(1, 3, 5)])

p <- composition %>% 
  ggplot(aes(x='', y=percent, fill=group)) +
    
    blank_theme +
    geom_col(width=1) +
    coord_polar('y', start=0) +
    facet_wrap(vars(get(group))) +
    scale_fill_manual(values=pal(n_colors)) + 
    ggtitle(label = paste(toupper(group), 'COMPOSITION FOR CLUSTERS'))

ggsave(plot=p, filename=paste0(root, '_', group, '_composition_pies.pdf'))
