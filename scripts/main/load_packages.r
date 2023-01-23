PKGS <- c('data.table', 'Seurat', 'dplyr', 'ggplot2', 'future', 'reticulate', 'parallel', 'splitstackshape', 'SoupX', 'mixtools', 'rootSolve')

invisible(sapply(PKGS, require, character.only=TRUE))
invisible(lapply(list.files('scripts/utility', full.names=TRUE, pattern='\\.r$'), source))


data_path <- 'data/'

setDTthreads(threads=1)

