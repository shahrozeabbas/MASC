pkgs <- c('data.table', 'Seurat', 'WGCNA', 'pheatmap')
invisible(sapply(pkgs, require, character.only = TRUE))


data <- fread('/PerturbSeqDemo/sequencing/exp1_5/exp1_5_MasterCellTable_Adamson2020_EnsemblIDs.csv')

gene_ids <- names(data)[21:ncol(data)]

filter_cols <- c('cell_barcode', gene_ids)
raw_counts <- data[, ..filter_cols]

barcodes <- raw_counts[, cell_barcode]
raw_counts[, cell_barcode := NULL]
setDF(raw_counts, rownames = barcodes)


my_data <- CreateSeuratObject(counts = raw_counts, min.cells = 200, project = 'exp1_5')
