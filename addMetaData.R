# load required packages
pkgs <- c('data.table', 'Seurat', 'pheatmap')
invisible(sapply(require, pkgs, character.only = TRUE))

# data directory
txg_dir <- 'Seq1_GEXnFBC_wNTG/SCAF1819_Dox/outs/'

# read data using data.table and Seurat
txg_data <- Read10X(data.dir = paste0(txg_dir, 'filtered_feature_bc_matrix'), unique.features = TRUE)
protospacers <- fread(input = paste0(txg_dir, 'crispr_analysis/protospacer_calls_per_cell.csv'))

# initialize Seurat object
data <- CreateSeuratObject(counts = txg_data[[1]])
data[['gRNA']] <- CreateAssayObject(counts = txg_data[[2]])

# parse protospacer file 
guides <- protospacers[, .(guide = unlist(strsplit(feature_call, '|', fixed = TRUE)), 
                 nUMI = unlist(strsplit(num_umis, '|', fixed = TRUE))), keyby = cell_barcode]

# read barcodes, guide labels, number of UMIs, and guide targets
# guides <- fread('~/Desktop/barcode_guide_nUMI_target_status_NT.csv')
# collect guides with the most UMI counts
mUMI_guides <- guides[guides[, .I[which.max(nUMI)], by = cell_barcode]$V1]

# store barcode names and metadata column names for Seurat object
barcodes <- mUMI_guides[, cell_barcode]
metadata_names <- c('guide_RNA', 'target_RNA', 'status_RNA', 'NT')
#add columns to guides table for Mixscape Analysis
mUMI_guides[, `:=` (target = tstrsplit(guide, '_')[[1]], 
               status = fifelse(guide %like% 'control', 'non-target', 'perturbed'),
               NT = fifelse(guide %like% 'control', 'NT', 'perturbed'))]

# Add meta data to Seurat object
md2add <- names(mUMI_guides)
for (m in seq_along(metadata_names)) {
  
  md <- mUMI_guides[, (md2add[m])]; names(md) <- barcodes
  data <- data %>% AddMetaData(metadata = as.factor(md), col.name = metadata_names[m])
  
}

g <- mUMI_guides[, guide]
names(g) <- barcodes
data <- data %>% AddMetaData(metadata = as.factor(g), col.name = 'guide_RNA')

t <- mUMI_guides[, target]
names(t) <- barcodes
data <- data %>% AddMetaData(metadata = as.factor(t), col.name = 'target_RNA')

s <- mUMI_guides[, status]
names(s) <- barcodes
data <- data %>% AddMetaData(metadata = as.factor(s), col.name = 'status_RNA')

nt <- mUMI_guides[, NT]
names(nt) <- barcodes
data <- data %>% AddMetaData(metadata = as.factor(nt), col.name = 'NT')






