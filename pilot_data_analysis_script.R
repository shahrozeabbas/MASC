# load required packages
pkgs <- c('data.table', 'Seurat', 'dplyr')
invisible(sapply(pkgs, require, character.only = TRUE))

# data directory
txg_dir <- 'Seq1_GEXnFBC_wNTG/SCAF1819_Dox/outs/'

# read data using data.table and Seurat
txg_data <- Read10X(data.dir = paste0(txg_dir, 'filtered_feature_bc_matrix_dox'), unique.features = TRUE)
protospacers <- fread(input = paste0(txg_dir, 'crispr_analysis/protospacer_calls_per_cell.csv'))

# initialize Seurat object
data <- CreateSeuratObject(counts = txg_data[[1]])
# data[['gRNA']] <- CreateAssayObject(counts = txg_data[[2]])

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
                    status = fifelse(guide %like% 'Control', 'non-target', 'perturbed'),
                    NT = fifelse(guide %like% 'Control', 'NT', 'perturbed'))]

mUMI_guides[, target := fifelse(target %chin% 'Human', 'control', target)]

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


# add mitochondrial percentage for cells
data[['percentMT']] <- PercentageFeatureSet(data, pattern = '^MT-')

# filter cells
data <- subset(data, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000 & percentMT < 10)

# transform before cell cycle scoring
data <- data %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData()

# determine cell cycle phase scores
data <- CellCycleScoring(data, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)

# regress cell cycle (go get some coffee)
all.genes <- rownames(data)
noise <- c('S.Score', 'G2M.Score')
data <- data %>% ScaleData(vars.to.regress = noise, features = all.genes) 

# perform PCA
data <- data %>% RunPCA(features = VariableFeatures(data))

# cluster the data
data <- data %>% FindNeighbors(dims = 1:20) %>% FindClusters(resolution = 1.2)

#dimensionality reduction for visualization
data <- data %>% RunUMAP(dims = 1:20) %>% RunTSNE(dims = 1:20)


# reset cell identities to match your group of interest like targets, guides, etc
data@active.ident <- as.factor(data$guide_RNA)
avg_exp <- as.data.frame(AverageExpression(data, features = rownames(data), slot = 'scale.data'))
fwrite(avg_exp, paste0('average_cell_expression_per_guide.csv'))

data@active.ident <- as.factor(data$target_RNA)
avg_exp <- as.data.frame(AverageExpression(data, features = rownames(data), slot = 'scale.data'))
fwrite(avg_exp, paste0('average_cell_expression_per_target.csv'))



# being Mixscape analysis

data <- data %>% CalcPerturbSig(
  assay = 'RNA', 
  slot = 'data',
  gd.class = 'target_RNA',
  nt.cell.class = 'control',
  reduction = 'pca', 
  ndims = 40, 
  num.neighbors = 20, 
  new.assay.name = 'PRTB')

DefaultAssay(object = data) <- 'PRTB'

VariableFeatures(object = data) <- VariableFeatures(object = data[['RNA']])
data <- data %>% ScaleData(do.scale = FALSE, do.center = TRUE) %>% RunPCA(reduction.key = 'prtbpca', reduction.name = 'prtbpca')

data <- data %>% RunUMAP(dims = 1:40, reduction = 'prtbpca', reduction.key = 'prtbumap', reduction.name = 'prtbumap')

custom_theme <- theme(
  plot.title = element_text(size=16, hjust = 0.5), 
  legend.key.size = unit(0.7, "cm"), 
  legend.text = element_text(size = 14))

dim_reduce <- 'umap'

q1 <- DimPlot(
  object = data, 
  group.by = 'guide_RNA', 
  reduction = dim_reduce, 
  pt.size = 0.1, cols = "Dark2", label = F, repel = T) +
  scale_color_manual(values = col_vector) +
  ggtitle("RNA Guides") +
  ylab("UMAP 2") +
  xlab("UMAP 1") +
  custom_theme + theme(legend.position = 'none')

q2 <- DimPlot(
  object = data, 
  group.by = 'Phase', 
  reduction = dim_reduce, 
  pt.size = 0.1, label = F, repel = T) +
  ggtitle("Cell Cycle Phase") +
  ylab("UMAP 2") +
  xlab("UMAP 1") + 
  custom_theme

q3 <- DimPlot(
  object = data,
  group.by = 'status_RNA',
  reduction = dim_reduce, 
  split.by = "NT", 
  ncol = 1, 
  pt.size = 0.1, 
  cols = c("grey39","goldenrod3", 'purple')) +
  ggtitle("Perturbation Status") +
  ylab("UMAP 2") +
  xlab("UMAP 1") + 
  custom_theme
  
(q1 / q2 + plot_layout(guides = 'auto') | q3)




Idents(data) <- 'control'
data <- RunMixscape(
  object = data, 
  assay = 'PRTB', 
  slot = "scale.data", 
  labels = 'target_RNA', 
  nt.class.name = 'control', 
  min.de.genes = 5, 
  iter.num = 10, 
  de.assay = 'RNA', 
  verbose = F)


df <- prop.table(table(data$mixscape_class.global, data$status_RNA), 2)

df2 <- reshape2::melt(df)
df2$Var2 <- as.character(df2$Var2)
test <- df2[which(df2$Var1 == "KO"),]
test <- test[order(test$value, decreasing = TRUE),]
new.levels <- test$Var2
df2$Var2 <- factor(df2$Var2, levels = new.levels )
df2$Var1 <- factor(df2$Var1, levels = c("control", "NP", "KO"))
df2$gene <- sapply(as.character(df2$Var2), function(x) strsplit(x, split = "g")[[1]][1])
df2$guide_number <- sapply(as.character(df2$Var2), 
                           function(x) strsplit(x, split = "g")[[1]][2])
df3 <- df2[-c(which(df2$gene == "non-tar")),]


































