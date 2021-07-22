#set wd and load packages
setwd("/data/ShernData/SWISNF_CrispriData2/02_PrimaryAnalysisOutput/00_FullCellrangerOutputs/")
getwd()

library(mlbench)
library(visdat)
library(naniar)
library(gtsummary)
library(labelled)
library(ggplot2)
library(limma)
library(stringr)
library(dplyr)
library(gtools)
library(tidyr)
library(tidyverse)
library(ggpubr)
library(kableExtra)
library(rio)
library(readxl)
library(fs)
library(ggpubr)
library(lubridate)
library(base)
library(Seurat)
library(SeuratDisk)
library(data.table)
library(clusterProfiler)
library(enrichplot)
library(pathview)
library(BiocGenerics)
library(readr)
library(stringr)


#set wd to SWISNF sample outs
#FIRST LANE
Pilot.data1 <- Read10X_h5("SCAF2095_SNF1/outs/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
#Pilot.data2 <- Read10X_h5("SNF2_outs/SCAF2096_SNF2_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
#Pilot.data3 <- Read10X_h5("SNF3_outs/SCAF2097_SNF3_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)


#The data is in a list form containing two separate matrices, one for gene expression (RNA data) and one for CRISPR guide identifier. Called these matrices individually below
exp = Pilot.data1[["Gene Expression"]]
exp2 = Pilot.data1[["CRISPR Guide Capture"]]

#Creating a Seurat object from the 10x data matrices
Pilot.data1 <- CreateSeuratObject(counts = exp, project = "Gene_Expression"  , min.cells = 0, min.features = 0)

#Add the CRISPR identifying matrix as a second assay object
Pilot.data1[['Guides']] <- CreateAssayObject(counts = exp2[,colnames(x=Pilot.data1)])
#Check the data so that it has been imported correctly
Pilot.data1


#read in guides for this lane
#guides <- fread('SCAF2095_SNF1_protospacer_calls_per_cell.csv')
guides <- fread('SCAF2095_SNF1/outs/crispr_analysis/protospacer_calls_per_cell.csv')
#remove cells with multiple barcode calls
guides <- guides[num_features == 1]
#add gene name column, identified as 'target', by splitting the guide identifying column (which is a 'number_gene' format)
guides <- guides[, target := tstrsplit(feature_call, '_')[1]]
##Subset cells with identified guides <- this is the data which will be progressed through the analysis
Pilot.data1 = subset(Pilot.data1, cells = guides$cell_barcode)


# store guide labels for single guides in a vector
g <- guides[, feature_call]

# check data type
class(g) # character

# name list by cell barcode
names(g) <- guides[, cell_barcode]

# change data type to factor for assignment purposes 
g <- as.factor(g) 

# check data type 
class(g) #factor

# add meta data to Seurat object (gyuide RNAs)
Pilot.data1 <- Pilot.data1 %>% AddMetaData(metadata = g, col.name = 'guide_RNA')

# Repeat for target meta data
t <- guides[, target]
names(t) <- guides[, cell_barcode]
t <- as.factor(t)
Pilot.data1 <- Pilot.data1 %>% AddMetaData(metadata = as.factor(t), col.name = 'genes')

##Add NT meta data
NT <- guides[, NT := target]
NT <- guides[NT %chin% 'NEG', NT := 'NT']
NT <-  NT[, NT]
NT <- as.factor(NT)
names(NT) <- guides[, cell_barcode]
Pilot.data1 <- Pilot.data1 %>% AddMetaData(metadata = NT, col.name = 'NT')

##Add Status meta data
status <- guides[target %chin% 'NEG', status := 'non-targeting']
status  <- guides[!target %chin% 'NEG', status := 'perturbed']
status  <-  guides[, status]
status  <- as.factor(status)
names(status) <- guides[, cell_barcode]
Pilot.data1 <- Pilot.data1 %>% AddMetaData(metadata = status, col.name = 'Status')


#SECOND LANE
#This is single cell sequencing data where each cell contains a CRISPR guide
#Pilot.data1 <- Read10X_h5("SNF1_outs/SCAF2095_SNF1_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
Pilot.data2 <- Read10X_h5("SCAF2096_SNF2/outs/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
#Pilot.data3 <- Read10X_h5("SNF3_outs/SCAF2097_SNF3_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)


#The data is in a list form containing two separate matrices, one for gene expression (RNA data) and one for CRISPR guide identifier. Called these matrices individually below
exp = Pilot.data2[["Gene Expression"]]
exp2 = Pilot.data2[["CRISPR Guide Capture"]]

#Creating a Seurat object from the 10x data matrices
Pilot.data2 <- CreateSeuratObject(counts = exp, project = "Gene_Expression"  , min.cells = 0, min.features = 0)

#Add the CRISPR identifying matrix as a second assay object
Pilot.data2[['Guides']] <- CreateAssayObject(counts = exp2[,colnames(x=Pilot.data2)])
#Check the data so that it has been imported correctly
Pilot.data2


#read in guides for this lane
guides <- fread('SCAF2096_SNF2/outs/crispr_analysis/protospacer_calls_per_cell.csv')
#remove cells with multiple barcode calls
guides <- guides[num_features == 1]
#add gene name column, identified as 'target', by splitting the guide identifying column (which is a 'number_gene' format)
guides <- guides[, target := tstrsplit(feature_call, '_')[1]]
##Subset cells with identified guides <- this is the data which will be progressed through the analysis
Pilot.data2 = subset(Pilot.data2, cells = guides$cell_barcode)


# store guide labels for single guides in a vector
g <- guides[, feature_call]

# check data type
class(g) # character

# name list by cell barcode
names(g) <- guides[, cell_barcode]

# change data type to factor for assignment purposes 
g <- as.factor(g) 

# check data type 
class(g) #factor

# add meta data to Seurat object (gyuide RNAs)
Pilot.data2 <- Pilot.data2 %>% AddMetaData(metadata = g, col.name = 'guide_RNA')

# Repeat for target meta data
t <- guides[, target]
names(t) <- guides[, cell_barcode]
t <- as.factor(t)
Pilot.data2 <- Pilot.data2 %>% AddMetaData(metadata = as.factor(t), col.name = 'genes')

##Add NT meta data
NT <- guides[, NT := target]
NT <- guides[NT %chin% 'NEG', NT := 'NT']
NT <-  NT[, NT]
NT <- as.factor(NT)
names(NT) <- guides[, cell_barcode]
Pilot.data2 <- Pilot.data2 %>% AddMetaData(metadata = NT, col.name = 'NT')

##Add Status meta data
status <- guides[target %chin% 'NEG', status := 'non-targeting']
status  <- guides[!target %chin% 'NEG', status := 'perturbed']
status  <-  guides[, status]
status  <- as.factor(status)
names(status) <- guides[, cell_barcode]
Pilot.data2 <- Pilot.data2 %>% AddMetaData(metadata = status, col.name = 'Status')


#THIRD LANE
#This is single cell sequencing data where each cell contains a CRISPR guide
#Pilot.data1 <- Read10X_h5("SNF1_outs/SCAF2095_SNF1_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
#Pilot.data2 <- Read10X_h5("SNF2_outs/SCAF2096_SNF2_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
Pilot.data3 <- Read10X_h5("SCAF2097_SNF3/outs/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)


#The data is in a list form containing two separate matrices, one for gene expression (RNA data) and one for CRISPR guide identifier. Called these matrices individually below
exp = Pilot.data3[["Gene Expression"]]
exp2 = Pilot.data3[["CRISPR Guide Capture"]]

#Creating a Seurat object from the 10x data matrices
Pilot.data3 <- CreateSeuratObject(counts = exp, project = "Gene_Expression"  , min.cells = 0, min.features = 0)

#Add the CRISPR identifying matrix as a second assay object
Pilot.data3[['Guides']] <- CreateAssayObject(counts = exp2[,colnames(x=Pilot.data3)])
#Check the data so that it has been imported correctly
Pilot.data3

#read in guides for this lane
guides <- fread('SCAF2097_SNF3/outs/crispr_analysis/protospacer_calls_per_cell.csv')
#remove cells with multiple barcode calls
guides <- guides[num_features == 1]
#add gene name column, identified as 'target', by splitting the guide identifying column (which is a 'number_gene' format)
guides <- guides[, target := tstrsplit(feature_call, '_')[1]]
##Subset cells with identified guides <- this is the data which will be progressed through the analysis
Pilot.data3 = subset(Pilot.data3, cells = guides$cell_barcode)


# store guide labels for single guides in a vector
g <- guides[, feature_call]

# check data type
class(g) # character

# name list by cell barcode
names(g) <- guides[, cell_barcode]

# change data type to factor for assignment purposes 
g <- as.factor(g) 

# check data type 
class(g) #factor

# add meta data to Seurat object (gyuide RNAs)
Pilot.data3 <- Pilot.data3 %>% AddMetaData(metadata = g, col.name = 'guide_RNA')

# Repeat for target meta data
t <- guides[, target]
names(t) <- guides[, cell_barcode]
t <- as.factor(t)
Pilot.data3 <- Pilot.data3 %>% AddMetaData(metadata = as.factor(t), col.name = 'genes')

##Add NT meta data
NT <- guides[, NT := target]
NT <- guides[NT %chin% 'NEG', NT := 'NT']
NT <-  NT[, NT]
NT <- as.factor(NT)
names(NT) <- guides[, cell_barcode]
Pilot.data3 <- Pilot.data3 %>% AddMetaData(metadata = NT, col.name = 'NT')

##Add Status meta data
status <- guides[target %chin% 'NEG', status := 'non-targeting']
status  <- guides[!target %chin% 'NEG', status := 'perturbed']
status  <-  guides[, status]
status  <- as.factor(status)
names(status) <- guides[, cell_barcode]
Pilot.data3 <- Pilot.data3 %>% AddMetaData(metadata = status, col.name = 'Status')



#merge these datasets
Pilot.data <- merge(Pilot.data1, y = c(Pilot.data2, Pilot.data3), add.cell.ids = c( "1", "2", "3"), project = "SWISNF")



# Setup custom theme for plotting.
custom_theme <- theme(
  plot.title = element_text(size=16, hjust = 0.5), 
  legend.key.size = unit(0.7, "cm"), 
  legend.text = element_text(size = 14))




#Once I have added all the required information as metadata, I need to do some standard QC

#Add another metadata column indicating the proportion of transcripts mapping to mitochondrial DNA
Pilot.data[["percent.mt"]] <- PercentageFeatureSet(Pilot.data, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Pilot.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Use QC metrics to subset higher quality cells from dataset
Pilot.data <- subset(Pilot.data, subset = nFeature_RNA > 1500 & nFeature_RNA < 7500 & percent.mt < 12)

#Re-visualise metrics after subset
VlnPlot(Pilot.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(Pilot.data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Pilot.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Find Variable Features
DefaultAssay(object = Pilot.data) <- 'RNA'

#Pilot.data <- NormalizeData(object = Pilot.data) %>% FindVariableFeatures() %>% ScaleData()
Pilot.data <- NormalizeData(object = Pilot.data, normalization.method = "LogNormalize", scale.factor = 10000) 
Pilot.data <- FindVariableFeatures(object = Pilot.data, selection.method = "vst", nfeatures = 5000) 

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Pilot.data), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Pilot.data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 
plot2

#Scaling the data
Pilot.data <- NormalizeData(object = Pilot.data) %>% FindVariableFeatures() %>% ScaleData()
Pilot.data <- RunPCA(object = Pilot.data)


# Run Uniform Manifold Approximation and Projection (UMAP) to visualize clustering in 2-D.
Pilot.data <- RunUMAP(object = Pilot.data, dims = 1:40)
DimPlot(Pilot.data, label = TRUE) + NoLegend()



##Cell Cycle Calling 
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  
# We can segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
cc.genes.updated.2019[!cc.genes.updated.2019 %in% rownames(Pilot.data)]
head(Pilot.data[[]])

Pilot.data <- CellCycleScoring(Pilot.data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(Pilot.data[[]])

Pilot.data <- RunPCA(Pilot.data, features = c(s.genes, g2m.genes))
DimPlot(Pilot.data)


##MIXSCAPE VIGNETTE

##Plot guides
p1 <- DimPlot(
  object = Pilot.data, 
  group.by = 'guide_RNA') +
  custom_theme
p1

##Plot phase
p2 <- DimPlot(
  object = Pilot.data, 
  group.by = 'Phase')+
  custom_theme
p2

p3 <- DimPlot(
  object = Pilot.data, 
  group.by = 'Status')+
  custom_theme
p3



##Cleaned up to this point!


##SAveRDS
#saveRDS(Pilot.data, file = "SWISNF_merged.data.rds")

#Pilot.data <- readRDS("Pilot.data.rds")

##Run perturbation signature calculation function 
Pilot.data<- CalcPerturbSig(
  object = Pilot.data, 
  assay = "RNA", 
  slot = "data", 
  gd.class ="NT", 
  nt.cell.class = "NT",
  reduction = "pca",
  ndims = 40, 
  num.neighbors = 20,
  split.by = NULL, 
  new.assay.name = "PRTB")

# Prepare PRTB assay for dimensionality reduction: 
# Normalize data, find variable features and center data.
DefaultAssay(object = Pilot.data) <- 'PRTB'

# Use variable features from RNA assay.
VariableFeatures(object = Pilot.data) <- VariableFeatures(object = Pilot.data[["RNA"]])
Pilot.data <- ScaleData(object = Pilot.data, do.scale = F, do.center = T)

# Run PCA to reduce the dimensionality of the data.
Pilot.data <- RunPCA(object = Pilot.data, reduction.key = 'prtbpca', reduction.name = 'prtbpca')

# Run UMAP to visualize clustering in 2-D.
Pilot.data <- RunUMAP(
  object = Pilot.data, 
  dims = 1:40, 
  reduction = 'prtbpca', 
  reduction.key = 'prtbumap', 
  reduction.name = 'prtbumap')

##save as RDS object
#saveRDS(Pilot.data, file = "Calc.Perturb.Sig.SWISNF.data.rds")
#Pilot.data <- readRDS("Calc.Perturb.Sig.Pilot.data.rds")

##Plot guides
q1 <- DimPlot(
  object = Pilot.data, 
  group.by = 'guide_RNA') +
  custom_theme
q1

##Plot phase
q2 <- DimPlot(
  object = Pilot.data, 
  group.by = 'Phase') +
  custom_theme
q2

##Plot phase
q3 <- DimPlot(
  object = Pilot.data, 
  group.by = 'Status')+
  custom_theme
q3

# Run mixscape.
Pilot.data <- RunMixscape(
  object = Pilot.data, 
  assay = "PRTB", 
  slot = "scale.data", 
  labels = "NT", 
  nt.class.name = "NT",
  new.class.name = "mixscape_class",
  min.de.genes = 2, 
  logfc.threshold = 0.1,
  iter.num = 10, 
  de.assay = "RNA",
  verbose = F
)


# Run mixscape using fine mode.
##Add Status meta data
#guides <- as.factor(Pilot.data@meta.data$guide_RNA)
#levels(guides) <- gsub('NEG_CTRL-221', 'NT', levels(guides))
#levels(guides) <- gsub('NEG_CTRL-222', 'NT', levels(guides)) 
#levels(guides) <- gsub('NEG_CTRL-223', 'NT', levels(guides))
#levels(guides) <- gsub('NEG_CTRL-224', 'NT', levels(guides))
#levels(guides) <- gsub('NEG_CTRL-225', 'NT', levels(guides))
#Pilot.data <- Pilot.data %>% AddMetaData(metadata = guides, col.name = 'guide_ID')

#Pilot.data <- RunMixscape(
#  object = Pilot.data, 
#  assay = "PRTB", 
#  slot = "scale.data", 
#  labels = "guide_ID", 
#  nt.class.name = "NT",
#  new.class.name = "mixscape_class",
#  min.de.genes = 2, 
#  logfc.threshold = 0.1,
#  iter.num = 10, 
# de.assay = "RNA",
#  verbose = F,
#  fine.mode = TRUE
#)

# Calculate percentage of KO cells for all target gene classes.
#If you want to do this guide by guide, you need to add metadata where all NT controls grouped as NT
#But rest of cells identified by individual guide labels
guides <- as.factor(Pilot.data@meta.data$guide_RNA)
levels(guides) <- gsub('NEG_CTRL-221', 'NT', levels(guides))
levels(guides) <- gsub('NEG_CTRL-222', 'NT', levels(guides)) 
levels(guides) <- gsub('NEG_CTRL-223', 'NT', levels(guides))
levels(guides) <- gsub('NEG_CTRL-224', 'NT', levels(guides))
levels(guides) <- gsub('NEG_CTRL-225', 'NT', levels(guides))
Pilot.data <- Pilot.data %>% AddMetaData(metadata = guides, col.name = 'guide_ID')

Idents(Pilot.data) <- 'guide_RNA'
df <- prop.table(table(Pilot.data$mixscape_class.global, Pilot.data$guide_ID),2)


df2 <- reshape2::melt(df)
df2$Var2 <- as.character(df2$Var2)
test <- df2[which(df2$Var1 == "KO"),]
test <- test[order(test$value, decreasing = T),]
new.levels <- test$Var2
df2$Var2 <- factor(df2$Var2, levels = new.levels )
df2$Var1 <- factor(df2$Var1, levels = c("NT", "NP", "KO"))
df2$gene <- sapply(as.character(df2$Var2), function(x) strsplit(x, split = "_")[[1]][1])
df2$guide_number <- sapply(as.character(df2$Var2), 
                           function(x) strsplit(x, split = "_")[[1]][2])
df3 <- df2[-c(which(df2$gene == "NT")),]

p1 <- ggplot(df3, aes(x = guide_number, y = value*100, fill= Var1)) +
  geom_bar(stat= "identity") +
  theme_classic()+
  scale_fill_manual(values = c("grey49", "grey79","coral1")) + 
  ylab("% of cells") +
  xlab("sgRNA")

p1 + theme(axis.text.x = element_text(size = 18, hjust = 1), 
           axis.text.y = element_text(size = 18), 
           axis.title = element_text(size = 16), 
           strip.text = element_text(size=16, face = "bold")) + 
  facet_wrap(vars(gene),ncol = 5, scales = "free") +
  labs(fill = "mixscape class") +theme(legend.title = element_text(size = 14),
                                       legend.text = element_text(size = 12))


# Remove non-perturbed cells and run LDA to reduce the dimensionality of the data.
Idents(Pilot.data) <- "mixscape_class.global"
sub <- subset(Pilot.data, idents = c("KO", "NT"))
sub 
#saveRDS(sub , file = "Subset_Perturbed_Cells.data.rds")
#sub <- readRDS("Subset_Perturbed_Cells.data.rds")
sub <- MixscapeLDA(object = sub, assay = "RNA", pc.assay = "PRTB", labels = "NT", nt.label = "NT", 
                   npcs = 10, logfc.threshold = 0.25, verbose = F)


# Use LDA results to run UMAP and visualize cells on 2-D.
#changed dims to 1:8 based on dims of sub seurat object

sub <- RunUMAP(sub, dims = 1:10,reduction = "lda", reduction.key = "ldaumap", reduction.name = "ldaumap")

# Visualize UMAP clustering results.
Idents(sub) <- "mixscape_class"
sub$mixscape_class <- as.factor(sub$mixscape_class)
p <- DimPlot(sub, reduction = "ldaumap", label = T, repel = T, label.size = 5)
p

DimPlot(sub, reduction = "ldaumap", split.by = "mixscape_class", ncol = 8, label = TRUE)

ob.list <- SplitObject(sub, split.by = "mixscape_class")
plot.list <- lapply(X = ob.list, FUN = function(x) {
  DimPlot(x, reduction = "umap", label = FALSE, label.size = 4)
})

#sub summary characteristics
Idents(sub) <- "mixscape_class"
All_sub <- count(sub@meta.data, vars = mixscape_class)
write.csv(All_sub, "~/Desktop/SWISNF_output_plots/sub_counts")


#Grab the differentially expressed genes identified between each cluster
Idents(object = sub) <- "NT"

sub@meta.data$NT <- as.factor((sub@meta.data$NT))
gene <- levels(sub$NT)
genes <- gene[gene != "NT"]

out <- vector("list")
for (i in genes){ 
  out[[i]]<- 
    FindMarkers(sub, ident.1 = 'NT', ident.2 = i, group.by = "NT")
  write.csv(out[[i]], file = sprintf("~/Desktop/SWISNF_output_plots/All_%s.csv",i))}
#all saved to DE_genes folder



##Filter datasets for only significant hits
out_reduc <- map(out, 
                 filter, p_val_adj < 0.05)

#write each of these files into an individual csv
iwalk(out_reduc, ~write.csv(.x, str_c("~/Desktop/SWISNF_output_plots/signif_", .y, ".csv"), row.names=TRUE))

#combine the significant hits to one file
out_select <- bind_rows(out_reduc, .id = "sample_id")
out_select <- out_select %>%
  rownames_to_column(var = 'gene')

write.csv(out_select, file = "~/Desktop/SWISNF_output_plots/combined_signif_hits.csv")

##