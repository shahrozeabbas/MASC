##Cleaned Code
##Load Packages
library(devtools)
library(Seurat)
library(dplyr)
library(Matrix)
library(base)
library(methods)
library(utils)
library(stats)
library(gdata)
library(graphics)
library(grDevices)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(scales)
library(reshape2)
#library(SeuratDisk)

# Setup custom theme for plotting.
custom_theme <- theme(
  plot.title = element_text(size=16, hjust = 0.5), 
  legend.key.size = unit(0.7, "cm"), 
  legend.text = element_text(size = 14))

##Read in data
setwd('/Users/abbass2/Desktop/ShernLab/PerturbSeq/PilotData_Seurat')
Pilot.data <- Read10X("dox/outs/filtered_feature_bc_matrix")
exp = Pilot.data[["Gene Expression"]]
exp2 = Pilot.data[["CRISPR Guide Capture"]]

##Read in protospacer_calls from sequencing core facility
##Guide_Name <- read.csv("crispr_analysis/protospacer_calls_per_cell.csv")


Pilot.data <- CreateSeuratObject(counts = exp, project = "Gene_Expression"  , min.cells = 0, min.features = 0)
##Pilot.data <-NormalizeData(Pilot.data)
Pilot.data[['Guides']] <- CreateAssayObject(counts = exp2[,colnames(x=Pilot.data)])
Pilot.data


library(data.table)
protospacers <- fread('ctl/metadata/protospacer_calls_per_cell.csv')
guides <- protospacers[,  .(guide = unlist(strsplit(feature_call, '|', fixed = TRUE)), nUMI = unlist(strsplit(num_umis, '|', fixed = TRUE))), keyby = cell_barcode]
#write.csv(guides,'barcode_guide_nUMI_target.csv')

# read barcodes, guide labels, number of UMIs, and guide targets
# guides <- fread('dox/metadata/barcode_guide_nUMI_target.csv')
guides <- guides[, target := tstrsplit(guide, '_')[1]]
# collect guides with the most UMI counts
max_umi_guides <- guides[guides[, .I[which.max(nUMI)], by = cell_barcode]$V1]

# store guide labels for max umi counts in a vector
g <- max_umi_guides[, guide]

# check data type
class(g) # character

# name list by cell barcode
names(g) <- max_umi_guides[, cell_barcode]

# change data type to factor for assignment purposes 
g <- as.factor(g) 

# check data type 
class(g) #factor

# add meta data to Seurat object (gyuide RNAs)
Pilot.data <- Pilot.data %>% AddMetaData(metadata = g, col.name = 'guide_RNA')

# Repeat for target meta data
t <- max_umi_guides[, target]
names(t) <- max_umi_guides[, cell_barcode]
t <- as.factor(t)
Pilot.data <- Pilot.data %>% AddMetaData(metadata = as.factor(t), col.name = 'genes')

##Subset cells with identified guides <- this is the data which will be progressed through the analysis
Pilot.data = subset(Pilot.data, cells = max_umi_guides$cell_barcode)

##Add NT meta data
NT <- max_umi_guides[, NT := target]
NT <- max_umi_guides[target %chin% 'Human', NT := 'NT']
NT <-  NT[, NT]
NT <- as.factor(NT)
names(NT) <- max_umi_guides[, cell_barcode]
Pilot.data <- Pilot.data %>% AddMetaData(metadata = NT, col.name = 'NT')

##Add Status meta data
guides <- max_umi_guides[target %chin% 'Human', status := 'non-targeting']
guides <- max_umi_guides[!target %chin% 'Human', status := 'perturbed']
guides <-  max_umi_guides[, status]
guides <- as.factor(guides)
names(guides) <- max_umi_guides[, cell_barcode]
Pilot.data <- Pilot.data %>% AddMetaData(metadata = guides, col.name = 'Status')



##Standard Pre-processing Workflow (Sajita lab Vignette: Guided Clustering Tutorial)
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Pilot.data[["percent.mt"]] <- PercentageFeatureSet(Pilot.data, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(Pilot.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Use QC metrics to subset higher quality cells from dataset
Pilot.data <- subset(Pilot.data, subset = nFeature_RNA > 1000 & nFeature_RNA < 3000 & percent.mt < 6)

#Re-visualise metrics after subset
VlnPlot(Pilot.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Pilot.data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Pilot.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

##Normalise and scale the cleaned data
Pilot.data <- NormalizeData(object = Pilot.data) %>% FindVariableFeatures() %>% ScaleData()


#Visualise
Pilot.data <- RunPCA(object = Pilot.data)
Pilot.data <- RunUMAP(object = Pilot.data, dims = 1:40)
DimPlot(Pilot.data)


##Cell Cycle Calling 
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cc.genes[!cc.genes %in% rownames(Pilot.data)]

# Assign Cell Cycle Scores
Pilot.data <- CellCycleScoring(Pilot.data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(Pilot.data[[]])

##MIXSCAPE VIGNETTE
# Prepare RNA assay for dimensionality reduction: 
# Normalize data, find variable features and scale data.
DefaultAssay(object = Pilot.data) <- 'RNA'
Pilot.data <- NormalizeData(object = Pilot.data) %>% FindVariableFeatures() %>% ScaleData()

# Run Principle Component Analysis (PCA) to reduce the dimensionality of the data.
Pilot.data <- RunPCA(object = Pilot.data)

# Run Uniform Manifold Approximation and Projection (UMAP) to visualize clustering in 2-D.
Pilot.data <- RunUMAP(object = Pilot.data, dims = 1:40)

DimPlot(Pilot.data, label = TRUE) + NoLegend()

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


##SAveRDS
#saveRDS(Pilot.data, file = "Pilot.data.rds")

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
  num.neighbors = 30,
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

Pilot.data <- RunTSNE(
  object = Pilot.data, 
  dims = 1:40, 
  reduction = 'prtbpca', 
  reduction.key = 'prtbtsne', 
  reduction.name = 'prtbtsne')

##save as RDS object
#saveRDS(Pilot.data, file = "Calc.Perturb.Sig.Pilot.data.rds")

##Plot guides
q1 <- DimPlot(
  object = Pilot.data, 
  group.by = 'guide_RNA', reduction = 'prtbtsne') +
  custom_theme
q1

##Plot phase
q2 <- DimPlot(
  object = Pilot.data, 
  group.by = 'Phase') +
  custom_theme
q2


# Run mixscape.
Idents(Pilot.data) <- 'NT'
Pilot.data <- RunMixscape(
  object = Pilot.data, 
  assay = "PRTB", 
  slot = "scale.data", 
  labels = "NT", 
  nt.class.name = "NT",
  new.class.name = "mixscape_class",
  min.de.genes = 1, 
  logfc.threshold = 0,
  iter.num = 10, 
  de.assay = "RNA", 
  verbose = F,
  fine.mode = TRUE,
  fine.mode.labels = "guide_RNA")


# Remove non-perturbed cells and run LDA to reduce the dimensionality of the data.
Idents(Pilot.data) <- "mixscape_class.global"
sub <- subset(Pilot.data, idents = c("KO", "NT"))

sub <- MixscapeLDA(object = sub, assay = "RNA", pc.assay = "PRTB", labels = "NT", nt.label = "NT", 
                   npcs = 10, logfc.threshold = 0, verbose = F)


# Use LDA results to run UMAP and visualize cells on 2-D.
#changed dims to 1:8 based on dims of sub seurat object

sub <- RunUMAP(sub, dims = 1:10,reduction = "lda", reduction.key = "ldaumap", reduction.name = "ldaumap")

# Visualize UMAP clustering results.
Idents(sub) <- "mixscape_class"
sub$mixscape_class <- as.factor(sub$mixscape_class)
p <- DimPlot(sub, reduction = "ldaumap", label = T, repel = T, label.size = 5)

col = setNames(object = hue_pal()(12), nm = levels(sub$mixscape_class))
names(col) <- c(names(col)[1:7], "NT", names(col)[9:12])
col[8] <- "grey39"

p + scale_color_manual(values = col, drop = FALSE) + ylab("UMAP 2") + xlab("UMAP 1") + custom_theme

##end of mixscape v3.1

#in mixscape v3.2
# Calculate percentage of KO cells for all target gene classes.
df <- prop.table(table(Pilot.data$mixscape_class.global, Pilot.data$NT),2)

df2 <- reshape2::melt(df)
df2$Var2 <- as.character(df2$Var2)
test <- df2[which(df2$Var1 == "KO"),]
test <- test[order(test$value, decreasing = T),]
new.levels <- test$Var2
df2$Var2 <- factor(df2$Var2, levels = new.levels )
df2$Var1 <- factor(df2$Var1, levels = c("NT", "NP", "KO"))
df2$gene <- sapply(as.character(df2$Var2), function(x) strsplit(x, split = "g")[[1]][1])
df2$guide_number <- sapply(as.character(df2$Var2), 
                           function(x) strsplit(x, split = "g")[[1]][2])
df3 <- df2[-c(which(df2$gene == "NT")),]

p1 <- ggplot(df3, aes(x = guide_number, y = value*100, fill= Var1)) +
  geom_bar(stat= "identity") +
  theme_classic()+
  scale_fill_manual(values = c("red", "grey79","blue")) + 
  ylab("% of cells") +
  xlab("sgRNA")

p1 + theme(axis.text.x = element_text(size = 18, hjust = 1), 
           axis.text.y = element_text(size = 18), 
           axis.title = element_text(size = 16), 
           strip.text = element_text(size=16, face = "bold")) + 
  facet_wrap(vars(gene),ncol = 5, scales = "free") +
  labs(fill = "mixscape class") +theme(legend.title = element_text(size = 14),
                                       legend.text = element_text(size = 12))


# Explore the perturbation scores of cells.
PlotPerturbScore(object = Pilot.data, 
                 target.gene.ident = "TMSB10",
                 group.by = "mixscape_class", 
                 col = "coral2") +labs(fill = "mixscape class")

# Inspect the posterior probability values in NP and KO cells.
VlnPlot(Pilot.data, "mixscape_class_p_ko", idents = c("NT", "TMSB10 KO", "TMSB10 NP")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) + 
  NoLegend() +
  ggtitle("mixscape posterior probabilities")



# Run DE analysis and visualize results on a heatmap ordering cells by their posterior 
# probability values.
Idents(object = Pilot.data) <- "NT"
MixscapeHeatmap(object = Pilot.data, 
                ident.1 = "NT", 
                ident.2 = "ACTL6A", 
                balanced = F, 
                assay = "RNA", 
                logfc.threshold = 0,
                max.genes = 20, angle = 0, 
                group.by = "mixscape_class", 
                max.cells.group = 300, 
                size=3.5) + NoLegend()

