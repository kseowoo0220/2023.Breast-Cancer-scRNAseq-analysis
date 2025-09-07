library(readr)
library(Matrix)
library(tidyverse)
library(scCustomize)

# load data
barcode <- read_tsv("Wu_etal_2021_BRCA_scRNASeq/count_matrix_barcodes.tsv", col_names = FALSE)
feature <- read_tsv("Wu_etal_2021_BRCA_scRNASeq/count_matrix_genes.tsv", col_names = FALSE)
counts <- readMM("Wu_etal_2021_BRCA_scRNASeq/count_matrix_sparse.mtx")

rownames(counts) <- feature$X1
colnames(counts) <- barcode$X1

project <- CreateSeuratObject(counts = counts)
rm(counts)

project <- saveRDS(project, file = "project.rds")

# pre-processing
library(Seurat)
metadata <- read.csv(
  file = "metadata.csv",
  header = TRUE,
  row.names = 1
)

project <- AddMetaData(object = project, metadata = metadata)

project[["percent.mt"]] <- PercentageFeatureSet(project, pattern = "^MT-")

VlnPlot(project, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

project <- subset(project, subset = nFeature_RNA > 250 & nCount_RNA > 500 & percent.mt < 20)

plot1 <- FeatureScatter(project, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(project, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

project <- NormalizeData(project, normalization.method = "LogNormalize", scale.factor = 10000)

project <- NormalizeData(project)

project <- FindVariableFeatures(project, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(project), 10)

project <- ScaleData(project)

project <- RunPCA(project, features = VariableFeatures(object = project), npcs = 100)


# run UMAP
set.seed(1234)
ElbowPlot(project)

project <- FindNeighbors(project, dims = 1:16)
project <- FindClusters(project, resolution = 0.6)

project <- RunUMAP(project, dims = 1:16)

DimPlot(project, reduction = "umap", label = TRUE)

project <- saveRDS(project, file = "project.rds")

# annotation
project <- readRDS("project.rds")
temp_colour_pal_2 <- data.frame(celltype = c("T-cells", "B-cells", "Plasmablasts", "Myeloid", "Epithelial",
                                             "Cycling", "Mesenchymal", "Endothelial"))
project@meta.data$celltype_major <- factor(project@meta.data$celltype_major,
                                           levels = unique(temp_colour_pal_2$celltype))
DimPlot(project, reduction = "umap", label.size = 3, pt.size = .01, group.by = "celltype_major",
        label = TRUE)


# marker genes for annotation
# 0, 1, 2, 22, 16
# T cell
cluster0.markers <- FindMarkers(project, only.pos = TRUE, ident.1 = 0, min.pct = 0.25,
                                logfc.threshold = 0.25)
# T cell
cluster1.markers <- FindMarkers(project, only.pos = TRUE, ident.1 = 1, min.pct = 0.25,
                                logfc.threshold = 0.25)
# T cell
cluster6.markers <- FindMarkers(project, only.pos = TRUE, ident.1 = 6, min.pct = 0.25,
                                logfc.threshold = 0.25)
# T cell
cluster30.markers <- FindMarkers(project, only.pos = TRUE, ident.1 = 30, min.pct = 0.25,
                                 logfc.threshold = 0.25)
# Epithelial
cluster20.markers <- FindMarkers(project, only.pos = TRUE, ident.1 = 20, min.pct = 0.25,
                                 logfc.threshold = 0.25)
# Basal_Myoepithelial
cluster21.markers <- FindMarkers(project, only.pos = TRUE, ident.1 = 21, min.pct = 0.25,
                                 logfc.threshold = 0.25)
# Epithelial
cluster22.markers <- FindMarkers(project, only.pos = TRUE, ident.1 = 22, min.pct = 0.25,
                                 logfc.threshold = 0.25)
# Cycling
cluster23.markers <- FindMarkers(project, only.pos = TRUE, ident.1 = 23, min.pct = 0.25,
                                 logfc.threshold = 0.25)
# cycling
cluster26.markers <- FindMarkers(project, only.pos = TRUE, ident.1 = 26, min.pct = 0.25,
                                 logfc.threshold = 0.25)
# Epithelial
cluster28.markers <- FindMarkers(project, only.pos = TRUE, ident.1 = 28, min.pct = 0.25,
                                 logfc.threshold = 0.25)
# Epithelial 
cluster29.markers <- FindMarkers(project, only.pos = TRUE, ident.1 = 29, min.pct = 0.25,
                                 logfc.threshold = 0.25)
# CAFs
cluster31.markers <- FindMarkers(project, only.pos = TRUE, ident.1 = 31, min.pct = 0.25,
                                 logfc.threshold = 0.25)
# Cycling
cluster32.markers <- FindMarkers(project, only.pos = TRUE, ident.1 = 32, min.pct = 0.25,
                                 logfc.threshold = 0.25)
# Epithelial
cluster33.markers <- FindMarkers(project, only.pos = TRUE, ident.1 = 33, min.pct = 0.25,
                                 logfc.threshold = 0.25)
# Epithelial
cluster5.markers <- FindMarkers(project, only.pos = TRUE, ident.1 = 5, min.pct = 0.25,
                                logfc.threshold = 0.25)
# Epithelial
cluster7.markers <- FindMarkers(project, only.pos = TRUE, ident.1 = 7, min.pct = 0.25,
                                logfc.threshold = 0.25)
# Epithelial
cluster20.markers <- FindMarkers(project, only.pos = TRUE, ident.1 = 20, min.pct = 0.25,
                                 logfc.threshold = 0.25)
# Epithelial
cluster13.markers <- FindMarkers(project, only.pos = TRUE, ident.1 = 13, min.pct = 0.25,
                                 logfc.threshold = 0.25)
# Myeloid
cluster27.markers <- FindMarkers(project, only.pos = TRUE, ident.1 = 27, min.pct = 0.25,
                                 logfc.threshold = 0.25)
# CAFs
cluster4.markers <- FindMarkers(project, only.pos = TRUE, ident.1 = 4, min.pct = 0.25,
                                logfc.threshold = 0.25)
# Epithelial
cluster8.markers <- FindMarkers(project, only.pos = TRUE, ident.1 = 8, min.pct = 0.25,
                                logfc.threshold = 0.25)
# Endothelial
cluster25.markers <- FindMarkers(project, only.pos = TRUE, ident.1 = 25, min.pct = 0.25,
                                 logfc.threshold = 0.25)
# Epithelial
cluster17.markers <- FindMarkers(project, only.pos = TRUE, ident.1 = 17, min.pct = 0.25,
                                 logfc.threshold = 0.25)
# CAFs
cluster10.markers <- FindMarkers(project, only.pos = TRUE, ident.1 = 10, min.pct = 0.25,
                                 logfc.threshold = 0.25)
# CAFs
cluster11.markers <- FindMarkers(project, only.pos = TRUE, ident.1 = 11, min.pct = 0.25,
                                 logfc.threshold = 0.25)
# Epithelial 
cluster15.markers <- FindMarkers(project, only.pos = TRUE, ident.1 = 15, min.pct = 0.25,
                                 logfc.threshold = 0.25)
# Epithelial 
cluster14.markers <- FindMarkers(project, only.pos = TRUE, ident.1 = 14, min.pct = 0.25,
                                 logfc.threshold = 0.25)
# Epithelial
cluster18.markers <- FindMarkers(project, only.pos = TRUE, ident.1 = 18, min.pct = 0.25,
                                 logfc.threshold = 0.25)
# Cycling 
cluster19.markers <- FindMarkers(project, only.pos = TRUE, ident.1 = 19, min.pct = 0.25,
                                 logfc.threshold = 0.25)

# manually annoated the clusteres
new.cluster.ids <- c("T-cells", "T-cells", "Endothelial", "Myeloid", "Mesenchymal",
                     "Epithelial", "T-cells", "Epithelial", "Epithelial", "B-cells",
                     "Mesenchymal", "Mesenchymal", "Plasmablasts", "Epithelial",
                     "Epithelial", "Epithelial", "Myeloid", "Epithelial", "Epithelial",
                     "Cycling", "Epithelial", "Epithelial", "Epithelial", "Cycling", 
                     "Plasmablasts", "Endothelial", "Cycling", "Myeloid", "Epithelial", 
                     "Epithelial", "T-cells", "Mesenchymal", "Cycling", "Epithelial")
names(new.cluster.ids) <- levels(project)
project <- RenameIdents(project, new.cluster.ids)
project@meta.data[["celltype_major"]] <- Idents(project)
DimPlot(project, label = TRUE, pt.size = 0.5) + NoLegend()


FeaturePlot(project, features = c("EPCAM", "MKI67", "CD3D", "CD68", "MS4A1", "JCHAIN", "PECAM1", "PDGFRB"))

# visualization
temp_colour_pal_2_stromal <- temp_colour_pal_2[!temp_colour_pal_2$celltype %in% c("Epithelial"),,drop=F]
temp_df <- project@meta.data[,colnames(project@meta.data) %in% c("nFeature_RNA", "nCount_RNA", "subtype", "orig.ident", "celltype_major_fig"),drop=F]
temp_df$barcode <- rownames(temp_df)
rownames(temp_df) <- NULL
temp_df$subtype <- gsub("\\+","",temp_df$subtype)
temp_df$subtype <- factor(temp_df$subtype,levels=c("TNBC", "HER2", "ER"))

project.list <- SplitObject(project, split.by = "celltype_major")
DimPlot(project.list$`Epithelial`, reduction = "umap", group.by = "orig.ident")
DimPlot(project.list$`Epithelial`, reduction = "umap", group.by = "subtype")

project.list <- saveRDS(project.list, file = "projectlist.rds")


# co-expression analysis
library(fcoex)
library(Seurat)
library(ggplot2)
library(SingleCellExperiment)
# Set up metadata as desired for aggregation and DE analysis
project <-SetIdent(project, value = project@meta.data$celltype_major)
metadata$cluster_id <- factor(project@active.ident)

# Extract raw counts and metadata to create SingleCellExperiment object
counts <- project@assays$RNA@counts 
View(counts)
metadata <- project@meta.data
View(metadata)

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), colData = metadata)

target <- colData(sce)
target <- target$celltype_major

exprs <- as.data.frame(assay(sce, 'counts'))


project.list <- readRDS("projectlist.rds")

# Extract a Seurat list of mesenchymal cells from the list and convert it to a Seurat object
mesenchymal <- project.list["Mesenchymal"]
mesenchymal <- mesenchymal[["Mesenchymal"]]

# Create a new fcoex object
exprs <- data.frame(GetAssayData(mesenchymal))
target <- Idents(mesenchymal)
fc <- new_fcoex(data.frame(exprs), target)

rm(project.list)

fc <- discretize(fc)
fc <- find_cbf_modules(fc,n_genes = 70, verbose = FALSE, is_parallel = FALSE)


library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
theme_set(theme_cowplot())
set.seed(12345)
allowWGCNAThreads(nThreads = 8)

project <- SetupForWGCNA(project,
                         gene_select = "fraction",
                         fraction = 0.05,
                         wgcna_name = "hdwgcna")

project <- MetacellsByGroups(seurat_obj = project,
                             group.by = c("celltype_major", "orig.ident"),
                             reduction = 'pca',
                             k = 25,
                             max_shared = 10,
                             ident.group = "celltype_major")

project <- NormalizeMetacells(project)

project <- SetDatExpr(
  project,
  group_name = "Mesenchymal",
  group.by = "celltype_major",
  assay = "RNA",
  slot = 'data'
)

# project <- SetDatExpr(  project, group_name = "T-cells",  group.by = "celltype_major",  assay = "RNA",  slot = 'data.immune')

project <- TestSoftPowers(
  project,
  networkType = 'signed'
)

plot_list <- PlotSoftPowers(project)

wrap_plots(plot_list, ncol = 2)

project <- ConstructNetwork(
  project, soft_power=6,
  setDatExpr=FALSE,
  tom_name = 'Mesenchymal',
  overwrite_tom = TRUE
)

tom <- read("TOM/Mesenchymal_TOM.rda")
PlotDendrogram(project, main='Mesenchymal hdWGCNA Dendrogram')

coexpression <- ScaleData(coexpression, features=VariableFeatures(coexpression))
coexpression <- ModuleEigengenes(
  coexpression,
  group.by.vars="orig.ident"
)

coexpression <- saveRDS(coexpression, file = "coexpression.rds")

coexpression <- readRDS("coexpression.rds")

hMEs <- GetMEs(coexpression)
MEs <- GetMEs(coexpression, harmonized = FALSE)



coexpression <- ModuleConnectivity(
  coexpression,
  group.by = 'celltype_major', group_name = 'Mesenchymal'
)


p <- PlotKMEs(coexpression,
              ncol = 5)

# get the module assignment table:
modules <- GetModules(coexpression)

# show the first 6 columns:
head(modules[,1:6])

# get hub genes
hub_df <- GetHubGenes(coexpression, n_hubs = 10)

head(hub_df)

# with Seurat method
coexpression <- ModuleExprScore(
  coexpression,
  n_genes = 25,
  method='Seurat'
)

# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
# library(UCell)
# coexpression <- ModuleExprScore(coexpression,n_genes = 25,method='UCell')

ModuleCorrelogram(coexpression)



# gene enrichment packages
# if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

# BiocManager::install("GeneOverlap")

library(enrichR)
library(GeneOverlap)

# dbs<-c('MSigDB_Oncogenic_Signatures','GO_Molecular_Function_2021')
#enrichr databases to test
# dbs <- c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021')
# dbs <- "MSigDB_Oncogenic_Signatures"

# dbs<-'WikiPathways_2019_Human'

dbs<-c('MSigDB_Hallmark_2020')
# perform enrichment tests
coexpression <- RunEnrichr(
  coexpression,
  dbs= dbs, # character vector of enrichr databases to test
  max_genes = 100 # number of genes per module to test
)

# retrieve the output table
# enrich_df <- GetEnrichrTable(coexpression)
# enrich_df1 <- GetEnrichrTable(coexpression)
# enrich_df2 <- GetEnrichrTable(coexpression)
# enrich_df3 <- GetEnrichrTable(coexpression)
# enrich_df4 <- GetEnrichrTable(coexpression)
enrich_df5 <- GetEnrichrTable(coexpression)

# enrichr dotplot
EnrichrDotPlot(
  coexpression,
  mods = "all", # use all modules (this is the default behavior)
  database = "MSigDB_Hallmark_2020", # this has to be one of the lists we used above!!!
  n_terms=1 # number of terms for each module
)

MEs <- GetMEs(coexpression, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
coexpression@meta.data <- cbind(coexpression@meta.data, MEs)

p1 <- DotPlot(coexpression, features=mods, group.by = 'celltype_major')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p1 <- p1 +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

# plot output
p1
# 



EnrichrDotPlot(
  coexpression,
  mods = "all", # use all modules (this is the default behavior)
  database = "MSigDB_Hallmark_2020", # this has to be one of the lists we used above!!!
  n_terms=1 # number of terms for each module
)

rename_list <- list(
  "salmon" = 'Oxidative Phosphorylation1', "turquoise" = 'Angiogenesis',"brown" = 'Epithelial Mesenchymal Transition1',"yellow" = 'Epithelial Mesenchymal Transition2', "blue" = 'Allograft Rejection', "black" = 'TNF-alpha signaling(NF-kB)1',"red" = 'TNF-alpha signaling(NF-kB)2',"magenta"='TNF-alpha signaling(NF-kB)3',"green"='Myc Targets V1_1',"tan"='Myc Targets V1_2',"lightcyan"='TGF-beta Signalling',"greenyellow"='DNA Repair',"midnightblue"='Estrogen Response Early',"pink"='Oxidative Phosphorylation2',"cyan"="UV Response Up","grey60"='Oxidative Phosphorylation3',"purple"='E2F Targets',"lightgreen"='Myc Targets V1_3')

coexpression <- ResetModuleNames(
  coexpression,
  new_name = rename_list
)

# print out the new module names
modules <- GetModules(coexpression)
print(levels(modules$module))

# enrichr dotplot

MEs <- GetMEs(coexpression, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
coexpression@meta.data <- cbind(coexpression@meta.data, MEs)

p1 <- DotPlot(coexpression, features = mods, group.by = 'celltype_major')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p1 <- p1 +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

# plot output
p1

