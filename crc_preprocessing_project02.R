#!/usr/bin/env Rscript

# ==============================================================================
# GSE294300 - Clean & Sequential Seurat v5 Analysis Pipeline
#
# Inputs:  GSE294300_RAW/crc_obj.rds
# Outputs: Saved in CRC_server_results/ (objects/, plots/, tables/, markers/)
# ==============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

set.seed(1234)
options(
  stringsAsFactors = FALSE,
  Seurat.object.assay.version = "v5"
)
# 
#Load the dataset
sample_dirs <- list.dirs(here("data"), full.names = TRUE)
sample_dirs

sample_dirs <- sample_dirs[-1] #Remove the parent directory

crc.data <- vector("list", length(sample_dirs)) #Create a list
names(crc.data) <- basename(sample_dirs) #Name the list according to the sample IDs
sample_id <- basename(sample_dirs)

for (i in seq_along(sample_dirs)) {
  crc.data[[i]] <- Read10X(data.dir = sample_dirs[[i]])
  crc.data[[i]]$SampleID <- sample_id
}
saveRDS(crc.data, "output/crc_data.rds") #Save the 10X object 

crc_obj <- CreateSeuratObject(crc_object = crc.data, project = "crc", min.cells = 3, min.features = 200)
# 0. DIRECTORY SETUP & PARAMETERS ----------------------------------------------

input_file <- "GSE294300_RAW/crc_obj.rds"
output_dir <- "CRC_server_results"

# Standard paths matching the workspace tree layout
obj_dir    <- file.path(output_dir, "objects")
plot_dir   <- file.path(output_dir, "plots")
table_dir  <- file.path(output_dir, "tables")
marker_dir <- file.path(output_dir, "markers")

dirs <- c(obj_dir, plot_dir, table_dir, marker_dir)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# Analysis parameters
MIN_FEATURES       <- 200
MAX_FEATURES       <- 6000
MAX_PERCENT_MT     <- 20
N_VARIABLE_FEATS   <- 3000
N_PCS_COMPUTED     <- 50
N_PCS_USED         <- 30
CLUSTER_RESOLUTION <- 0.5

# 36 CRC samples: 18 Normal and 18 Tumor
normal_accessions <- sprintf("GSM%d", 8901920:8901937)
tumor_accessions  <- sprintf("GSM%d", 8901938:8901955)

sample_info <- data.frame(
  sample_id = c(normal_accessions, tumor_accessions),
  condition = factor(rep(c("Normal", "Tumor"), each = 18), levels = c("Normal", "Tumor"))
)

write.csv(sample_info, file.path(table_dir, "sample_metadata.csv"), row.names = FALSE)

# 1. REBUILD SAMPLES WITH UNIQUE CELL BARCODES ---------------------------------

message("--> Step 1: Loading raw object and split-rebuilding sample layers...")
source_object <- readRDS(input_file)
count_layers  <- Layers(source_object[["RNA"]], search = "^counts")

sample_list <- vector("list", length(count_layers))

for (i in seq_along(count_layers)) {
  layer     <- count_layers[[i]]
  sample_id <- sub("^counts\\.", "", layer)
  
  # Prefix cell barcodes to avoid duplicated barcode conflicts across layers
  counts <- LayerData(source_object, assay = "RNA", layer = layer)
  colnames(counts) <- paste(sample_id, colnames(counts), sep = "_")
  
  # Instantiate Seurat object
  obj <- CreateSeuratObject(counts = counts, project = sample_id)
  obj$sample_id <- sample_id
  obj$condition <- as.character(sample_info$condition[sample_info$sample_id == sample_id])
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  # 2. SAMPLE-LEVEL QUALITY CONTROL & FILTERING --------------------------------
  
  # QC Violin plot per sample
  p_violin <- VlnPlot(
    obj, 
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
    ncol = 3, 
    pt.size = 0
  ) + plot_annotation(title = paste("QC Metrics -", sample_id))
  
  ggsave(
    filename = file.path(plot_dir, paste0(sample_id, "_QC_violin.png")),
    plot = p_violin, width = 12, height = 4, dpi = 180
  )
  
  # QC Scatter plots
  p_scat1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  p_scat2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  ggsave(
    filename = file.path(plot_dir, paste0(sample_id, "_QC_scatter.png")),
    plot = p_scat1 + p_scat2, width = 10, height = 4.5, dpi = 180
  )
  
  # Subset high-quality cells
  obj <- subset(
    obj,
    subset = nFeature_RNA > MIN_FEATURES &
             nFeature_RNA < MAX_FEATURES &
             percent.mt < MAX_PERCENT_MT
  )
  
  sample_list[[i]] <- obj
}

rm(source_object)
gc()

# Merge all filtered sample objects
message("--> Step 2: Merging filtered samples into single Seurat object...")
crc <- merge(
  x = sample_list[[1]],
  y = sample_list[-1],
  project = "GSE294300_CRC"
)

rm(sample_list)
gc()

# Checkpoint 1: Save Merged QC Object
qc_obj_path <- file.path(obj_dir, "01_QC_filtered_merged.rds")
message("Saving Checkpoint 1: ", qc_obj_path)
saveRDS(crc, qc_obj_path)

# 3. NORMALIZATION, FEATURE SELECTION, SCALING & PCA ---------------------------

message("--> Step 3: Normalizing data and identifying variable features...")
DefaultAssay(crc) <- "RNA"

crc <- NormalizeData(crc, normalization.method = "LogNormalize", scale.factor = 10000)
crc <- FindVariableFeatures(crc, selection.method = "vst", nfeatures = N_VARIABLE_FEATS)

# Export variable features metadata
variable_genes <- VariableFeatures(crc)
write.csv(
  data.frame(gene = variable_genes), 
  file.path(table_dir, "variable_features.csv"), 
  row.names = FALSE
)

# Plot Variable Features
p_var <- VariableFeaturePlot(crc)
p_var_top10 <- LabelPoints(plot = p_var, points = head(variable_genes, 10), repel = TRUE)

ggsave(
  filename = file.path(plot_dir, "variable_features.png"),
  plot = p_var + p_var_top10, width = 12, height = 5, dpi = 200
)

# Scale variable genes only to optimize runtime and memory usage
crc <- ScaleData(crc, features = variable_genes)

# Compute PCA
message("--> Step 4: Computing PCA embeddings...")
crc <- RunPCA(crc, features = variable_genes, npcs = N_PCS_COMPUTED, seed.use = 1234)

# PCA Plots
p_elbow <- ElbowPlot(crc, ndims = N_PCS_COMPUTED)
ggsave(file.path(plot_dir, "PCA_elbow.png"), p_elbow, width = 7, height = 5, dpi = 200)

p_pca_sample <- DimPlot(crc, reduction = "pca", group.by = "sample_id", raster = TRUE)
p_pca_cond   <- DimPlot(crc, reduction = "pca", group.by = "condition", raster = TRUE)
ggsave(
  filename = file.path(plot_dir, "PCA_sample_condition.png"),
  plot = p_pca_sample + p_pca_cond, width = 14, height = 6, dpi = 200
)

# PCA Heatmap
png(file.path(plot_dir, "PCA_DimHeatmap_PC1_15.png"), width = 2400, height = 1800, res = 180)
DimHeatmap(crc, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

# Checkpoint 2: Save PCA & Scaled Object
pca_obj_path <- file.path(obj_dir, "02_normalized_scaled_PCA.rds")
message("Saving Checkpoint 2: ", pca_obj_path)
saveRDS(crc, pca_obj_path)

# 4. NEIGHBOR GRAPH, CLUSTERING & UMAP -----------------------------------------

message("--> Step 5: Building SNN graph, clustering, and computing UMAP...")
crc <- FindNeighbors(crc, reduction = "pca", dims = 1:N_PCS_USED)
crc <- FindClusters(crc, resolution = CLUSTER_RESOLUTION, random.seed = 1234)
crc <- RunUMAP(crc, reduction = "pca", dims = 1:N_PCS_USED, seed.use = 1234)

# UMAP Plots
p_umap_clust <- DimPlot(crc, reduction = "umap", label = TRUE, repel = TRUE, raster = TRUE) + NoLegend()
p_umap_sample <- DimPlot(crc, reduction = "umap", group.by = "sample_id", raster = TRUE)
p_umap_cond   <- DimPlot(crc, reduction = "umap", group.by = "condition", raster = TRUE)

ggsave(
  filename = file.path(plot_dir, "UMAP_cluster_sample_condition.png"),
  plot = p_umap_clust + p_umap_sample + p_umap_cond,
  width = 18, height = 5.5, dpi = 220
)

# Export Cluster Distribution Statistics
cluster_counts <- as.data.frame(table(cluster = Idents(crc)))
write.csv(cluster_counts, file.path(table_dir, "cluster_sizes.csv"), row.names = FALSE)

cluster_by_sample <- as.data.frame(table(cluster = Idents(crc), sample_id = crc$sample_id))
write.csv(cluster_by_sample, file.path(table_dir, "cluster_by_sample.csv"), row.names = FALSE)

# Checkpoint 3: Save Clustered UMAP Object
umap_obj_path <- file.path(obj_dir, "03_clustered_UMAP.rds")
message("Saving Checkpoint 3: ", umap_obj_path)
saveRDS(crc, umap_obj_path)

# 5. MARKERS ANALYSIS & ANNOTATION PREPARATION ---------------------------------

message("--> Step 6: Identifying marker genes...")
crc <- JoinLayers(crc, assay = "RNA")

# Differential expression test
markers <- FindAllMarkers(
  crc,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.25,
  max.cells.per.ident = 1500,
  random.seed = 1234
)

markers <- markers %>% arrange(cluster, p_val_adj, desc(avg_log2FC))
write.csv(markers, file.path(marker_dir, "all_cluster_markers.csv"), row.names = FALSE)

# Filter top 10 markers per cluster
top10_markers <- markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10, with_ties = FALSE) %>%
  ungroup()

write.csv(top10_markers, file.path(marker_dir, "top10_markers_per_cluster.csv"), row.names = FALSE)

# Top 10 Markers Heatmap Plot
heatmap_genes <- intersect(unique(top10_markers$gene), rownames(crc))
if (length(heatmap_genes) > 0) {
  crc <- ScaleData(crc, features = heatmap_genes, verbose = FALSE)
  p_heatmap <- DoHeatmap(crc, features = heatmap_genes, raster = TRUE) + NoLegend()
  
  ggsave(
    filename = file.path(plot_dir, "top10_cluster_markers_heatmap.png"),
    plot = p_heatmap, width = 14, height = 10, dpi = 200
  )
}

# CRC Expression DotPlot (Epithelial, Immune, Stromal markers)
crc_panel <- c("EPCAM", "KRT8", "CD3D", "CD4", "CD8A", "MS4A1", "COL1A1", "PECAM1")
crc_panel <- intersect(crc_panel, rownames(crc))

if (length(crc_panel) > 0) {
  p_dot <- DotPlot(crc, features = crc_panel) + RotatedAxis()
  ggsave(
    filename = file.path(plot_dir, "CRC_annotation_DotPlot.png"),
    plot = p_dot, width = 10, height = 6, dpi = 200
  )
}

# Export Template for Cell Type Annotation
annotation_template <- data.frame(
  cluster = levels(Idents(crc)),
  cell_type = "TO_ANNOTATE",
  notes = ""
)
write.csv(annotation_template, file.path(table_dir, "cluster_annotation_template.csv"), row.names = FALSE)

# Checkpoint 4: Save Final Pre-annotation Object
final_obj_path <- file.path(obj_dir, "04_markers_preannotation.rds")
message("Saving Checkpoint 4: ", final_obj_path)
saveRDS(crc, final_obj_path)

message("--> Pipeline completed successfully!")