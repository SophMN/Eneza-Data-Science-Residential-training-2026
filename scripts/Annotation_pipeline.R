#!/usr/bin/env Rscript


# CELL ANNOTATION PIPELINE: SINGLER + STROMAL MARKERS
 

library(Seurat)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)

set.seed(1234)

# 1. LOAD FULL SEURAT OBJECT 

seurat_obj <- readRDS("crc_results/objects/04_markers_preannotation.rds")

# Ensure Seurat v5 layers are joined
seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")
DefaultAssay(seurat_obj) <- "RNA"

# 2. STEP A: SUBSET REPRESENTATIVE CELLS FOR SPEED 

singler_cells <- unlist(lapply(levels(Idents(seurat_obj)), function(cl) {
  cluster_cells <- which(Idents(seurat_obj) == cl)
  if (length(cluster_cells) > 1500) {
    sample(cluster_cells, 1500)
  } else {
    cluster_cells
  }
}))

seurat_sub <- seurat_obj[, singler_cells]

# Save SingleR input subset checkpoint
dir.create("crc_results/objects", showWarnings = FALSE, recursive = TRUE)
saveRDS(list(expression = GetAssayData(seurat_sub, layer = "counts"), 
             clusters = Idents(seurat_sub)), 
        "crc_results/objects/SingleR_input_1500_cells_per_cluster.rds")

# 3. STEP B: RUN AUTOMATED ANNOTATION (SINGLER) 

sce_sub <- as.SingleCellExperiment(seurat_sub, assay = "RNA")

ref_hpca <- HumanPrimaryCellAtlasData()
ref_blue <- BlueprintEncodeData()
ref_mona <- MonacoImmuneData()

# Prediction at cluster level using the subsetted object and original identities
pred_hpca <- SingleR(test = sce_sub, ref = ref_hpca, labels = ref_hpca$label.main, clusters = Idents(seurat_sub))
pred_blue <- SingleR(test = sce_sub, ref = ref_blue, labels = ref_blue$label.main, clusters = Idents(seurat_sub))
pred_mona <- SingleR(test = sce_sub, ref = ref_mona, labels = ref_mona$label.main, clusters = Idents(seurat_sub))

# Consolidate prediction results into a summary data frame
annotations_summary <- data.frame(
  cluster = rownames(pred_hpca),
  HPCA = pred_hpca$labels,
  BlueprintEncode = pred_blue$labels,
  MonacoImmune = pred_mona$labels
)

# 4. MAP ANNOTATIONS BACK TO FULL SEURAT OBJECT --------------------------------

# Preserve original cluster identifiers
seurat_obj$seurat_clusters_orig <- Idents(seurat_obj)

# Create lookup vectors (cluster_id -> annotation)
hpca_map <- setNames(pred_hpca$labels, rownames(pred_hpca))
blue_map <- setNames(pred_blue$labels, rownames(pred_blue))
mona_map <- setNames(pred_mona$labels, rownames(pred_mona))

# Assign metadata per cell based on its cluster ID (unname() avoids Seurat barcode mismatch)
seurat_obj$SingleR_HPCA      <- unname(hpca_map[as.character(Idents(seurat_obj))])
seurat_obj$SingleR_Blueprint <- unname(blue_map[as.character(Idents(seurat_obj))])
seurat_obj$SingleR_Monaco    <- unname(mona_map[as.character(Idents(seurat_obj))])

# Set default identity to HPCA annotations
Idents(seurat_obj) <- "SingleR_HPCA"

# 5. EVALUATE CAF AND GLIAL STROMAL MARKERS 

caf_markers <- list(
  General_CAF         = c("DCN", "COL1A1", "COL1A2", "LUM"),
  myCAF               = c("ACTA2", "TAGLN", "MYH11"),
  iCAF                = c("IL6", "CXCL12", "CFD", "PDGFRA"),
  vCAF                = c("RGS5", "MCAM"),
  proinflammatory_CAF = c("INHBA", "SULF1", "CTHRC1")
)

glial_markers <- c("SOX10", "PLP1", "MPZ", "S100B", "MBP", "CDH19", "PTPRZ1")

# Filter panels to retain only genes present in the dataset
valid_caf_genes   <- intersect(unique(unlist(caf_markers)), rownames(seurat_obj))
valid_glial_genes <- intersect(glial_markers, rownames(seurat_obj))

# Compute per-cell module scores
seurat_obj <- AddModuleScore(seurat_obj, features = list(valid_caf_genes), name = "CAF_Score_")
seurat_obj <- AddModuleScore(seurat_obj, features = list(valid_glial_genes), name = "Glial_Score_")

# 5bis. MANUALLY ANNOTATE CLUSTERS & EXPORT TEMPLATE 

annotation_template <- data.frame(
  cluster = levels(seurat_obj$seurat_clusters_orig),
  SingleR_HPCA = hpca_map[levels(seurat_obj$seurat_clusters_orig)],
  manual_cell_type = "TO_ANNOTATE",
  notes = ""
)

dir.create("crc_results/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(
  annotation_template, 
  "crc_results/tables/cluster_annotation_template.csv", 
  row.names = FALSE
)

# 6. EXPORT SUMMARY TABLES 

write.csv(
  annotations_summary, 
  "crc_results/tables/SingleR_cluster_annotations.csv", 
  row.names = FALSE
)

# Compute average expression per cluster on original clusters
cluster_stromal_expression <- AverageExpression(
  seurat_obj,
  features = c(valid_caf_genes, valid_glial_genes),
  group.by = "seurat_clusters_orig",
  layer = "data"
)$RNA

write.csv(
  cluster_stromal_expression,
  "crc_results/tables/cluster_stromal_markers_average.csv"
)

# 7. VISUALIZATION AND PLOTTING 

# A. UMAP overlay of primary SingleR annotations
p_umap <- DimPlot(seurat_obj, reduction = "umap", group.by = "SingleR_HPCA", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP - SingleR Annotations (HPCA)")

# B. FeaturePlot overlays for CAF and Glial module scores
p_scores <- FeaturePlot(
  seurat_obj,
  features = c("CAF_Score_1", "Glial_Score_1"),
  reduction = "umap",
  cols = c("grey", "red")
)

# C. DotPlot of detailed marker gene expression across original clusters
p_dot_stromal <- DotPlot(
  seurat_obj,
  features = list(CAF = valid_caf_genes, Glial = valid_glial_genes),
  group.by = "seurat_clusters_orig"
) + RotatedAxis() + ggtitle("CAF & Glial Marker Expression across Clusters")

# Save plots to disk
dir.create("crc_results/plots", showWarnings = FALSE, recursive = TRUE)
ggsave("crc_results/plots/UMAP_annotated_SingleR.png", p_umap, width = 10, height = 7)
ggsave("crc_results/plots/UMAP_CAF_and_Glial_scores.png", p_scores, width = 12, height = 5)
ggsave("crc_results/plots/DotPlot_CAF_Glial_markers.png", p_dot_stromal, width = 14, height = 6)

# 8. SAVE ANNOTATED SEURAT OBJECT 

saveRDS(seurat_obj, "crc_results/objects/04_seurat_annotated.rds")