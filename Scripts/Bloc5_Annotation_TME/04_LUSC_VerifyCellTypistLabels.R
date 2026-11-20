#!/usr/bin/env Rscript
# ============================================================
# 04_LUSC_VerifyCellTypistLabels.R
# GSE243013 — Bloc5_Annotation_TME, Script 04 (LUSC)
#
# Purpose: verify the CellTypist majority_voting labels within each of the
# 20 clusters. Every label subgroup is checked individually against its
# canonical markers, regardless of its proportion within the parent cluster.
# This avoids the subjectivity of choosing a grouping level (e.g. CD4 vs CD8,
# Treg vs Th) for a percentage-based threshold.
#
# Method: for each cluster, cells are grouped by CellTypist_majority label
# (labels with fewer than MIN_CELLS cells are excluded as noise), then
# FindAllMarkers is run between the label subgroups.
#
# Input:  Objects/03_LUSC_seu_celltypist.rds
# Output: Results/Bloc5_Annotation_TME/Tables/04_LUSC_cluster<N>_markers_all.tsv
#         Results/Bloc5_Annotation_TME/Tables/04_LUSC_cluster<N>_markers_top10.tsv
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(dplyr)
  library(readr)
})

IN_OBJ  <- here("Objects", "03_LUSC_seu_celltypist.rds")
OUT_TAB <- here("Results", "Bloc5_Annotation_TME", "Tables")
dir.create(OUT_TAB, showWarnings = FALSE, recursive = TRUE)

message("Loading annotated Seurat object...")
seu <- readRDS(IN_OBJ)
message("Cells: ", ncol(seu), " | Clusters: ", length(levels(seu$seurat_clusters)))

MIN_CELLS <- 30

# All 20 clusters
TARGET_CLUSTERS <- as.character(0:19)

for (cl in TARGET_CLUSTERS) {
  
  message("=== Resolving cluster ", cl, " ===")
  
  cl_cells <- colnames(seu)[seu$seurat_clusters == cl]
  label_counts <- table(seu$CellTypist_majority[cl_cells])
  keep_labels <- names(label_counts[label_counts >= MIN_CELLS])
  
  message("Cluster ", cl, ": ", length(cl_cells), " cells, ",
          length(keep_labels), " labels retained (>= ", MIN_CELLS, " cells)")
  print(sort(label_counts[keep_labels], decreasing = TRUE))
  
  target_cells <- cl_cells[seu$CellTypist_majority[cl_cells] %in% keep_labels]
  seu_sub <- subset(seu, cells = target_cells)
  Idents(seu_sub) <- seu_sub$CellTypist_majority
  
  message("Running FindAllMarkers on cluster ", cl, " subgroups...")
  markers_sub <- FindAllMarkers(seu_sub, only.pos = TRUE,
                                min.pct = 0.25, logfc.threshold = 0.25)
  message("Total markers found: ", nrow(markers_sub))
  
  if (nrow(markers_sub) == 0) {
    message("No markers found for cluster ", cl, " — skipping top10 extraction.")
    next
  }
  top10_sub <- markers_sub %>%
    group_by(cluster) %>%
    slice_head(n = 10) %>%
    select(cluster, gene, avg_log2FC) %>%
    ungroup()
  
  write_tsv(markers_sub, file.path(OUT_TAB,
                                   paste0("04_LUSC_cluster", cl, "_markers_all.tsv")))
  write_tsv(top10_sub, file.path(OUT_TAB,
                                 paste0("04_LUSC_cluster", cl, "_markers_top10.tsv")))
  
  message("Saved: cluster ", cl, " markers tables")
}

message("DONE Script 04 — LUSC CellTypist label verification")