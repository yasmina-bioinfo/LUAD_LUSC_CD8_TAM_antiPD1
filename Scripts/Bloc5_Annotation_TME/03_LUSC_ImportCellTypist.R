#!/usr/bin/env Rscript
# ============================================================
# 03_LUSC_ImportCellTypist.R
# GSE243013 — Bloc5_Annotation_TME, Script 03 (LUSC)
#
# Purpose: import the CellTypist predictions produced in Python (Script 02),
# re-attach them to the Seurat object by cell barcode, and tabulate the
# majority_voting label distribution within EVERY cluster (no cluster
# skipped based on how clean its top10 markers look).
#
# Input:  Objects/01_LUSC_seu_clustered.rds
#         Objects/celltypist_export/LUSC_celltypist_predictions.csv
# Output: Objects/03_LUSC_seu_celltypist.rds
#         Results/Bloc5_Annotation_TME/Tables/03_LUSC_celltypist_per_cluster.csv
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(dplyr)
  library(data.table)
})

IN_OBJ  <- here("Objects", "01_LUSC_seu_clustered.rds")
IN_PRED <- here("Objects", "celltypist_export", "LUSC_celltypist_predictions.csv")
OUT_TAB <- here("Results", "Bloc5_Annotation_TME", "Tables")
dir.create(OUT_TAB, showWarnings = FALSE, recursive = TRUE)

message("Loading clustered Seurat object...")
seu <- readRDS(IN_OBJ)
seu <- JoinLayers(seu)
message("Cells: ", ncol(seu))

# -----------------------------
# 1) Read CellTypist predictions and match to Seurat object by barcode
# -----------------------------
message("Loading CellTypist predictions...")
pred <- fread(IN_PRED)
message("Predictions loaded: ", nrow(pred), " cells")

# Sanity check: same cell count, same barcodes (order may differ)
stopifnot(
  "Cell count mismatch between Seurat object and predictions" = 
    nrow(pred) == ncol(seu),
  "Barcode mismatch — check export/import alignment" = 
    all(sort(pred$V1) == sort(colnames(seu)))
)
message("Barcode alignment check passed.")

# Match order explicitly by barcode (never assume row order is preserved)
pred_matched <- pred[match(colnames(seu), pred$V1), ]

seu$CellTypist_predicted <- pred_matched$predicted_labels
seu$CellTypist_majority  <- pred_matched$majority_voting

message("CellTypist labels attached to Seurat object.")

# -----------------------------
# 2) Systematic per-cluster distribution check — ALL clusters,
# no exception, regardless of how clean the top10 markers looked
# -----------------------------
message("Checking CellTypist majority_voting distribution per cluster...")

celltypist_summary <- seu@meta.data %>%
  count(seurat_clusters, CellTypist_majority, name = "Freq") %>%
  arrange(seurat_clusters, desc(Freq))

fwrite(celltypist_summary, file.path(OUT_TAB, "03_LUSC_celltypist_per_cluster.csv"))
message("Saved: 03_LUSC_celltypist_per_cluster.csv")

# Quick view: top label(s) per cluster
for (cl in levels(seu$seurat_clusters)) {
  cat("\n=== Cluster", cl, "===\n")
  tab <- celltypist_summary %>% filter(seurat_clusters == cl) %>% head(8)
  print(tab)
}

# -----------------------------
# 3) Save object with CellTypist labels attached
# -----------------------------
saveRDS(seu, here("Objects", "03_LUSC_seu_celltypist.rds"))
message("Saved: Objects/03_LUSC_seu_celltypist.rds")
message("DONE Script 03 — LUSC CellTypist import")