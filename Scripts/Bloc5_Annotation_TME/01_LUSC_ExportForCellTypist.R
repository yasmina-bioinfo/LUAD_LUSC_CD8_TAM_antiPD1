#!/usr/bin/env Rscript
# ============================================================
# 01_LUSC_ExportForCellTypist.R
# GSE243013 — Bloc5_Annotation_TME, Script 01 (LUSC)
#
# Purpose: export the log-normalized expression matrix of the clustered
# LUSC object in MatrixMarket format, as input for CellTypist (run in
# Python by Script 02).
#
# Rationale: a manual export was chosen over an automated Seurat ->
# AnnData converter (e.g. SeuratDisk), to keep every step inspectable
# rather than relying on a black-box conversion. The log-normalized
# "data" layer is exported (NormalizeData, scale.factor = 1e4),
# matching the normalization CellTypist requires for AnnData input
# (log1p, 10,000 counts per cell).
#
# Input:  Objects/01_LUSC_seu_clustered.rds
# Output: Objects/celltypist_export/LUSC_matrix.mtx
#         Objects/celltypist_export/LUSC_barcodes.tsv
#         Objects/celltypist_export/LUSC_features.tsv
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(Matrix)
})

IN_OBJ  <- here("Objects", "01_LUSC_seu_clustered.rds")
OUT_DIR <- here("Objects", "celltypist_export")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

message("Loading clustered Seurat object...")
seu <- readRDS(IN_OBJ)
seu <- JoinLayers(seu)
message("Cells: ", ncol(seu), " | Genes: ", nrow(seu))

# -----------------------------
# 1) Extract log-normalized data (genes x cells, sparse)
# -----------------------------
data_matrix <- GetAssayData(seu, layer = "data")
message("Matrix class: ", class(data_matrix)[1])

# Force materialization to an in-memory sparse matrix if the layer is
# still backed by a lazy/disk-linked object (e.g. class "RenameDims"),
# which lacks a max() method and would fail the sanity check below.
if (!inherits(data_matrix, "dgCMatrix")) {
  message("Converting to dgCMatrix...")
  data_matrix <- as(data_matrix, "dgCMatrix")
}

message("Matrix dims: ", nrow(data_matrix), " genes x ", ncol(data_matrix), " cells")

# Sanity check before export: confirm this is log-normalized (not raw
# counts) — max value should be a modest float, not large integers
message("Max value in matrix (sanity check, expect a small float, not raw counts): ",
        round(max(data_matrix), 2))

# -----------------------------
# 2) Write MatrixMarket format (genes x cells) + barcode/feature lists
# -----------------------------
writeMM(data_matrix, file.path(OUT_DIR, "LUSC_matrix.mtx"))
write.table(colnames(data_matrix), file.path(OUT_DIR, "LUSC_barcodes.tsv"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(rownames(data_matrix), file.path(OUT_DIR, "LUSC_features.tsv"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

message("Saved: LUSC_matrix.mtx, LUSC_barcodes.tsv, LUSC_features.tsv")
message("DONE Script 01 — LUSC export for CellTypist")