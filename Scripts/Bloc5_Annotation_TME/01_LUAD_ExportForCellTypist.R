#!/usr/bin/env Rscript
# ============================================================
# 01_LUAD_ExportForCellTypist.R
# GSE243013 — Bloc5_Annotation_TME, Script 01 (LUAD)
#
# Purpose: export the log-normalized expression matrix of the clustered
# LUAD object in MatrixMarket format, as input for CellTypist (run in
# Python by Script 02).
#
# Rationale: a manual export was chosen over an automated Seurat ->
# AnnData converter (e.g. SeuratDisk), to keep every step inspectable
# rather than relying on a black-box conversion. The log-normalized
# "data" layer is exported (NormalizeData, scale.factor = 1e4),
# matching the normalization CellTypist requires for AnnData input
# (log1p, 10,000 counts per cell).
#
# Input:  Objects/01_LUAD_seu_clustered.rds
# Output: Objects/celltypist_export/LUAD_matrix.mtx
#         Objects/celltypist_export/LUAD_barcodes.tsv
#         Objects/celltypist_export/LUAD_features.tsv
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(Matrix)
})

IN_OBJ  <- here("Objects", "01_LUAD_seu_clustered.rds")
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
writeMM(data_matrix, file.path(OUT_DIR, "LUAD_matrix.mtx"))
write.table(colnames(data_matrix), file.path(OUT_DIR, "LUAD_barcodes.tsv"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(rownames(data_matrix), file.path(OUT_DIR, "LUAD_features.tsv"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

message("Saved: LUAD_matrix.mtx, LUAD_barcodes.tsv, LUAD_features.tsv")
message("DONE Script 01 — LUAD export for CellTypist")