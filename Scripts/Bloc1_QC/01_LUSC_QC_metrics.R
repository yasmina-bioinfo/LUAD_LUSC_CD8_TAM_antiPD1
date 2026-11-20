#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Bloc1_QC, Script 01: Create Seurat object + QC metrics (LUSC)
# Builds the LUSC Seurat object from the BPCells matrix (genes kept
# if detected in at least 3 cells, cells kept if at least 200 genes
# detected, applied by CreateSeuratObject), adds sampleID and
# pathological_response, computes percent.mt and plots QC violins.
# Stops if the sampleIDs differ from the Script 00 reference list.
# Input:  Objects/01_LUSC_matrix_BPCells/
#         Objects/00_LUSC_reference_ids.rds
#         Data/GSE243013_NSCLC_immune_scRNA_metadata.csv.gz
#         Results/Bloc0_DataPrep/Tables/02_LUSC_clinical_metadata.csv
# Output: Objects/01_LUSC_seu_raw.rds
#         Results/Bloc1_QC/Figures/01_LUSC_QC_violin.png
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(BPCells)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

MAT_DIR   <- here("Objects", "01_LUSC_matrix_BPCells")
IDS_FILE  <- here("Objects", "00_LUSC_reference_ids.rds")
CLIN_FILE <- here("Results", "Bloc0_DataPrep", "Tables", "02_LUSC_clinical_metadata.csv")
OUT_OBJ   <- here("Objects")
OUT_FIG   <- here("Results", "Bloc1_QC", "Figures")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)

message("Opening BPCells matrix...")
mat <- open_matrix_dir(MAT_DIR)

message("Creating Seurat object...")
seu <- CreateSeuratObject(counts = mat, project = "GSE243013_LUSC",
                          min.cells = 3, min.features = 200)
message("Seurat object: ", ncol(seu), " cells, ", nrow(seu), " genes")

message("Adding clinical metadata...")
meta_map <- fread(here("Data", "GSE243013_NSCLC_immune_scRNA_metadata.csv.gz"),
                  select = c("cellID", "sampleID"))
meta_map <- meta_map[cellID %in% colnames(seu)]
meta_clin <- fread(CLIN_FILE)
meta_join <- merge(meta_map, meta_clin, by = "sampleID", all.x = TRUE)
meta_join <- meta_join[match(colnames(seu), meta_join$cellID), ]
seu$sampleID <- meta_join$sampleID
lusc_ids <- readRDS(IDS_FILE)
if (!setequal(unique(seu$sampleID), lusc_ids$sampleID)) {
  stop("sampleID in Seurat object do not match Script 00 reference list.")
}
message("sampleID match Script 00 reference list. OK.")
seu$pathological_response <- meta_join$pathological_response

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

p_vln <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                 group.by = "pathological_response", ncol = 3, pt.size = 0) & theme_bw()
ggsave(file.path(OUT_FIG, "01_LUSC_QC_violin.png"), p_vln, width = 12, height = 4, dpi = 300, bg = "white")

saveRDS(seu, file.path(OUT_OBJ, "01_LUSC_seu_raw.rds"))
message("DONE Script 01 — LUSC QC metrics")