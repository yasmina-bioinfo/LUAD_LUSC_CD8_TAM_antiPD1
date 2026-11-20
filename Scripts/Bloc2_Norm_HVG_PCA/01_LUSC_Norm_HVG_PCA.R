#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Bloc2_Norm_HVG_PCA, Script 01: Normalization + HVG + PCA (LUSC)
# LogNormalize (scale factor 1e4), 3000 highly variable genes (vst),
# scaling of the HVGs, PCA on 50 components, and elbow plot used to
# choose the number of PCs for Harmony (Bloc3). The scaled data are
# removed before saving to reduce the object size.
# Input:  Objects/03_LUSC_seu_resp.rds
# Output: Objects/01_LUSC_seu_pca.rds
#         Results/Bloc2_Norm_HVG_PCA/Figures/01_LUSC_ElbowPlot.png
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(ggplot2)
})

IN_OBJ  <- here("Objects", "03_LUSC_seu_resp.rds")
OUT_OBJ <- here("Objects")
OUT_FIG <- here("Results", "Bloc2_Norm_HVG_PCA", "Figures")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)

message("Loading QC Seurat object...")
seu <- readRDS(IN_OBJ)
DefaultAssay(seu) <- "RNA"
message("Cells: ", ncol(seu), " | Genes: ", nrow(seu))

message("Normalizing...")
seu <- NormalizeData(seu, normalization.method = "LogNormalize",
                     scale.factor = 1e4, verbose = FALSE)

# nfeatures = 3000: same as LUAD, see the justification in the LUAD script.

message("Selecting highly variable genes (nfeatures = 3000)...")
seu <- FindVariableFeatures(seu, selection.method = "vst",
                            nfeatures =3000, verbose = FALSE)
message("HVGs selected: ", length(VariableFeatures(seu)))

message("Scaling HVGs...")
seu <- ScaleData(seu, features = VariableFeatures(seu), verbose = FALSE)

message("Running PCA...")
seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = 50, verbose = FALSE)

p_elbow <- ElbowPlot(seu, ndims = 50) + theme_bw() +
  labs(title = "PCA Elbow Plot — GSE243013 LUSC",
       x = "Principal Component", y = "Standard Deviation") +
  theme(plot.title = element_text(size = 12, face = "bold"))
ggsave(file.path(OUT_FIG, "01_LUSC_ElbowPlot.png"), p_elbow,
       width = 6, height = 4, dpi = 300, bg = "white")


seu[["RNA"]]$scale.data <- NULL
gc()

saveRDS(seu, file.path(OUT_OBJ, "01_LUSC_seu_pca.rds"))
message("DONE Script 01 — LUSC Normalization/HVG/PCA")