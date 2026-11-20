#!/usr/bin/env Rscript
# ============================================================
# GSE243013 — Bloc4_Clustering_UMAP, Script 01: Clustering + UMAP (LUSC)
# Builds the neighbour graph, clusters (resolution 0.5) and computes
# the UMAP, all on the Harmony embedding (first 20 dimensions).
# Plots the UMAP by cluster and by patient (integration check).
# Input:  Objects/01_LUSC_seu_harmony.rds
# Output: Objects/01_LUSC_seu_clustered.rds
#         Results/Bloc4_Clustering_UMAP/Figures/01_LUSC_UMAP_clusters.png
#         Results/Bloc4_Clustering_UMAP/Figures/01_LUSC_UMAP_patient.png
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(ggplot2)
})

IN_OBJ  <- here("Objects", "01_LUSC_seu_harmony.rds")
OUT_OBJ <- here("Objects")
OUT_FIG <- here("Results", "Bloc4_Clustering_UMAP", "Figures")
dir.create(OUT_FIG, showWarnings = FALSE, recursive = TRUE)

message("Loading Harmony object...")
seu <- readRDS(IN_OBJ)

message("Finding neighbors (Harmony embedding, 20 dims)...")
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:20, verbose = FALSE)

message("Clustering (resolution = 0.5)...")
seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)
message("Number of clusters: ", length(levels(seu$seurat_clusters)))

message("Running UMAP (Harmony embedding, 20 dims)...")
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:20, verbose = FALSE)

p_clust <- DimPlot(seu, reduction = "umap", label = TRUE, repel = TRUE) +
  theme_bw() + labs(title = "Global TME clustering — GSE243013 LUSC")

p_patient <- DimPlot(seu, reduction = "umap", group.by = "sampleID") +
  theme_bw() + labs(title = "By patient (Harmony integration check)") +
  theme(legend.position = "none")

ggsave(file.path(OUT_FIG, "01_LUSC_UMAP_clusters.png"), p_clust,
       width = 7, height = 5, dpi = 300, bg = "white")
ggsave(file.path(OUT_FIG, "01_LUSC_UMAP_patient.png"), p_patient,
       width = 7, height = 5, dpi = 300, bg = "white")

saveRDS(seu, file.path(OUT_OBJ, "01_LUSC_seu_clustered.rds"))
message("DONE Script 01 — LUSC clustering/UMAP. Clusters: ",
        length(levels(seu$seurat_clusters)))