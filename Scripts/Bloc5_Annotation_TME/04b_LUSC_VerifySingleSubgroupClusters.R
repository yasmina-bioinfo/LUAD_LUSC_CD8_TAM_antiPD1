#!/usr/bin/env Rscript
# ============================================================
# 04b_LUSC_VerifySingleSubgroupClusters.R
# GSE243013 — Bloc5_Annotation_TME, Script 04b (LUSC)
#
# Purpose: complete the label verification of Script 04 for the clusters in
# which only one CellTypist_majority subgroup reaches MIN_CELLS. In these
# clusters no comparison group exists, so FindAllMarkers cannot be run and
# Script 04 writes no marker table for them. The identity of the subgroup is
# verified instead by direct expression: for each official marker, the
# proportion of cells of the subgroup with non-zero log-normalised expression.
# No threshold is applied; the proportions are reported as they are.
#
# Marker rule:
#  - labels with markers given in the CellTypist papers use those markers
#    (Plasma cells: Dominguez Conde et al.; Tem/Trm cytotoxic T cells:
#    Mullan et al. 2023, Table 3);
#  - labels without markers in either paper (Mast cells, pDC) use the 10
#    markers returned by the CellTypist model (Immune_All_Low.pkl) through
#    extract_CellTypist_markers.py, all of them, in the order returned.
#
# Input:  Objects/03_LUSC_seu_celltypist.rds
# Output: Results/Bloc5_Annotation_TME/Tables/04b_LUSC_single_subgroup_marker_check.csv
# ============================================================
suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(data.table)
})

IN_OBJ  <- here("Objects", "03_LUSC_seu_celltypist.rds")
OUT_TAB <- here("Results", "Bloc5_Annotation_TME", "Tables")
dir.create(OUT_TAB, showWarnings = FALSE, recursive = TRUE)

MIN_CELLS <- 30   # same threshold as Script 04

# -----------------------------
# Markers (see rule in the header)
# Model top 10: output of extract_CellTypist_markers.py (Immune_All_Low.pkl),
# run on 2026-11-01, N_MARKERS = 10.
# -----------------------------
marker_table <- bind_rows(
  tibble(label = "Plasma cells",
         gene = c("CD38", "XBP1", "SDC1"),
         source = "Dominguez Conde et al. (CellTypist paper)",
         role = "marker"),
  tibble(label = "Plasma cells",
         gene = "MKI67",
         source = "Dominguez Conde et al. (CellTypist paper)",
         role = "control (expected low in plasma cells)"),
  tibble(label = "Tem/Trm cytotoxic T cells",
         gene = c("GZMK", "CD8A", "CCL5"),
         source = "Mullan et al. 2023, Table 3",
         role = "marker"),
  tibble(label = "Mast cells",
         gene = c("AL157895.1", "TPSAB1", "CPA3", "IL1RL1", "SIGLEC6",
                  "IL1RAPL1", "CLC", "CCN3", "MIR600HG", "TMEM233"),
         source = "CellTypist model, top 10",
         role = "marker"),
  tibble(label = "pDC",
         gene = c("SCT", "AC007381.1", "LRRC26", "C3orf86", "LCNL1",
                  "AC006115.2", "AC022126.1", "TRBV10-2", "MIR600HG", "CD207"),
         source = "CellTypist model, top 10",
         role = "marker")
)

message("Loading CellTypist-annotated Seurat object...")
seu <- readRDS(IN_OBJ)

# -----------------------------
# 1) Clusters with a single subgroup >= MIN_CELLS (detected, not chosen)
# -----------------------------
cluster_label_counts <- table(
  cluster = as.character(seu$seurat_clusters),
  label   = seu$CellTypist_majority
)

single_clusters <- rownames(cluster_label_counts)[
  rowSums(cluster_label_counts >= MIN_CELLS) == 1
]
message("Clusters with a single subgroup >= ", MIN_CELLS, " cells: ",
        paste(sort(as.integer(single_clusters)), collapse = ", "))
stopifnot(
  "Unexpected set of single-subgroup clusters" =
    setequal(single_clusters, c("9", "11", "14", "19"))
)

# -----------------------------
# 2) Direct expression of the markers
# -----------------------------
res <- lapply(single_clusters, function(cl) {
  counts <- cluster_label_counts[cl, ]
  lab <- names(counts)[counts >= MIN_CELLS]
  mk  <- marker_table %>% filter(label == lab)
  stopifnot("No markers defined for this label" = nrow(mk) > 0)
  
  cells <- colnames(seu)[as.character(seu$seurat_clusters) == cl &
                           seu$CellTypist_majority == lab]
  in_object <- mk$gene %in% rownames(seu)
  pct <- rep(NA_real_, nrow(mk))
  if (any(in_object)) {
    expr <- FetchData(seu, vars = mk$gene[in_object], cells = cells, layer = "data")
    pct[in_object] <- round(100 * unname(colMeans(expr > 0)), 2)
  }
  mk %>% mutate(cluster = cl, n_cells = length(cells),
                in_object = in_object, pct_expressing = pct)
}) %>%
  bind_rows() %>%
  select(cluster, label, n_cells, gene, source, role, in_object, pct_expressing) %>%
  arrange(as.integer(cluster))

print(as.data.frame(res), row.names = FALSE)

fwrite(res, file.path(OUT_TAB, "04b_LUSC_single_subgroup_marker_check.csv"))
message("Saved: 04b_LUSC_single_subgroup_marker_check.csv")
message("DONE Script 04b — LUSC single-subgroup cluster verification")