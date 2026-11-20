#!/usr/bin/env Rscript
# ============================================================
# 01_LUSC_MyeloidCoverage.R
# GSE243013 — Bloc7_TAM, Script 01 (LUSC)
#
# Purpose: Check the composition and per-patient coverage of the myeloid
#          compartment before its classification with scGate (Script 03):
#          which CellTypist labels carry lineage_TAM / lineage_TAM_related,
#          how many TAM cells each patient contributes, and how many
#          patients are retained at candidate cell thresholds.
#          Run once per histology.
#
# Input:  Objects/05_LUSC_seu_annotated.rds
# Output: Results/Bloc7_TAM/Tables/01_LUSC_TAM_cells_per_patient.csv
# ============================================================

library(here); library(Seurat); library(dplyr); library(tibble)

HISTOLOGY <- "LUSC"
seu <- readRDS(here("Objects", "05_LUSC_seu_annotated.rds"))
stopifnot(
  "Cohort mismatch: expected 100 LUSC patients" =
    n_distinct(seu$sampleID) == 100
)

# ---- 1. Composition of the TAM lineage flags ---------------------------
cat("=== annotation_final labels flagged as TAM ===\n")
print(table(seu$annotation_final[seu$lineage_TAM], useNA = "ifany"))

cat("\n=== labels flagged as TAM_related (not TAM) ===\n")
print(table(seu$annotation_final[seu$lineage_TAM_related & !seu$lineage_TAM],
            useNA = "ifany"))

# ---- 2. TAM cells per patient ------------------------------------------
md <- seu@meta.data %>%
  rownames_to_column("cellID") %>%
  filter(lineage_TAM) %>%
  select(cellID, sampleID, pathological_response)

cov <- md %>%
  count(sampleID, pathological_response, name = "n_tam") %>%
  as_tibble() %>%
  arrange(n_tam)
print(cov, n = Inf)

dir.create(here("Results", "Bloc7_TAM", "Tables"), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(cov, here("Results", "Bloc7_TAM", "Tables", "01_LUSC_TAM_cells_per_patient.csv"))

cat("\n=== patients retained at candidate thresholds ===\n")
for (thr in c(10, 20, 30, 50, 100)) {
  kept <- cov %>% filter(n_tam >= thr) %>% count(pathological_response)
  cat("threshold", thr, ":", paste(kept$pathological_response, kept$n,
                                   collapse = " | "), "\n")
}

cat("\nScript 01 (myeloid coverage check,", HISTOLOGY, ") complete.\n")