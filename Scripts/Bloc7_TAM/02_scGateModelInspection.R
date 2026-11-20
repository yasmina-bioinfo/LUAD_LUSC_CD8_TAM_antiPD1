#!/usr/bin/env Rscript
# ============================================================
# 02_scGateModelInspection.R
# GSE243013 — Bloc7_TAM, Script 02 (LUAD and LUSC)
#
# Purpose: Inspect the scGate TME_HiRes Macrophage and Monocyte models
#          before applying them in Script 03, to confirm which markers
#          define each population, which populations are excluded and at
#          which hierarchical level, and that the two models do not
#          retain the same cells. Both definitions are exported.
#
# Rationale: scGate is the identity-deciding step for the myeloid
#          compartment (Script 03). Applying a model without first
#          reading its definition would make it impossible to interpret
#          which cells are retained or discarded, and why.
#
# Histology-independent: run once, applies to both LUAD and LUSC.
#
# Input:  none (queries the scGate model database)
# Output: Results/Bloc7_TAM/Tables/
#           02_scGate_model_Macrophage.csv
#           02_scGate_model_Monocyte.csv
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(scGate)
  library(dplyr)
  library(readr)
})

OUT_TAB <- here("Results", "Bloc7_TAM", "Tables")
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

models <- get_scGateDB()

# ---- 1. Which model collections and cell types are available? ----------
cat("=== available model collections (human) ===\n")
print(names(models$human))

cat("\n=== cell types in TME_HiRes ===\n")
print(names(models$human$TME_HiRes))

# ---- 2. Definition of the Macrophage model -----------------------------
# Read to confirm: (a) which markers define the population, (b) which
# populations are actively excluded, (c) at which hierarchical level.
cat("\n=== TME_HiRes / Macrophage ===\n")
print(models$human$TME_HiRes$Macrophage)

# ---- 3. Definition of the Monocyte model -------------------------------
# Checked to confirm the two models do not retain the same cells.
cat("\n=== TME_HiRes / Monocyte ===\n")
print(models$human$TME_HiRes$Monocyte)

# ---- 4. Export both definitions for the record -------------------------
write_csv(models$human$TME_HiRes$Macrophage,
          file.path(OUT_TAB, "02_scGate_model_Macrophage.csv"))
write_csv(models$human$TME_HiRes$Monocyte,
          file.path(OUT_TAB, "02_scGate_model_Monocyte.csv"))

cat("\nScript02 (scGate model inspection) complete.\n")