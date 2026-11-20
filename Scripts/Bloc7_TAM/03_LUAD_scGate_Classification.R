#!/usr/bin/env Rscript
# ============================================================
# 03_LUAD_scGate_Classification.R
# GSE243013 — Bloc7_TAM, Script 03 (LUAD)
#
# Purpose: Classify the myeloid compartment (lineage_TAM |
#          lineage_TAM_related, CellTypist, Bloc 5) with the scGate
#          TME_HiRes Macrophage and Monocyte models into
#          Macrophage_pure, Monocyte_pure, Ambiguous_double_pure
#          (retained by both models) and scGate_unresolved, then cross
#          this classification with the CellTypist status
#          (Confirmed_Myeloid_Bloc5 / Unresolved_Bloc5 / Other).
#          Run once per histology (see 03_LUSC_scGate_Classification.R).
#
# Design notes:
#   - scGate is run only on lineage_TAM | lineage_TAM_related, not on the
#     full object. For TAM there is no downstream tool (no ProjecTILs
#     equivalent) that requires scGate to stay blind to the CellTypist
#     myeloid gate, unlike CD8 (Bloc 6) where scGate ran on the full
#     object because ProjecTILs was applied afterward. Restricting the
#     input avoids diluting "unresolved" with non-myeloid cells.
#   - No batching: the myeloid subset is small enough to run in one pass.
#   - Both TME_HiRes models are run (Macrophage and Monocyte).
#   - scGate_unresolved is not characterised with additional marker
#     panels: no consensus reference atlas exists for TAM states the way
#     ProjecTILs exists for CD8, so further ad hoc marker checks would add
#     another unarbitrated opinion on top of CellTypist and scGate. This
#     is recorded as a methodological limitation.
#   - Confirmed_Myeloid_Bloc5 / Unresolved_Bloc5 is named verbatim so it
#     is unambiguous which pipeline stage the variable comes from.
#
# Input:  Objects/05_LUAD_seu_annotated.rds
# Output: Objects/03_LUAD_seu_TAM_classified.rds
#         Results/Bloc7_TAM/Tables/03_LUAD_scGate_TAM_class.csv
#         Results/Bloc7_TAM/Tables/03_LUAD_crosstab_Bloc5_vs_scGate.csv
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(scGate)
  library(dplyr)
  library(readr)
  library(tibble)
})

set.seed(42)

HISTOLOGY <- "LUAD"
IN_OBJ    <- here("Objects", "05_LUAD_seu_annotated.rds")
OUT_TAB   <- here("Results", "Bloc7_TAM", "Tables")
OUT_OBJ   <- here("Objects", "03_LUAD_seu_TAM_classified.rds")
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

# Expected size of the lineage_TAM | lineage_TAM_related subset (README, Bloc 5).
EXPECTED_TAM_SUBSET_CELLS <- 33530

# ---- 1. Load and subset ---------------------------------------------------
message("Loading ", HISTOLOGY, " annotated object... ", Sys.time())
seu <- readRDS(IN_OBJ)
message("Total cells in object: ", ncol(seu))

seu_tam <- subset(seu, subset = lineage_TAM | lineage_TAM_related)
message("TAM | TAM_related subset: ", ncol(seu_tam), " cells")
stopifnot(ncol(seu_tam) == EXPECTED_TAM_SUBSET_CELLS)

rm(seu); gc()

# ---- 2. Confirmed_Myeloid_Bloc5 / Unresolved_Bloc5 -----------------------
# Intermediate macrophages is never "confirmed" here: Bloc5 had no
# validated marker set to confirm this state directly, so cells
# labelled Intermediate macrophages are treated as unresolved and
# go through scGate below like any other unresolved myeloid cell.
confirmed_labels <- c("Macrophages", "Classical monocytes", "Non-classical monocytes")

seu_tam$Confirmed_Myeloid_Bloc5 <- case_when(
  grepl("Intermediate macrophages", seu_tam$annotation_final)          ~ "Unresolved_Bloc5",
  seu_tam$annotation_final %in% confirmed_labels                       ~ "Confirmed_Myeloid_Bloc5",
  grepl("^Unresolved", seu_tam$annotation_final)                       ~ "Unresolved_Bloc5",
  seu_tam$annotation_final == "Other (mixed, minor)"                   ~ "Other",
  TRUE ~ NA_character_
)
stopifnot(!anyNA(seu_tam$Confirmed_Myeloid_Bloc5))

cat("\n=== Confirmed_Myeloid_Bloc5 / Unresolved_Bloc5 status ===\n")
print(table(seu_tam$Confirmed_Myeloid_Bloc5))

# ---- 3. scGate, restricted to this subset, no batching -------------------
# scGate is run on a slimmed copy of the object (DietSeurat) to limit memory
# use: the full object was too heavy for the available RAM. The copy keeps the
# harmony reduction (dimreducs = "harmony"); keeping none (NULL) broke
# downstream steps in Bloc 6. scGate is called with its default
# reduction = "calculate": it computes its own PCA for the kNN smoothing
# (k.param = 30) and does not use the harmony embedding.
message("\nRunning scGate (Macrophage + Monocyte models) on the TAM subset... ", Sys.time())

models <- get_scGateDB()
seu_diet <- DietSeurat(seu_tam, dimreducs = "harmony")

seu_diet <- scGate(seu_diet,
                   model  = models$human$TME_HiRes$Macrophage,
                   ncores = 4)
seu_tam$scGate_Macrophage <- seu_diet$is.pure[colnames(seu_tam)]

seu_diet <- scGate(seu_diet,
                   model  = models$human$TME_HiRes$Monocyte,
                   ncores = 4)
seu_tam$scGate_Monocyte <- seu_diet$is.pure[colnames(seu_tam)]

rm(seu_diet); gc()

cat("\n=== scGate raw calls ===\n")
cat("Macrophage model:\n"); print(table(seu_tam$scGate_Macrophage, useNA = "ifany"))
cat("Monocyte model:\n");   print(table(seu_tam$scGate_Monocyte,   useNA = "ifany"))

# ---- 4. scGate_TAM_class ---------------------------------------------------
# Macrophage_pure / Monocyte_pure / scGate_unresolved (neither model retains
# the cell) / Ambiguous_double_pure (both models retain the cell). The double
# case is kept separate rather than folded into scGate_unresolved, because it
# is a different outcome (two identities instead of none).
seu_tam$scGate_TAM_class <- case_when(
  seu_tam$scGate_Macrophage == "Pure" & seu_tam$scGate_Monocyte == "Pure" ~ "Ambiguous_double_pure",
  seu_tam$scGate_Macrophage == "Pure"                                     ~ "Macrophage_pure",
  seu_tam$scGate_Monocyte   == "Pure"                                     ~ "Monocyte_pure",
  TRUE                                                                    ~ "scGate_unresolved"
)

cat("\n=== scGate_TAM_class ===\n")
print(table(seu_tam$scGate_TAM_class))

# ---- 5. Formal crosstab: Confirmed_Myeloid_Bloc5 x scGate_TAM_class -------
crosstab <- table(Confirmed_Myeloid_Bloc5 = seu_tam$Confirmed_Myeloid_Bloc5,
                  scGate_TAM_class       = seu_tam$scGate_TAM_class)

cat("\n=== crosstab: Confirmed_Myeloid_Bloc5 x scGate_TAM_class ===\n")
print(crosstab)

crosstab_df <- as.data.frame.matrix(crosstab) %>%
  rownames_to_column("Confirmed_Myeloid_Bloc5")
write_csv(crosstab_df, file.path(OUT_TAB, "03_LUAD_crosstab_Bloc5_vs_scGate.csv"))

# ---- 6. Per-cell export ----------------------------------------------------
percell <- seu_tam@meta.data %>%
  rownames_to_column("cellID") %>%
  select(cellID, sampleID, pathological_response, annotation_final,
         Confirmed_Myeloid_Bloc5, scGate_Macrophage, scGate_Monocyte, scGate_TAM_class)

write_csv(percell, file.path(OUT_TAB, "03_LUAD_scGate_TAM_class.csv"))

# ---- 7. Save classified object ------------------------------------------
saveRDS(seu_tam, OUT_OBJ)

message("\nScript03 (scGate classification, ", HISTOLOGY, ") complete. ", Sys.time())
