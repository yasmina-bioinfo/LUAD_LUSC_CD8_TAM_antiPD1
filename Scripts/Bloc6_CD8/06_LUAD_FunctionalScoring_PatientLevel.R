#!/usr/bin/env Rscript
# ============================================================
# 06_LUAD_FunctionalScoring_PatientLevel.R
# GSE243013 — Bloc6_CD8, Script 06 (LUAD)
#
# Purpose: Patient-level comparison of CD8.TEX functional scores
#          (cytotoxicity, exhaustion, stemness) between MPR and non-MPR.
#          Patients with fewer than 10 CD8.TEX cells are excluded
#          (muscat default filtering, Crowell et al. 2020). Per-patient
#          median scores are compared with a Wilcoxon rank-sum test,
#          BH-corrected across the three signatures.
#
# Input:  Objects/05_LUAD_seu_functional_scores.rds
# Output: Results/Bloc6_CD8/Tables/
#           06_LUAD_TEX_patient_level_scores.csv
#           06_LUAD_TEX_patientlevel_wilcoxon.csv
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(scales)
})

IN_SCORES <- here("Objects", "05_LUAD_seu_functional_scores.rds")
OUT_TAB   <- here("Results", "Bloc6_CD8", "Tables")
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

# CD8.TEX cells expected from Script 05 (guard against a wrong input file)
EXPECTED_TEX_CELLS <- 14865  

message("Loading LUAD functional-score object... ", Sys.time())
seu_ex <- readRDS(IN_SCORES)
message("Cells loaded: ", ncol(seu_ex))
stopifnot(ncol(seu_ex) == EXPECTED_TEX_CELLS)

# ---- Patient inclusion: >=10 CD8.TEX cells per patient ----
retained_patients <- seu_ex@meta.data %>%
  filter(functional.cluster == "CD8.TEX") %>%
  count(sampleID) %>%
  filter(n >= 10) %>%
  pull(sampleID)

message("Patients retained (>=10 CD8.TEX cells): ", length(retained_patients))

# ---- Restrict to CD8.TEX cells of retained patients ----
meta <- seu_ex@meta.data %>%
  filter(sampleID %in% retained_patients,
         functional.cluster == "CD8.TEX")

message("CD8.TEX cells retained: ", nrow(meta))

score_cols <- c("Tcell.cytotoxicity", "Tcell.exhaustion", "Tcell.stemness")

# ---- Per-patient median score ----
patient_scores <- meta %>%
  group_by(sampleID, pathological_response) %>%
  summarise(across(all_of(score_cols), median, .names = "median_{.col}"),
            n_cells = n(), .groups = "drop")

print(count(patient_scores, pathological_response))

write_csv(patient_scores,
          file.path(OUT_TAB, "06_LUAD_TEX_patient_level_scores.csv"))

# ---- Patient-level Wilcoxon test, BH-corrected across the three signatures ----
tex_results <- lapply(score_cols, function(sig) {
  col <- paste0("median_", sig)
  mpr_vals    <- patient_scores[[col]][patient_scores$pathological_response == "MPR"]
  nonmpr_vals <- patient_scores[[col]][patient_scores$pathological_response == "non-MPR"]
  test <- wilcox.test(nonmpr_vals, mpr_vals)
  data.frame(
    signature = sig,
    n_MPR = length(mpr_vals), n_nonMPR = length(nonmpr_vals),
    median_MPR = median(mpr_vals), median_nonMPR = median(nonmpr_vals),
    p_value = test$p.value
  )
})
tex_results <- bind_rows(tex_results)
tex_results$p_adj_BH  <- p.adjust(tex_results$p_value, method = "BH")
tex_results$p_display <- scales::pvalue(tex_results$p_adj_BH, accuracy = 0.001, add_p = TRUE)

message("CD8.TEX patient-level results:")
print(tex_results)
write_csv(tex_results,
          file.path(OUT_TAB, "06_LUAD_TEX_patientlevel_wilcoxon.csv"))

message("Script 06 (patient-level functional scoring, LUAD) complete. ", Sys.time())