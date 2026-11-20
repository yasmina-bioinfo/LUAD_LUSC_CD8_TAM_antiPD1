#!/usr/bin/env Rscript
# ============================================================
# 06_LUAD_CollecTRI_Overlap.R
# GSE243013 — Bloc7_TAM, Script 06 (LUAD)
#
# Purpose: Identify which transcription factors are significant in both
#          myeloid populations (Macrophage_pure and Monocyte_pure) and
#          which in one only, from the output of Script 05.
#          Run once per histology (see 06_LUSC_CollecTRI_Overlap.R).
#
# Rationale: Script 05 tests each population separately. The question
#          here is whether the significant factors form a common set,
#          which would point to a shared response to the microenvironment,
#          or differ between the two populations. Factors specific to a
#          single population are listed separately, since they are what
#          distinguishes one population from the other.
#
# No new statistical test: the script reads the results of Script 05 and counts
# the factors significant in one or both populations.
#
# Input:  Results/Bloc7_TAM/Tables/05_LUAD_CollecTRI_myeloid.csv
# Output: Results/Bloc7_TAM/Tables/06_LUAD_CollecTRI_overlap.csv
#         (one row per significant factor: number of populations,
#         populations, direction)
# ============================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(readr); library(tibble)
})

HISTOLOGY <- "LUAD"
GROUPS    <- c("Macrophage_pure", "Monocyte_pure")

acts <- read_csv(here("Results", "Bloc7_TAM", "Tables",
                      paste0("05_", HISTOLOGY, "_CollecTRI_myeloid.csv")),
                 show_col_types = FALSE)

sig <- acts %>% filter(p_adj < 0.05)

# ---- 1. In how many groups is each factor significant? -----------------
# The direction column records whether the factor is more active in MPR
# or non-MPR; a factor significant in several groups could in principle
# differ in direction between them, which this makes visible.
shared <- sig %>%
  group_by(source) %>%
  summarise(n_groups  = n_distinct(group),
            groups    = paste(sort(unique(group)), collapse = " + "),
            direction = paste(unique(ifelse(score > 0, "MPR", "nonMPR")),
                              collapse = "/"),
            .groups = "drop")

write_csv(shared, here("Results", "Bloc7_TAM", "Tables",
                       paste0("06_", HISTOLOGY, "_CollecTRI_overlap.csv")))

cat("=== factors by number of groups in which they are significant ===\n")

# ---- 2. The set shared by both populations -------------------------
core <- shared %>% filter(n_groups == 2) %>% pull(source)
cat("\n=== significant in BOTH groups (", length(core), ") ===\n", sep = "")
print(sort(core))

cat("\ndirections within the shared set:\n")
print(shared %>% filter(n_groups == 2) %>% dplyr::count(direction))

# ---- 3. What is specific to each group ---------------------------------
for (g in GROUPS) {
  only <- shared %>% filter(n_groups == 1, groups == g) %>% pull(source)
  cat("\n=== specific to ", g, " (", length(only), ") ===\n", sep = "")
  print(sort(only))
}

# ---- 4. Factors more active in responding tumors -----------------------
# Reported separately because they are few: the asymmetry between response
# groups is itself part of the result.
cat("\n=== more active in MPR ===\n")
sig %>%
  filter(score > 0) %>%
  select(group, source, score, p_adj) %>%
  arrange(group, p_adj) %>%
  print(n = Inf)

cat("\n=== counts by direction and group ===\n")
print(sig %>% dplyr::count(group, direction = ifelse(score > 0, "MPR", "nonMPR")))

cat("\nScript06 (overlap, ", HISTOLOGY, ") complete.\n", sep = "")