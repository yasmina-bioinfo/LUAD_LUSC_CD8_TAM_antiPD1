#!/usr/bin/env Rscript
# ============================================================
# 04_LUSC_FilterInteractions.R
# GSE243013 — Bloc8_Communication, Script 04 (LUSC)
#
# Purpose: Select the prioritised interactions of the MultiNicheNet output
#          (get_top_n_lr_pairs, top_n = 10000, i.e. no effective cut-off)
#          and restrict them to CD8.TEX-myeloid interactions.
#          Run once per histology (see 04_LUSC_FilterInteractions.R).
#
# Filters:
#   - Pairs with the same sender and receiver cell type are removed: in
#     dissociated data, a signal sent to the same cell cannot be
#     distinguished from one sent to a neighbouring cell of the same type.
#   - Only interactions with CD8.TEX as sender or receiver are kept. This
#     also removes macrophage-monocyte pairs, which are outside the
#     question (communication between CD8.TEX and the myeloid compartment).
#
# Input:  Objects/03_LUSC_multinichenet_output.rds
# Output: Results/Bloc8_Communication/Tables/
#           04_LUSC_filtered_interactions_FULL.csv
# ============================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(readr); library(multinichenetr); library(tidyr)
})

HISTOLOGY <- "LUSC"

OUT_TAB <- here("Results", "Bloc8_Communication", "Tables")

# ---- 1. Load the MultiNicheNet output ---------------------------------
f_out <- here("Objects", paste0("03_", HISTOLOGY, "_multinichenet_output.rds"))
stopifnot(file.exists(f_out))
output <- readRDS(f_out)

stopifnot(all(c("group", "top_group", "sender", "receiver", "ligand", "receptor",
                "id", "prioritization_score",
                "fraction_expressing_ligand_receptor") %in%
                colnames(output$prioritization_tables$group_prioritization_tbl)),
          setequal(unique(output$prioritization_tables$group_prioritization_tbl$group),
                   c("MPR", "non.MPR")))

# ---- 2. Select and filter ----------------------------------------------
all_pairs <- get_top_n_lr_pairs(output$prioritization_tables, top_n = 10000)

# Remove same-cell-type pairs and keep only CD8.TEX <-> myeloid interactions (research question)
filtered <- all_pairs %>% ungroup() %>%
  filter(sender != receiver,
         sender == "CD8.TEX" | receiver == "CD8.TEX")

full <- filtered %>%
  group_by(group) %>%
  arrange(desc(prioritization_score), .by_group = TRUE) %>%
  mutate(rank_non_autocrine = row_number()) %>%
  ungroup()

cat("all prioritised                      :", nrow(all_pairs), "\n")
cat("same sender and receiver (autocrine) :", sum(all_pairs$sender == all_pairs$receiver), "\n")
cat("after removing autocrine pairs       :", sum(all_pairs$sender != all_pairs$receiver), "\n")
cat("CD8.TEX-myeloid                      :", nrow(full), "\n")
print(table(full$group))

# ---- 3. Direction of the dialogue -------------------------------------
cat("\n=== interactions by direction ===\n")
print(as.data.frame(
  full %>%
    dplyr::count(group, sender, receiver) %>%
    tidyr::pivot_wider(names_from = group, values_from = n, values_fill = 0) %>%
    arrange(desc(non.MPR))
))

# ---- 4. Export ----------------------------------------------------------
write_csv(full, file.path(OUT_TAB,
                          paste0("04_", HISTOLOGY, "_filtered_interactions_FULL.csv")))

cat("\nScript 04 (filter interactions,", HISTOLOGY, ") complete.\n")