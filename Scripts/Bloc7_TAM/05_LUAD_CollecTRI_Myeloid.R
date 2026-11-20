#!/usr/bin/env Rscript
# ============================================================
# 05_LUAD_CollecTRI_Myeloid.R
# GSE243013 — Bloc7_TAM, Script 05 (LUAD)
#
# Purpose: Infer transcription factor activity in the myeloid compartment,
#          comparing MPR and non-MPR within each of the two populations
#          retained by scGate (Script 03): Macrophage_pure and
#          Monocyte_pure.
#
# Method:  Pseudobulk (raw counts summed per patient; patients with at
#          least 10 cells in the population). Genes are filtered with
#          edgeR's filterByExpr (Chen et al., 2016), the recommended
#          strategy for pseudobulk: it keeps genes whose counts per million
#          exceed a threshold in a sufficient number of samples, which
#          avoids favouring patients contributing more cells. Differential
#          expression MPR vs non-MPR with DESeq2 (Love et al., 2014); the
#          Wald statistic of each gene is the input of decoupleR run_ulm
#          (CollecTRI network built as in Bloc6 Script 08, minsize = 5).
#
# Input:  Objects/03_LUAD_seu_TAM_classified.rds
#         Data/interactions_collectri_9606.tsv.gz  (see README)
# Output: Results/Bloc7_TAM/Tables/05_LUAD_CollecTRI_myeloid.csv
# ============================================================

suppressPackageStartupMessages({
  library(here); library(Seurat); library(dplyr); library(tibble)
  library(readr); library(stringr); library(DESeq2); library(decoupleR)
  library(edgeR)
})

set.seed(42)

HISTOLOGY <- "LUAD"
OUT_TAB   <- here("Results", "Bloc7_TAM", "Tables")
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

MIN_CELLS    <- 10   # per patient, within each group
MIN_PATIENTS <- 4    # per response group, to run the contrast

# ---- 1. Build the two-group column -------------------------------------

message("Loading TAM-classified object (Bloc7)... ", Sys.time())
seu <- readRDS(here("Objects", paste0("03_", HISTOLOGY, "_seu_TAM_classified.rds")))

seu$myeloid_group <- NA_character_
seu$myeloid_group[seu$scGate_TAM_class == "Macrophage_pure"] <- "Macrophage_pure"
seu$myeloid_group[seu$scGate_TAM_class == "Monocyte_pure"]   <- "Monocyte_pure"

cat("\n=== cells per group and response ===\n")
print(table(seu$myeloid_group, seu$pathological_response, useNA = "no"))

# ---- 2. CollecTRI network ----------------------------------------------
# The network is not downloaded by this script: see README (CollecTRI file).
IN_TCRI <- here("Data", "interactions_collectri_9606.tsv.gz")
if (!file.exists(IN_TCRI)) {
  stop("CollecTRI file not found: ", IN_TCRI, ". Download it as described in the README: ",
       "http://no-tls.static.omnipathdb.org/resources/interactions_collectri_9606.tsv.gz")
}
collectri_raw <- read_tsv(IN_TCRI, show_col_types = FALSE)

# Same network as Bloc6, Script 08 (Step 5 there gives the sources): complexes grouped
# as AP1 (JUN/FOS) or NFKB (REL/NFKB), no interaction excluded, sign +1 if
# is_stimulation == 1 and -1 otherwise, duplicate regulator-target pairs removed.
net <- collectri_raw %>%
  arrange(str_detect(source, "COMPLEX")) %>%
  mutate(regulator = case_when(
    !str_detect(source, "COMPLEX") ~ source_genesymbol,
    str_detect(source_genesymbol, "JUN") | str_detect(source_genesymbol, "FOS") ~ "AP1",
    str_detect(source_genesymbol, "REL") | str_detect(source_genesymbol, "NFKB") ~ "NFKB")) %>%
  distinct(regulator, target_genesymbol, .keep_all = TRUE) %>%
  transmute(source = regulator,
            target = target_genesymbol,
            mor    = case_when(is_stimulation == 1 ~ 1, is_stimulation == 0 ~ -1))
stopifnot(!anyNA(net$source), !anyNA(net$mor))

cat("\nCollecTRI: ", n_distinct(net$source), " TFs, ",
    nrow(net), " interactions\n", sep = "")

# ---- 3. One analysis per group -----------------------------------------
run_group <- function(group_name) {
  
  cat("\n\n########## ", group_name, " ##########\n")
  
  cells <- colnames(seu)[!is.na(seu$myeloid_group) &
                           seu$myeloid_group == group_name]
  md <- seu@meta.data[cells, c("sampleID", "pathological_response")]
  
  keep <- md %>%
    rownames_to_column("cell") %>%
    dplyr::count(sampleID, pathological_response) %>%
    filter(n >= MIN_CELLS)
  
  cat("patients retained: ")
  print(table(keep$pathological_response))
  
  if (min(table(keep$pathological_response)) < MIN_PATIENTS) {
    cat("!! fewer than ", MIN_PATIENTS,
        " patients in one group, contrast not run\n", sep = "")
    return(NULL)
  }
  
  cells <- cells[md$sampleID %in% keep$sampleID]
  
  # pseudobulk: raw counts summed per patient
  counts <- GetAssayData(seu, assay = "RNA", layer = "counts")[, cells]
  pb <- sapply(keep$sampleID, function(s) {
    idx <- which(seu$sampleID[match(cells, colnames(seu))] == s)
    Matrix::rowSums(counts[, idx, drop = FALSE])
  })
  colnames(pb) <- keep$sampleID
  
  condition <- factor(keep$pathological_response,
                      levels = c("non-MPR", "MPR"))
  coldata <- data.frame(row.names = keep$sampleID, condition = condition)
  
  # ---- gene filtering with filterByExpr ----
  # A DGEList is built solely to run the filter; the retained genes are
  # then passed to DESeq2, which performs the test itself.
  dge  <- DGEList(counts = pb, group = condition)
  keep_genes <- filterByExpr(dge, group = condition)
  pb <- pb[keep_genes, ]
  cat("genes retained (filterByExpr): ", nrow(pb), "\n", sep = "")
  
  dds <- DESeqDataSetFromMatrix(round(pb), coldata, ~ condition)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", "MPR", "non-MPR"))
  
  stat_mat <- as.matrix(data.frame(stat = res$stat,
                                   row.names = rownames(res)))
  stat_mat <- stat_mat[!is.na(stat_mat[, 1]), , drop = FALSE]
  
  acts <- run_ulm(mat = stat_mat, net = net,
                  .source = "source", .target = "target", .mor = "mor",
                  minsize = 5)
  
  acts <- acts %>%
    mutate(group = group_name, histology = HISTOLOGY,
           p_adj = p.adjust(p_value, method = "BH")) %>%
    arrange(p_adj)
  
  cat("\nTFs with adjusted P < 0.05: ", sum(acts$p_adj < 0.05), "\n", sep = "")
  cat("  higher in MPR:     ", sum(acts$p_adj < 0.05 & acts$score > 0), "\n", sep = "")
  cat("  higher in non-MPR: ", sum(acts$p_adj < 0.05 & acts$score < 0), "\n", sep = "")
  
  print(acts %>% filter(p_adj < 0.05) %>% head(15))
  
  acts
}

results_all <- bind_rows(
  run_group("Macrophage_pure"),
  run_group("Monocyte_pure")
)

write_csv(results_all,
          file.path(OUT_TAB, paste0("05_", HISTOLOGY, "_CollecTRI_myeloid.csv")))

cat("\nScript05 (CollecTRI myeloid, ", HISTOLOGY, ") complete.\n", sep = "")