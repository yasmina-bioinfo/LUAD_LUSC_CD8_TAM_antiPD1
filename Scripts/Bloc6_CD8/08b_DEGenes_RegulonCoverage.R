#!/usr/bin/env Rscript
# ============================================================
# 08b_DEGenes_RegulonCoverage.R
# GSE243013 — Bloc6_CD8, Script 08b (LUAD and LUSC)
#
# Purpose: Two checks on the CD8.TEX pseudobulk comparison of Script 08,
#          testing whether the absence of significant factors in LUSC
#          follows from a weaker differential signal or from regulon
#          coverage rather than from statistical power: (1) number of
#          differentially expressed genes (DESeq2, MPR vs non-MPR) and
#          median absolute Wald statistic; (2) number of regulons
#          evaluable in CollecTRI and their median size.
#          Same patients and same gene filter (filterByExpr) as
#          Script 08, so the numbers correspond to what run_ulm
#          actually received.
#
# Input:  Objects/03_{LUAD,LUSC}_seu_projectils.rds
#         Data/interactions_collectri_9606.tsv.gz (see README)
# Output: console table; Results/Bloc6_CD8/Tables/08b_DEGenes_RegulonCoverage.csv
# ============================================================

suppressPackageStartupMessages({
  library(here); library(Seurat); library(DESeq2); library(edgeR)
  library(dplyr); library(tibble); library(readr); library(stringr)
})

MIN_CELLS   <- 10
MIN_TARGETS <- 5

# Patients expected from Script 06/08 (guard against a wrong input file)
EXPECTED_PATIENTS <- c(LUAD = 47, LUSC = 89)

# Same network as Script 08 (see Step 5 there for sources): complexes grouped as AP1
# (JUN/FOS) or NFKB (REL/NFKB), no interaction excluded, duplicate pairs removed.
# The sign is not needed here: only the regulator-target pairs are counted.
net <- read_tsv(here("Data", "interactions_collectri_9606.tsv.gz"),
                show_col_types = FALSE) %>%
  mutate(regulator = case_when(
    !str_detect(source, "COMPLEX") ~ source_genesymbol,
    str_detect(source_genesymbol, "JUN") | str_detect(source_genesymbol, "FOS") ~ "AP1",
    str_detect(source_genesymbol, "REL") | str_detect(source_genesymbol, "NFKB") ~ "NFKB")) %>%
  transmute(source = regulator, target = target_genesymbol) %>%
  distinct()
stopifnot(!anyNA(net$source))

check_tex <- function(histology) {
  
  seu <- readRDS(here("Objects", paste0("03_", histology, "_seu_projectils.rds")))
  
  retained <- seu@meta.data %>%
    filter(functional.cluster == "CD8.TEX") %>%
    dplyr::count(sampleID) %>%
    filter(n >= MIN_CELLS) %>%
    pull(sampleID)
  
  message(histology, " — patients retained (>=", MIN_CELLS, " CD8.TEX cells): ", length(retained))
  stopifnot(length(retained) == EXPECTED_PATIENTS[[histology]])
  
  seu_tex <- subset(seu, subset = functional.cluster == "CD8.TEX" &
                      sampleID %in% retained)
  rm(seu); gc()
  
  pb <- AggregateExpression(seu_tex, group.by = "sampleID", assays = "RNA",
                            slot = "counts", return.seurat = FALSE)$RNA
  
  meta <- seu_tex@meta.data %>% as_tibble() %>%
    distinct(sampleID, pathological_response) %>%
    filter(sampleID %in% colnames(pb)) %>%
    arrange(match(sampleID, colnames(pb))) %>%
    column_to_rownames("sampleID")
  
  condition <- factor(meta$pathological_response,
                      levels = c("non-MPR", "MPR"))
  dge <- DGEList(counts = pb, group = condition)
  pb_filt <- pb[filterByExpr(dge, group = condition), ]
  
  dds <- DESeqDataSetFromMatrix(round(pb_filt), meta, ~ pathological_response)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("pathological_response", "MPR", "non-MPR"))
  res <- res[!is.na(res$padj), ]
  
  sizes <- net %>%
    filter(target %in% rownames(pb_filt)) %>%
    dplyr::count(source) %>%
    filter(n >= MIN_TARGETS)
  
  rm(seu_tex); gc()
  
  tibble(histology = histology,
         n_patients = ncol(pb_filt),
         n_genes_kept = nrow(pb_filt),
         n_genes_tested = nrow(res),
         n_padj_05 = sum(res$padj < 0.05),
         n_padj_01 = sum(res$padj < 0.01),
         median_abs_stat = round(median(abs(res$stat), na.rm = TRUE), 3),
         n_regulons_evaluable = nrow(sizes),
         median_regulon_size = median(sizes$n))
}

res_all <- bind_rows(check_tex("LUAD"), check_tex("LUSC"))

cat("\n=== CD8.TEX — differential signal and regulon coverage ===\n")
print(as.data.frame(t(column_to_rownames(res_all, "histology"))))

write_csv(res_all, here("Results", "Bloc6_CD8", "Tables", "08b_DEGenes_RegulonCoverage.csv"))