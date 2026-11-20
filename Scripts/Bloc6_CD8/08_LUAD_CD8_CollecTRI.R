#!/usr/bin/env Rscript
# ============================================================
# 08_LUAD_CD8_CollecTRI.R
# GSE243013 — Bloc6_CD8, Script 08 (LUAD)
#
# Purpose: Identify transcription factors (TF) whose activity differs
#          between MPR and non-MPR in CD8.TEX. Raw counts are summed
#          per patient (pseudobulk), genes are filtered with edgeR
#          filterByExpr (Chen et al., 2016) and compared with DESeq2.
#          The resulting per-gene Wald statistic is the input of a
#          single decoupleR run_ulm call with the CollecTRI network.
#
# Population: same patients as Script 06 (CD8.TEX cells, patients with
#          >= 10 CD8.TEX cells). Contrast: MPR vs non-MPR, so a positive
#          score means higher activity in responding tumors (same gene
#          filter and convention as Bloc7).
#
# Input:  Objects/03_LUAD_seu_projectils.rds
#         Data/interactions_collectri_9606.tsv.gz (OmniPath table of
#         CollecTRI interactions, human; see README)
# Output: Results/Bloc6_CD8/Tables/
#           08_LUAD_TEX_CollecTRI_pseudobulk_DESeq2.csv
#           08_LUAD_CollecTRI_network_summary.csv
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(DESeq2)
  library(decoupleR)
  library(edgeR)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

IN_OBJ  <- here("Objects", "03_LUAD_seu_projectils.rds")
IN_TCRI <- here("Data", "interactions_collectri_9606.tsv.gz")
OUT_TAB <- here("Results", "Bloc6_CD8", "Tables")
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

# Patients expected from Script 06 (same threshold, same population)
EXPECTED_TEX_PATIENTS <- 47 

# ---- Step 1: load CD8 object, restrict to CD8.TEX cells of patients
# with >=10 TEX cells (identical threshold and patient set as Script 07,
# for consistency between the two functional analyses) ----

message("Loading LUAD CD8 object... ", Sys.time())
seu <- readRDS(IN_OBJ)

retained_patients <- seu@meta.data %>%
  filter(functional.cluster == "CD8.TEX") %>%
  dplyr::count(sampleID) %>%
  filter(n >= 10) %>%
  pull(sampleID)

message("Patients retained (>=10 CD8.TEX cells): ", length(retained_patients))
stopifnot(length(retained_patients) == EXPECTED_TEX_PATIENTS)

seu_tex <- subset(seu, subset = functional.cluster == "CD8.TEX" &
                    sampleID %in% retained_patients)
rm(seu); gc()
message("CD8.TEX cells retained: ", ncol(seu_tex))

# ---- Step 2: pseudobulk aggregation -- sum RAW counts per patient ----
# AggregateExpression() sums raw counts (not normalized data) across all
# cells belonging to the same sampleID; this is exactly what DESeq2's
# statistical model requires as input (Love et al., 2014).

message("Aggregating raw counts per patient (pseudobulk)...")
pseudobulk <- AggregateExpression(
  seu_tex, group.by = "sampleID", assays = "RNA", slot = "counts",
  return.seurat = FALSE
)$RNA
message("Pseudobulk matrix: ", nrow(pseudobulk), " genes x ", ncol(pseudobulk), " patients")

# ---- Step 3: DESeq2 differential analysis, MPR vs non-MPR ----
# Build per-patient metadata (one row per pseudobulk column)
patient_meta <- seu_tex@meta.data %>%
  as_tibble() %>%
  distinct(sampleID, pathological_response) %>%
  filter(sampleID %in% colnames(pseudobulk)) %>%
  arrange(match(sampleID, colnames(pseudobulk))) %>%
  column_to_rownames("sampleID")

# Gene filtering with edgeR's filterByExpr (Chen et al., 2016), the
# recommended strategy for pseudobulk: genes are retained if their counts
# per million exceed a threshold in a number of samples set by the design.
# Using CPM rather than raw counts avoids favouring patients contributing
# more cells. Same filter as Bloc7, so the two blocs are comparable.

condition <- factor(patient_meta$pathological_response,
                    levels = c("non-MPR", "MPR"))
dge <- edgeR::DGEList(counts = pseudobulk, group = condition)
pseudobulk_filt <- pseudobulk[edgeR::filterByExpr(dge, group = condition), ]
message("Genes retained after low-count filtering: ", nrow(pseudobulk_filt),
        " / ", nrow(pseudobulk))

dds <- DESeqDataSetFromMatrix(
  countData = round(pseudobulk_filt),
  colData   = patient_meta,
  design    = ~ pathological_response
)
dds <- DESeq(dds)

# Wald statistic per gene: direction + significance in one number,
# exactly the kind of input the decoupleR authors recommend for run_ulm
# (equivalent to limma's t-value in their bulk RNA-seq vignette).
res <- results(dds, contrast = c("pathological_response", "MPR", "non-MPR"))
message("DESeq2 differential analysis done. ", Sys.time())

deg_stat <- as.data.frame(res) %>%
  filter(!is.na(stat)) %>%
  tibble::rownames_to_column("gene") %>%
  select(gene, stat)

message("Genes with a valid differential statistic: ", nrow(deg_stat))

# ---- Step 4: load the CollecTRI network ----
# The network is not downloaded by this script: see README (CollecTRI file).
if (!file.exists(IN_TCRI)) {
  stop("CollecTRI file not found: ", IN_TCRI, ". Download it as described in the README: ",
       "http://no-tls.static.omnipathdb.org/resources/interactions_collectri_9606.tsv.gz")
}
message("Loading CollecTRI network from ", IN_TCRI)
collectri_raw <- read_tsv(IN_TCRI, show_col_types = FALSE)

# ---- Step 5: CollecTRI network ----
# Network built following the CollecTRI authors and decoupleR::get_collectri():
#  - complexes (COMPLEX in `source`) are kept and renamed AP1 (JUN/FOS) or NFKB (REL/NFKB):
#    saezlab/CollecTRI, scripts/CollecTRI/03_pypath_collectri.py; decoupleR get_collectri(),
#    R/utils-omnipath.R, split_complexes = FALSE (default);
#  - sign: +1 if is_stimulation == 1, -1 if is_stimulation == 0 (get_collectri); no interaction
#    is excluded, those annotated as both stimulation and inhibition receive +1;
#  - duplicate regulator-target pairs removed, first occurrence kept (get_collectri).

net_full <- collectri_raw %>%
  arrange(str_detect(source, "COMPLEX")) %>%
  mutate(regulator = case_when(
    !str_detect(source, "COMPLEX") ~ source_genesymbol,
    str_detect(source_genesymbol, "JUN") | str_detect(source_genesymbol, "FOS") ~ "AP1",
    str_detect(source_genesymbol, "REL") | str_detect(source_genesymbol, "NFKB") ~ "NFKB"))

net <- net_full %>%
  distinct(regulator, target_genesymbol, .keep_all = TRUE) %>%
  transmute(
    source = regulator,
    target = target_genesymbol,
    mor    = case_when(is_stimulation == 1 ~ 1, is_stimulation == 0 ~ -1)
  )
stopifnot(!anyNA(net$source), !anyNA(net$mor))

# Network checks (printed and saved)
net_summary <- bind_rows(
  collectri_raw %>% dplyr::count(is_stimulation, is_inhibition) %>%
    transmute(section = "flags",
              detail = paste0("is_stimulation=", is_stimulation, ";is_inhibition=", is_inhibition),
              n),
  net_full %>% filter(str_detect(source, "COMPLEX")) %>%
    dplyr::count(source, source_genesymbol, regulator) %>%
    transmute(section = "complex_mapping",
              detail = paste0(source, " (", source_genesymbol, ") -> ", regulator),
              n),
  tibble(section = "network",
         detail  = c("interactions in raw table", "interactions after duplicate removal", "regulators"),
         n       = c(nrow(collectri_raw), nrow(net), n_distinct(net$source)))
)
print(net_summary, n = Inf)
write_csv(net_summary, file.path(OUT_TAB, "08_LUAD_CollecTRI_network_summary.csv"))

# ---- Step 6: run_ulm ONCE on the differential statistic vector ----
# mat here is a single-column matrix (genes x 1 "sample", which is really
# the MPR-vs-non-MPR contrast already resolved by DESeq2)

stat_mat <- deg_stat %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()
colnames(stat_mat) <- "MPR_vs_nonMPR"

tf_results <- run_ulm(mat = stat_mat, net = net,
                      .source = 'source', .target = 'target', .mor = 'mor',
                      minsize = 5) %>%
  arrange(p_value) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

message("TFs tested (passed minsize filter): ", nrow(tf_results))
message("Top 10 TF by p-value:")
print(head(tf_results, 10))

write_csv(tf_results, file.path(OUT_TAB, "08_LUAD_TEX_CollecTRI_pseudobulk_DESeq2.csv"))

message("Script 08 (CollecTRI pseudobulk, LUAD) complete. ", Sys.time())