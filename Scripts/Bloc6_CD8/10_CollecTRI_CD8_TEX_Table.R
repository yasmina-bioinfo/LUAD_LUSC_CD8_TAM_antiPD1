#!/usr/bin/env Rscript
# ============================================================
# 10_CollecTRI_CD8_TEX_Table.R
# GSE243013 — Bloc6_CD8, Script 10 (LUAD and LUSC)
#
# Purpose: Build the table of transcription factors with the smallest P
#          values in CD8.TEX cells (MPR vs non-MPR), from the CollecTRI
#          results of Script 08.
#          The twelve factors with the smallest P values are listed per
#          histology, grouped by direction (higher in non-MPR first, then
#          higher in MPR) and ordered by P value within each group;
#          factors with adjusted P < 0.05 are marked with an asterisk.
#          The script checks that every significant factor is listed and
#          that the formatted table matches the source results.
#
# Input:  Results/Bloc6_CD8/Tables/
#           08_{LUAD,LUSC}_TEX_CollecTRI_pseudobulk_DESeq2.csv
# Output: Results/Bloc6_CD8/Tables/10_CollecTRI_CD8_TEX_table.csv
# ============================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(readr)
})

OUT_TAB  <- here("Results", "Bloc6_CD8", "Tables")
N_LISTED <- 12

# Expected from Script 08 (guard against a wrong input file)
EXPECTED_TESTED <- c(LUAD = 479, LUSC = 461)
EXPECTED_SIG    <- c(LUAD = 12,  LUSC = 0)

make_table <- function(h) {
  
  f <- file.path(OUT_TAB, paste0("08_", h, "_TEX_CollecTRI_pseudobulk_DESeq2.csv"))
  if (!file.exists(f)) stop("Input not found: ", f, ". Run Script 08 first.")
  d <- read_csv(f, show_col_types = FALSE)
  
  if (!all(c("source", "score", "p_value", "p_adj_BH") %in% names(d)))
    stop("Unexpected columns in ", f)
  if (anyDuplicated(d$source)) stop("Duplicated transcription factors in ", f)
  
  n_tested <- nrow(d)
  n_sig    <- sum(d$p_adj_BH < 0.05)
  if (n_tested != EXPECTED_TESTED[[h]] || n_sig != EXPECTED_SIG[[h]])
    stop(h, ": ", n_tested, " factors tested and ", n_sig, " significant; expected ",
         EXPECTED_TESTED[[h]], " and ", EXPECTED_SIG[[h]])
  if (n_sig > N_LISTED)
    stop(h, ": more significant factors than listed rows; raise N_LISTED")
  
  # Selection: the N_LISTED smallest P values, then grouped by direction
  listed <- d %>%
    arrange(p_value) %>%
    slice_head(n = N_LISTED) %>%
    mutate(higher_in = ifelse(score < 0, "non-MPR", "MPR")) %>%
    arrange(factor(higher_in, levels = c("non-MPR", "MPR")), p_value)
  
  # Every significant factor must be in the table
  stopifnot(all(d$source[d$p_adj_BH < 0.05] %in% listed$source))
  
  # Formatting: 2 decimals for the score, 3 significant digits for the P value
  fp <- function(x) sapply(x, function(v) format(signif(v, 3), trim = TRUE))
  out <- listed %>%
    transmute(histology            = h,
              transcription_factor = paste0(source, ifelse(p_adj_BH < 0.05, "*", "")),
              score                = sprintf("%.2f", score),
              adjusted_P           = fp(p_adj_BH),
              higher_in)
  
  # Round trip: the formatted table must match the source results
  back <- out %>%
    transmute(source   = sub("\\*$", "", transcription_factor),
              score    = as.numeric(score),
              p_adj_BH = as.numeric(adjusted_P))
  chk <- inner_join(back, d %>% select(source, score, p_adj_BH),
                    by = "source", suffix = c("_tab", "_csv"))
  stopifnot(nrow(chk) == N_LISTED,
            all(abs(chk$score_tab - chk$score_csv) <= 0.0051),
            all(abs(chk$p_adj_BH_tab - chk$p_adj_BH_csv) <= 0.0051 * chk$p_adj_BH_csv),
            sum(grepl("\\*$", out$transcription_factor)) == n_sig)
  
  cat("\n##########", h, "##########\n")
  cat(n_tested, "factors tested,", n_sig, "significant,", nrow(out), "listed\n")
  print(as.data.frame(out %>% select(-histology)))
  out
}

res <- bind_rows(make_table("LUAD"), make_table("LUSC"))

write_csv(res, file.path(OUT_TAB, "10_CollecTRI_CD8_TEX_table.csv"))

cat("\nScript 10 (CollecTRI CD8.TEX table) complete.\n")