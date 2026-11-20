#!/usr/bin/env Rscript
# ============================================================
# 09_CollecTRI_Myeloid_Table.R
# GSE243013 — Bloc7_TAM, Script 09 (LUAD and LUSC)
#
# Purpose: Build the table of transcription factors significantly
#          associated with response (MPR vs non-MPR) in Macrophage_pure
#          and Monocyte_pure cells, from the
#          CollecTRI results of Script 05. Only factors with adjusted
#          P < 0.05 are listed, with score and adjusted P shown for the
#          population(s) where they are significant (– otherwise).
#          Categories: Shared (significant in both populations),
#          Macrophage_pure only, Monocyte_pure only. Within a category,
#          factors are ordered by adjusted P (for shared factors, by the
#          larger of the two adjusted P values).
#          The script checks that every significant result is listed and
#          that the formatted table matches the source results.
#
# Input:  Results/Bloc7_TAM/Tables/05_{LUAD,LUSC}_CollecTRI_myeloid.csv
# Output: Results/Bloc7_TAM/Tables/09_CollecTRI_Myeloid_table.csv
# ============================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(readr)
})

OUT_TAB <- here("Results", "Bloc7_TAM", "Tables")
GROUPS  <- c("Macrophage_pure", "Monocyte_pure")
CATS    <- c("Shared", "Macrophage_pure only", "Monocyte_pure only")

# Expected number of listed factors (guard against a wrong input file)
EXPECTED_ROWS <- c(LUAD = 72, LUSC = 15)

fs <- function(x) ifelse(is.na(x), "–", sprintf("%.2f", x))
fp <- function(x) ifelse(is.na(x), "–",
                         sapply(x, function(v) format(signif(v, 3), trim = TRUE)))

make_table <- function(h) {
  
  f <- file.path(OUT_TAB, paste0("05_", h, "_CollecTRI_myeloid.csv"))
  if (!file.exists(f)) stop("Input not found: ", f, ". Run Script 05 first.")
  d <- read_csv(f, show_col_types = FALSE)
  
  if (!all(c("source", "score", "group", "histology", "p_adj") %in% names(d)))
    stop("Unexpected columns in ", f)
  if (!all(d$histology == h)) stop("Cohort mismatch in ", f)
  if (!all(GROUPS %in% d$group)) stop("Missing population in ", f)
  if (anyDuplicated(d[, c("group", "source")])) stop("Duplicated factor within a population in ", f)
  
  sig <- d %>% filter(group %in% GROUPS, p_adj < 0.05)
  
  mac <- sig %>% filter(group == "Macrophage_pure") %>%
    select(source, score_mac = score, p_mac = p_adj)
  mon <- sig %>% filter(group == "Monocyte_pure") %>%
    select(source, score_mon = score, p_mon = p_adj)
  
  listed <- full_join(mac, mon, by = "source") %>%
    mutate(category = case_when(
      !is.na(score_mac) & !is.na(score_mon) ~ "Shared",
      !is.na(score_mac)                     ~ "Macrophage_pure only",
      TRUE                                  ~ "Monocyte_pure only"),
      p_order = pmax(p_mac, p_mon, na.rm = TRUE)) %>%
    arrange(factor(category, levels = CATS), p_order, source)
  
  if (nrow(listed) != EXPECTED_ROWS[[h]])
    stop(h, ": ", nrow(listed), " factors listed; expected ", EXPECTED_ROWS[[h]])
  
  out <- listed %>%
    transmute(histology            = h,
              transcription_factor = source,
              Macrophage_pure_score = fs(score_mac),
              Macrophage_pure_p_adj = fp(p_mac),
              Monocyte_pure_score   = fs(score_mon),
              Monocyte_pure_p_adj   = fp(p_mon),
              category)
  
  # Round trip: formatted table -> long format -> compare with source
  back <- bind_rows(
    out %>% filter(Macrophage_pure_score != "–") %>%
      transmute(source = transcription_factor, group = "Macrophage_pure",
                score = as.numeric(Macrophage_pure_score),
                p_adj = as.numeric(Macrophage_pure_p_adj)),
    out %>% filter(Monocyte_pure_score != "–") %>%
      transmute(source = transcription_factor, group = "Monocyte_pure",
                score = as.numeric(Monocyte_pure_score),
                p_adj = as.numeric(Monocyte_pure_p_adj)))
  
  chk <- inner_join(back, sig %>% select(source, group, score, p_adj),
                    by = c("source", "group"), suffix = c("_tab", "_csv"))
  stopifnot(nrow(back) == nrow(sig),                       # nothing missing, nothing extra
            nrow(chk)  == nrow(sig),
            all(abs(chk$score_tab - chk$score_csv) <= 0.0051),
            all(abs(chk$p_adj_tab - chk$p_adj_csv) <= 0.0051 * chk$p_adj_csv))
  
  cat("\n##########", h, "##########\n")
  cat("Significant (adjusted P < 0.05): Macrophage_pure", nrow(mac),
      "| Monocyte_pure", nrow(mon), "| listed factors", nrow(out), "\n")
  print(table(out$category))
  print(as.data.frame(out %>% select(-histology)))
  out
}

res <- bind_rows(make_table("LUAD"), make_table("LUSC"))

write_excel_csv(res, file.path(OUT_TAB, "09_CollecTRI_Myeloid_table.csv"))

cat("\nScript 09 (CollecTRI myeloid table) complete.\n")