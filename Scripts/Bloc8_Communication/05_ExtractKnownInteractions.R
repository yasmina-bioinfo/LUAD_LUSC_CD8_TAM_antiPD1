#!/usr/bin/env Rscript
# ============================================================
# 05_ExtractKnownInteractions.R
# GSE243013 — Bloc8_Communication, Script 05 (LUAD and LUSC)
#
# Purpose: Extract from interaction tables the interactions that are
#          searched for (here, ligand-receptor axes described in the
#          literature), and report for each axis whether it is
#          Present or Absent in each histology and response group.
#          An axis is Present in a group if at least one row of the
#          table has a ligand of the axis, a receptor of the axis, and
#          the sender and receiver cell types given for the axis.
#          Presence does not indicate strength or significance.
#
#          Generic: for another project, change only the parameters
#          below and provide another axes file.
#
# Input:  Inputs/literature_axes.csv
#           columns: axis, ligands, receptors, sender, receiver
#           (several genes separated by ";")
#         Results/Bloc8_Communication/Tables/
#           04_{LUAD,LUSC}_filtered_interactions_FULL.csv
#           columns: group, sender, receiver, ligand, receptor
# Output: Results/Bloc8_Communication/Tables/
#           05_known_axes_interactions.csv  (rows extracted, unchanged,
#                                             with the axis name added)
#           05_known_axes_presence.csv      (Present/Absent per axis)
# ============================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(readr); library(stringr)
})

# ---- Parameters ----
OUT_TAB <- here("Results", "Bloc8_Communication", "Tables")
INTERACTION_FILES <- c(
  LUAD = file.path(OUT_TAB, "04_LUAD_filtered_interactions_FULL.csv"),
  LUSC = file.path(OUT_TAB, "04_LUSC_filtered_interactions_FULL.csv"))
AXES_FILE     <- here("Inputs", "literature_axes.csv")
GROUP_ORDER   <- c("MPR", "non.MPR")
CELL_ALIASES  <- list(myeloid = c("Macrophage_pure", "Monocyte_pure"))
SEP           <- ";"

split_genes  <- function(x) str_trim(str_split(x, fixed(SEP))[[1]])
expand_cells <- function(x) {
  x <- str_trim(x)
  if (x %in% names(CELL_ALIASES)) CELL_ALIASES[[x]] else x
}

# ---- Axes file ----
if (!file.exists(AXES_FILE)) stop("Axes file not found: ", AXES_FILE)
axes <- read_csv(AXES_FILE, col_types = cols(.default = "c"))
need_axes <- c("axis", "ligands", "receptors", "sender", "receiver")
if (!all(need_axes %in% names(axes))) stop("Axes file must have columns: ", paste(need_axes, collapse = ", "))
if (anyNA(axes[need_axes])) stop("Empty cell in the axes file")
if (anyDuplicated(axes$axis)) stop("Duplicated axis name in the axes file")

# ---- Interaction tables ----
need_int <- c("group", "sender", "receiver", "ligand", "receptor")
tabs <- lapply(INTERACTION_FILES, function(f) {
  if (!file.exists(f)) stop("Input not found: ", f)
  d <- read_csv(f, show_col_types = FALSE)
  if (!all(need_int %in% names(d))) stop("Unexpected columns in ", f)
  if (!all(d$group %in% GROUP_ORDER)) stop("Unexpected group values in ", f)
  d
})

# Typo guards: requested cell types and genes must exist somewhere in the tables
cells_in_tables <- unique(unlist(lapply(tabs, function(d) c(d$sender, d$receiver))))
cells_asked <- unique(unlist(lapply(c(axes$sender, axes$receiver), expand_cells)))
if (!all(cells_asked %in% cells_in_tables))
  stop("Cell type(s) not found in the tables: ",
       paste(setdiff(cells_asked, cells_in_tables), collapse = ", "))

genes_in_tables <- unique(unlist(lapply(tabs, function(d) c(d$ligand, d$receptor))))
genes_asked <- unique(unlist(lapply(c(axes$ligands, axes$receptors), split_genes)))
genes_never <- setdiff(genes_asked, genes_in_tables)
if (length(genes_never))
  message("Note: gene(s) never present in the interaction tables (check spelling): ",
          paste(genes_never, collapse = ", "))

# ---- Extraction ----
extracted <- list()
for (h in names(tabs)) {
  for (i in seq_len(nrow(axes))) {
    ax <- axes[i, ]
    L <- split_genes(ax$ligands);  R  <- split_genes(ax$receptors)
    S <- expand_cells(ax$sender);  Rc <- expand_cells(ax$receiver)
    d <- tabs[[h]]
    
    e <- d %>%
      filter(ligand %in% L, receptor %in% R, sender %in% S, receiver %in% Rc) %>%
      mutate(axis = ax$axis, histology = h, .before = 1)
    
    # Independent check: same count through explicit combination keys
    keys_wanted <- do.call(paste, c(expand.grid(L, R, S, Rc, stringsAsFactors = FALSE), sep = "|"))
    keys_table  <- paste(d$ligand, d$receptor, d$sender, d$receiver, sep = "|")
    stopifnot(nrow(e) == sum(keys_table %in% keys_wanted))
    
    extracted[[paste(h, i)]] <- e
  }
}
lines <- bind_rows(extracted)

# ---- Presence / absence ----
counts <- lines %>% count(axis, histology, group, name = "n")
presence <- axes %>% select(axis)
for (h in names(tabs)) {
  for (g in GROUP_ORDER) {
    n <- counts %>% filter(histology == h, group == g)
    v <- n$n[match(axes$axis, n$axis)]
    v[is.na(v)] <- 0
    presence[[paste(h, g, sep = "_")]] <- ifelse(v > 0, "Present", "Absent")
  }
}
# Consistency: Present cells = (axis, histology, group) combinations found in the extracted rows
stopifnot(sum(presence[-1] == "Present") == nrow(counts))

write_csv(lines,    file.path(OUT_TAB, "05_known_axes_interactions.csv"))
write_csv(presence, file.path(OUT_TAB, "05_known_axes_presence.csv"))

cat("\n", nrow(lines), "interactions extracted for", nrow(axes), "axes\n\n")
print(as.data.frame(presence))
cat("\nInteractions per axis, histology and group:\n")
print(as.data.frame(counts), row.names = FALSE)
cat("\nScript 05 (known interactions extraction) complete.\n")