#!/usr/bin/env python3
# ============================================================
# extract_CellTypist_markers.py
# GSE243013 — Bloc5_Annotation_TME, utility script
#
# Not part of the sequential pipeline (no input or output in Objects/ or
# Results/). Queries the trained CellTypist model (Immune_All_Low.pkl)
# directly for the top marker genes it uses internally to recognise a given
# cell type. Used in the label verification (Script 04) when the source
# papers (Mullan et al. 2023; Domínguez Conde et al. 2022) do not list
# explicit markers for a specific label (e.g. Mast cells, pDC).
#
# Usage: change CELL_TYPE below to the CellTypist label to query.
# ============================================================

import celltypist
from celltypist import models

CELL_TYPE = "pDC"  # change to the label being verified
N_MARKERS = 10

models.download_models(model = 'Immune_All_Low.pkl')
model = models.Model.load(model = 'Immune_All_Low.pkl')
top_markers = model.extract_top_markers(CELL_TYPE, N_MARKERS)
print(top_markers)
