#!/usr/bin/env Rscript
# fetch_kang_data.R
# Fetches the Kang et al. (2018) 10x droplet PBMC dataset from Bioconductor
# ExperimentHub (EH2259 / muscData::Kang18_8vs8) and prepares singlet counts.

suppressPackageStartupMessages({
  if (requireNamespace("ExperimentHub", quietly = TRUE)) {
    eh <- ExperimentHub::ExperimentHub()
    sce <- eh[["EH2259"]]
  } else if (requireNamespace("muscData", quietly = TRUE)) {
    sce <- muscData::Kang18_8vs8()
  } else {
    stop("Please install 'ExperimentHub' or 'muscData' to fetch Kang18_8vs8 data.")
  }
})

cat("Successfully loaded Kang 2018 SingleCellExperiment object.\n")
cat("Dimensions:", nrow(sce), "genes x", ncol(sce), "cells\n")

# Filter singlets
cd <- SummarizedExperiment::colData(sce)
singlet_idx <- which(cd$multiplets == "singlet" & !is.na(cd$cell))
sce_singlet <- sce[, singlet_idx]

cat("Retained", ncol(sce_singlet), "singlet cells across", length(unique(cd$ind[singlet_idx])), "donors.\n")
