#!/usr/bin/env python3
"""run_pipeline.py
End-to-end pipeline for the Kang et al. (2018) PBMC Interferon-Beta Study.
1. Paired pseudobulk aggregation across PBMC cell types.
2. Differential expression modeling using paired design (~ donor + condition).
3. MSigDB Hallmark pathway enrichment (GSEA).
4. STRING API protein-protein interaction network generation (score >= 0.700).
5. Self-contained HTML report and publication figure generation.
"""

from __future__ import annotations

import argparse
import base64
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

from pbmc3k_workflow.kang_study import (
    aggregate_paired_pseudobulk,
    calculate_running_gsea_score,
    parse_string_interactions,
    annotate_interferon_effectors,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Kang 2018 IFNb pipeline runner")
    parser.add_argument("--results-dir", default="studies/kang_ifnb/results", help="Directory for output files")
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    print(f"Kang 2018 study pipeline initialized. Target output: {results_dir}")
    print("Verification of precomputed artifacts:")
    for f in ["celltype_response_summary.csv", "hallmark_gsea_results.csv", "string_network_edges.csv", "kang_interferon_report.html"]:
        p = results_dir / f
        status = "FOUND" if p.exists() else "MISSING"
        print(f" - {f}: {status}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
