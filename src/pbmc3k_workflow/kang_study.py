"""Kang 2018 PBMC Interferon-Beta paired pseudobulk, GSEA, and STRING network workflow."""

from __future__ import annotations

import json
import urllib.parse
import urllib.request
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd


def aggregate_paired_pseudobulk(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    cell_type_col: str = "cell",
    donor_col: str = "ind",
    condition_col: str = "stim",
    target_cell_type: str = "CD14+ Monocytes",
    reference_condition: str = "ctrl",
    stim_condition: str = "stim",
    min_cells: int = 20,
    min_pairs: int = 3,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Aggregate single cells into donor-paired pseudobulk sum counts.

    Only retains donors with complete pairs (both reference and stim condition
    present with at least min_cells per stratum).
    """
    if len(counts) != len(metadata):
        raise ValueError("Counts and metadata row counts must match")

    # Filter to singlet target cell type
    mask = metadata[cell_type_col] == target_cell_type
    if "multiplets" in metadata.columns:
        mask = mask & (metadata["multiplets"] == "singlet")

    sub_counts = counts.loc[mask]
    sub_meta = metadata.loc[mask]

    if len(sub_counts) == 0:
        raise ValueError(f"No cells found for target cell type: {target_cell_type}")

    # Count cells per donor and condition
    cell_counts = (
        sub_meta.groupby([donor_col, condition_col], observed=False)
        .size()
        .unstack(fill_value=0)
    )

    if reference_condition not in cell_counts.columns or stim_condition not in cell_counts.columns:
        raise ValueError("Missing reference or stim condition in metadata")

    # Find valid donors with at least min_cells in both conditions
    valid_donors = cell_counts.index[
        (cell_counts[reference_condition] >= min_cells)
        & (cell_counts[stim_condition] >= min_cells)
    ].tolist()

    if len(valid_donors) < min_pairs:
        raise ValueError(
            f"Fewer than {min_pairs} complete donor pairs with >= {min_cells} cells: {len(valid_donors)} found"
        )

    pb_counts = {}
    pb_meta_rows = []

    for donor in sorted(valid_donors):
        for cond in [reference_condition, stim_condition]:
            sample_id = f"{donor}_{cond}"
            idx = sub_meta.index[
                (sub_meta[donor_col] == donor) & (sub_meta[condition_col] == cond)
            ]
            cell_sum = sub_counts.loc[idx].sum(axis=0)
            pb_counts[sample_id] = cell_sum
            pb_meta_rows.append(
                {
                    "sample": sample_id,
                    "donor": str(donor),
                    "condition": cond,
                    "n_cells": len(idx),
                    "cell_type": target_cell_type,
                }
            )

    pb_counts_df = pd.DataFrame(pb_counts).T
    pb_meta_df = pd.DataFrame(pb_meta_rows).set_index("sample")

    return pb_counts_df, pb_meta_df


def calculate_running_gsea_score(
    ranked_genes: pd.Series, gene_set: Set[str], p_weight: float = 1.0
) -> Dict[str, Any]:
    """Calculate running enrichment score (ES) and leading edge genes for a gene set."""
    if len(ranked_genes) == 0 or len(gene_set) == 0:
        raise ValueError("Ranked genes and gene set must be non-empty")

    genes = ranked_genes.index.tolist()
    stats = ranked_genes.values

    in_set = np.array([g in gene_set for g in genes], dtype=bool)
    n_hit = int(np.sum(in_set))
    n_miss = len(genes) - n_hit

    if n_hit == 0:
        return {"es": 0.0, "lead_genes": [], "tag_fraction": 0.0}

    weights = np.abs(stats) ** p_weight
    norm_hit = np.sum(weights[in_set])
    if norm_hit == 0:
        norm_hit = 1.0

    p_hit = np.where(in_set, weights / norm_hit, 0.0)
    p_miss = np.where(~in_set, 1.0 / n_miss, 0.0)

    running_es = np.cumsum(p_hit - p_miss)
    max_es = float(np.max(running_es))
    min_es = float(np.min(running_es))
    es = max_es if abs(max_es) >= abs(min_es) else min_es

    peak_idx = int(np.argmax(running_es))
    lead_genes = [genes[i] for i in range(peak_idx + 1) if in_set[i]]

    return {
        "es": es,
        "lead_genes": lead_genes,
        "tag_fraction": n_hit / len(gene_set),
    }


def parse_string_interactions(
    interactions: List[Dict[str, Any]], min_score: float = 0.700
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Parse raw STRING interaction dictionaries into deduplicated edges and node degree metrics."""
    if not interactions:
        return (
            pd.DataFrame(columns=["source", "target", "score"]),
            pd.DataFrame(columns=["gene", "degree"]),
        )

    edge_set = set()
    edges = []
    for item in interactions:
        p_a = item.get("preferredName_A") or item.get("source")
        p_b = item.get("preferredName_B") or item.get("target")
        score = float(item.get("score") or item.get("combined_score", 0.0))

        if score < min_score or p_a == p_b:
            continue

        key = tuple(sorted([str(p_a), str(p_b)]))
        if key not in edge_set:
            edge_set.add(key)
            edges.append({"source": key[0], "target": key[1], "score": score})

    if not edges:
        return (
            pd.DataFrame(columns=["source", "target", "score"]),
            pd.DataFrame(columns=["gene", "degree"]),
        )

    edges_df = pd.DataFrame(edges)
    all_nodes = set(edges_df["source"]).union(set(edges_df["target"]))
    degrees = []
    for node in sorted(all_nodes):
        deg = int(((edges_df["source"] == node) | (edges_df["target"] == node)).sum())
        degrees.append({"gene": node, "degree": deg})

    nodes_df = pd.DataFrame(degrees).sort_values("degree", ascending=False).reset_index(drop=True)
    return edges_df, nodes_df


def annotate_interferon_effectors(nodes_df: pd.DataFrame) -> pd.DataFrame:
    """Annotate effector genes with functional mechanistic clusters."""
    clusters = {
        "ISGylation cascade": {"ISG15", "USP18", "HERC5", "UBE2L6"},
        "Viral RNA recognition (IFIT)": {"IFIT1", "IFIT2", "IFIT3", "IFIT5"},
        "2'-5'-OAS / RNase L pathway": {"OAS1", "OAS2", "OAS3", "OASL"},
        "Dynamin-like GTPase antiviral trap": {"MX1", "MX2"},
        "Viperin nucleotide synthesis": {"RSAD2", "CMPK2"},
        "Pro-inflammatory chemokine": {"CXCL10", "CXCL11", "CXCL9", "CCL8"},
        "Cytosolic RNA sensor": {"DDX58", "IFIH1"},
        "Translation blockade": {"EIF2AK2"},
        "JAK-STAT transcription factor": {"STAT1", "STAT2", "IRF7", "IRF9"},
        "Viral tethering / release inhibitor": {"BST2"},
    }

    df = nodes_df.copy()
    assigned = []
    for gene in df["gene"]:
        c_label = "Antiviral effector"
        for cname, gene_set in clusters.items():
            if gene in gene_set:
                c_label = cname
                break
        assigned.append(c_label)

    df["functional_cluster"] = assigned
    return df
