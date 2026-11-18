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


# ---------------------------------------------------------------------------
# Machine-Readable Study Manifest Generation and Verification
# ---------------------------------------------------------------------------

import hashlib
from datetime import datetime, timezone
from pathlib import Path
from typing import Union

DEFAULT_FILE_DESCRIPTIONS: Dict[str, str] = {
    "celltype_response_summary.csv": "Cross-lineage summary table of complete donor pairs, genes tested, DEG counts, and canonical ISG log2 fold changes",
    "deseq2_B_cells.csv": "Differential expression Wald test statistics and p-values for B cells (~ donor + condition)",
    "deseq2_CD14__Monocytes.csv": "Differential expression Wald test statistics and p-values for CD14+ Monocytes (~ donor + condition)",
    "deseq2_CD4_T_cells.csv": "Differential expression Wald test statistics and p-values for CD4 T cells (~ donor + condition)",
    "deseq2_CD8_T_cells.csv": "Differential expression Wald test statistics and p-values for CD8 T cells (~ donor + condition)",
    "deseq2_Dendritic_cells.csv": "Differential expression Wald test statistics and p-values for Dendritic cells (~ donor + condition)",
    "deseq2_FCGR3A__Monocytes.csv": "Differential expression Wald test statistics and p-values for FCGR3A+ Monocytes (~ donor + condition)",
    "deseq2_NK_cells.csv": "Differential expression Wald test statistics and p-values for NK cells (~ donor + condition)",
    "gsea_hallmark_dotplot.png": "Raster visualization (PNG) of MSigDB Hallmark pathway enrichment in CD14+ Monocytes",
    "gsea_hallmark_dotplot.svg": "Vector graphic (SVG) of MSigDB Hallmark pathway enrichment in CD14+ Monocytes",
    "hallmark_gsea_results.csv": "MSigDB Hallmark Gene Set Enrichment Analysis results (NES, nominal p, FDR q, leading edge genes)",
    "kang_interferon_report.html": "Self-contained standalone HTML interactive study report with embedded figures and summary statistics",
    "string_network.png": "High-resolution raster visualization (PNG) of the 20-node STRING protein-protein interaction network",
    "string_network.svg": "Scalable vector graphic (SVG) of the 20-node STRING protein-protein interaction network",
    "string_network_edges.csv": "Deduplicated high-confidence STRING interaction edges (combined score >= 0.700)",
    "string_network_nodes.csv": "STRING network node degree and functional mechanistic annotations for top interferon effectors",
}


def compute_file_sha256(file_path: Union[str, Path]) -> str:
    """Compute the hexadecimal SHA-256 digest of a file."""
    path = Path(file_path)
    if not path.is_file():
        raise FileNotFoundError(f"File not found: {file_path}")
    hasher = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(65536):
            hasher.update(chunk)
    return hasher.hexdigest()


def build_study_manifest(
    results_dir: Optional[Union[str, Path]] = None,
    custom_params: Optional[Dict[str, Any]] = None,
    timestamp: Optional[str] = None,
) -> Dict[str, Any]:
    """Build a machine-readable manifest recording GEO accession, software versions, parameters, seeds, and output files."""
    manifest = {
        "manifest_version": "1.0.0",
        "generated_at": timestamp or "2026-11-18T17:23:00Z",
        "study_metadata": {
            "study_name": "Kang et al. (2018) PBMC Interferon-Beta Response Analysis",
            "geo_accession": "GSE96583",
            "bioconductor_experimenthub_id": "EH2259",
            "bioconductor_package": "muscData::Kang18_8vs8",
            "organism": "Homo sapiens",
            "taxonomy_id": 9606,
            "platform": "10x Genomics Chromium Single Cell 3'",
            "reference_publication": {
                "title": "Multiplexed droplet single-cell RNA-sequencing using natural genetic variation",
                "authors": "Kang HM, Subramaniam M, Targ S, Nguyen M, Maliskova L, McCarthy E, Wan E, Wong S, Byrnes C, Lanata CM, Gate RE, Mostafavi S, Lum PY, Criswell LA, Raychaudhuri S, Ye CJ",
                "journal": "Nature Biotechnology",
                "volume": 36,
                "pages": "89–94",
                "year": 2018,
                "doi": "10.1038/nbt.4042",
                "pmid": "29227470",
            },
            "cohort": {
                "description": "PBMC aliquots from 8 systemic lupus erythematosus (SLE) patient donors stimulated with recombinant IFN-beta or vehicle control",
                "conditions": ["ctrl", "stim"],
                "treatment": "Recombinant human IFN-beta (100 U/mL) for 6 hours",
                "total_cells_sequenced": 29065,
                "filtered_singlet_cells": 24679,
                "total_features": 35635,
                "donor_count": 8,
                "donor_ids": ["101", "1015", "1016", "1039", "107", "1244", "1256", "1488"],
                "evaluated_cell_types": [
                    "CD14+ Monocytes",
                    "FCGR3A+ Monocytes",
                    "CD4 T cells",
                    "B cells",
                    "Dendritic cells",
                    "NK cells",
                    "CD8 T cells",
                ],
            },
        },
        "package_versions": {
            "python_environment": {
                "python": ">=3.10",
                "scanpy": ">=1.10,<2",
                "anndata": ">=0.10,<1",
                "pydeseq2": "0.4.11",
                "gseapy": "1.1.4",
                "pandas": ">=2.0.0",
                "numpy": ">=1.24.0",
                "scipy": ">=1.11.0",
                "matplotlib": ">=3.8,<4",
                "networkx": ">=3.2",
            },
            "r_bioconductor_environment": {
                "r_base": ">=4.3.0",
                "DESeq2": ">=1.40.0",
                "SingleCellExperiment": ">=1.22.0",
                "SummarizedExperiment": ">=1.30.0",
                "ExperimentHub": ">=2.8.0",
                "muscData": ">=1.14.0",
                "clusterProfiler": ">=4.8.0",
                "msigdbr": ">=7.5.1",
                "Matrix": ">=1.6.0",
            },
            "external_databases": {
                "string_db": {
                    "api_version": "12.0",
                    "api_url": "https://string-db.org/api/json/network",
                    "species": 9606,
                    "min_combined_score": 0.700,
                },
                "msigdb": {
                    "collection": "MSigDB Hallmark Gene Sets (2020)",
                    "pathway_count": 50,
                    "source": "https://www.gsea-msigdb.org/gsea/msigdb",
                },
            },
        },
        "analysis_parameters": {
            "pseudobulk_aggregation": {
                "aggregation_function": "sum",
                "cell_type_col": "cell",
                "donor_col": "ind",
                "condition_col": "stim",
                "reference_condition": "ctrl",
                "stimulated_condition": "stim",
                "min_cells_per_donor_condition": 20,
                "min_complete_donor_pairs": 3,
                "multiplet_filter": "singlet_only",
            },
            "differential_expression": {
                "statistical_framework": "Negative Binomial Generalized Linear Model (GLM)",
                "design_formula": "~ donor + condition",
                "contrast": ["condition", "stim", "ctrl"],
                "test_type": "Wald test",
                "multiple_testing_correction": "Benjamini-Hochberg (BH) FDR",
                "significance_threshold_padj": 0.05,
                "significance_threshold_abs_log2fc": 1.0,
                "gene_filter_min_counts": 10,
                "gene_filter_min_samples": 3,
                "primary_discovery_cell_type": "CD14+ Monocytes",
            },
            "pathway_enrichment": {
                "algorithm": "Pre-ranked Gene Set Enrichment Analysis (GSEA)",
                "ranking_metric": "Paired Wald test statistic (z-score) from CD14+ Monocytes",
                "gene_set_database": "MSigDB Hallmark (50 pathways)",
                "weight_parameter_p": 1.0,
                "min_gene_set_size": 15,
                "max_gene_set_size": 500,
                "permutations": 1000,
            },
            "protein_interaction_network": {
                "source": "STRING Database v12.0",
                "organism_taxid": 9606,
                "score_threshold": 0.700,
                "confidence_tier": "High confidence (>= 0.700)",
                "effector_genes": [
                    "ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "MX2",
                    "OAS1", "OAS2", "OAS3", "OASL", "RSAD2", "USP18",
                    "HERC5", "BST2", "DDX58", "IFIH1", "CXCL10",
                    "CXCL11", "STAT1", "EIF2AK2",
                ],
                "node_count": 20,
                "edge_count": 179,
            },
        },
        "random_seeds": {
            "python_random_seed": 42,
            "numpy_random_seed": 42,
            "gsea_permutation_seed": 42,
            "network_graph_layout_seed": 42,
            "test_harness_seed": 42,
        },
        "output_files": {},
    }

    if custom_params:
        for section, values in custom_params.items():
            if section in manifest and isinstance(manifest[section], dict) and isinstance(values, dict):
                manifest[section].update(values)
            else:
                manifest[section] = values

    if results_dir:
        res_path = Path(results_dir)
        if res_path.exists():
            for p in sorted(res_path.glob("*")):
                # Exclude manifest files to prevent circular hash dependencies
                if p.is_file() and not p.name.endswith("manifest.json"):
                    fmt = (
                        "text/csv"
                        if p.suffix == ".csv"
                        else (
                            "image/png"
                            if p.suffix == ".png"
                            else (
                                "image/svg+xml"
                                if p.suffix == ".svg"
                                else (
                                    "text/html"
                                    if p.suffix == ".html"
                                    else "application/octet-stream"
                                )
                            )
                        )
                    )
                    manifest["output_files"][p.name] = {
                        "file_name": p.name,
                        "relative_path": f"studies/kang_ifnb/results/{p.name}",
                        "format": fmt,
                        "size_bytes": p.stat().st_size,
                        "sha256": compute_file_sha256(p),
                        "description": DEFAULT_FILE_DESCRIPTIONS.get(
                            p.name, "Kang study benchmark artifact"
                        ),
                    }

    return manifest


def verify_study_manifest(
    manifest: Union[Dict[str, Any], str, Path],
    results_dir: Optional[Union[str, Path]] = None,
) -> Tuple[bool, List[str]]:
    """Verify the structure, required fields, and artifact checksums of a study manifest."""
    if isinstance(manifest, (str, Path)):
        p = Path(manifest)
        if not p.is_file():
            return False, [f"Manifest file does not exist: {manifest}"]
        with open(p, "r", encoding="utf-8") as f:
            manifest_dict = json.load(f)
    elif isinstance(manifest, dict):
        manifest_dict = manifest
    else:
        return False, ["Manifest must be a dict, filepath, or Path object"]

    errors = []
    required_sections = [
        "manifest_version",
        "study_metadata",
        "package_versions",
        "analysis_parameters",
        "random_seeds",
        "output_files",
    ]

    for sec in required_sections:
        if sec not in manifest_dict:
            errors.append(f"Missing required section: '{sec}'")

    # Check GEO accession
    study_meta = manifest_dict.get("study_metadata", {})
    if study_meta.get("geo_accession") != "GSE96583":
        errors.append(
            f"Expected GEO accession 'GSE96583', found '{study_meta.get('geo_accession')}'"
        )

    # Check random seeds
    seeds = manifest_dict.get("random_seeds", {})
    for seed_key in ["python_random_seed", "numpy_random_seed"]:
        if seed_key not in seeds:
            errors.append(f"Missing expected random seed: '{seed_key}'")

    # Check output files if results_dir is provided
    if results_dir:
        res_path = Path(results_dir)
        output_files = manifest_dict.get("output_files", {})
        if not output_files:
            errors.append("Manifest output_files section is empty")

        for fname, finfo in output_files.items():
            fpath = res_path / fname
            if not fpath.exists():
                errors.append(f"Output artifact missing on disk: {fname}")
                continue

            expected_size = finfo.get("size_bytes")
            if expected_size is not None and fpath.stat().st_size != expected_size:
                errors.append(
                    f"File size mismatch for {fname}: expected {expected_size}, got {fpath.stat().st_size}"
                )

            expected_sha = finfo.get("sha256")
            if expected_sha is not None:
                actual_sha = compute_file_sha256(fpath)
                if actual_sha != expected_sha:
                    errors.append(
                        f"SHA256 checksum mismatch for {fname}: expected {expected_sha}, got {actual_sha}"
                    )

    return len(errors) == 0, errors


def save_study_manifest(
    manifest: Dict[str, Any], output_path: Union[str, Path]
) -> Path:
    """Save manifest dictionary to a formatted JSON file."""
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
    return out
