"""Download, preprocess, cluster, and summarize the PBMC 3k dataset."""
from __future__ import annotations

import argparse
import csv
import json
import tarfile
import urllib.request
from pathlib import Path

DATA_URL = "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"


def validate_parameters(min_genes: int, max_mito_pct: float, n_top_genes: int, resolution: float) -> None:
    if min_genes < 0:
        raise ValueError("min_genes must be non-negative")
    if not 0 <= max_mito_pct <= 100:
        raise ValueError("max_mito_pct must be between 0 and 100")
    if n_top_genes < 1:
        raise ValueError("n_top_genes must be at least 1")
    if resolution <= 0:
        raise ValueError("resolution must be positive")


def download_data(data_dir: Path) -> Path:
    data_dir.mkdir(parents=True, exist_ok=True)
    archive = data_dir / "pbmc3k_filtered_gene_bc_matrices.tar.gz"
    matrix_dir = data_dir / "filtered_gene_bc_matrices" / "hg19"
    if (matrix_dir / "matrix.mtx").exists():
        return matrix_dir
    if not archive.exists():
        print(f"Downloading public 10x data from {DATA_URL}")
        urllib.request.urlretrieve(DATA_URL, archive)
    with tarfile.open(archive, "r:gz") as tar:
        root = data_dir.resolve()
        for member in tar.getmembers():
            destination = (data_dir / member.name).resolve()
            if root not in destination.parents and destination != root:
                raise ValueError(f"Unsafe archive path: {member.name}")
            if not (member.isdir() or member.isfile()):
                raise ValueError(f"Unsupported archive member type: {member.name}")
        for member in tar.getmembers():
            destination = (data_dir / member.name).resolve()
            if member.isdir():
                destination.mkdir(parents=True, exist_ok=True)
            elif member.isfile():
                destination.parent.mkdir(parents=True, exist_ok=True)
                source = tar.extractfile(member)
                if source is None:
                    raise OSError(f"Unable to read archive member: {member.name}")
                with source, destination.open("wb") as target:
                    target.write(source.read())
    if not (matrix_dir / "matrix.mtx").exists():
        raise FileNotFoundError(f"Expected matrix not found under {matrix_dir}")
    return matrix_dir



HIGH_SCORE = 0.5
HIGH_MARGIN = 0.25
MODERATE_SCORE = 0.2
MODERATE_MARGIN = 0.1
MIN_MARKERS_DETECTED = 2


def classify_marker_scores(scores: dict[str, float], detected: dict[str, int]) -> dict[str, object]:
    """Assign a heuristic label; scores are mean log1p-normalized marker expression, not probabilities."""
    if not scores:
        return {"cell_type": "Unassigned", "confidence": "Unassigned", "top_score": 0.0, "runner_up": "", "margin": 0.0}
    ranked = sorted(scores.items(), key=lambda item: (-item[1], item[0]))
    label, top = ranked[0]
    runner, second = ranked[1] if len(ranked) > 1 else ("", 0.0)
    margin = top - second
    n_detected = detected.get(label, 0)
    if top < 0.1 or n_detected < MIN_MARKERS_DETECTED:
        cell_type, confidence = "Unassigned", "Unassigned"
    elif margin < 0.1:
        cell_type, confidence = "Ambiguous", "Low"
    elif top >= HIGH_SCORE and margin >= HIGH_MARGIN and n_detected >= 3:
        cell_type, confidence = label, "High"
    elif top >= MODERATE_SCORE and margin >= MODERATE_MARGIN and n_detected >= MIN_MARKERS_DETECTED:
        cell_type, confidence = label, "Moderate"
    else:
        cell_type, confidence = label, "Low"
    return {"cell_type": cell_type, "confidence": confidence, "top_scoring_cell_type": label, "top_score": float(top), "runner_up": runner, "margin": float(margin)}


def marker_set_definitions() -> dict:
    """Load the documented marker gene sets shipped with this package."""
    path = Path(__file__).with_name("marker_sets.json")
    return json.loads(path.read_text())


def score_marker_sets(adata, marker_sets: dict) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Score marker sets in adata.raw and return cell annotations and gene-set audit rows."""
    import numpy as np

    raw = adata.raw
    gene_lookup = {str(name).upper(): i for i, name in enumerate(raw.var_names)}
    scores_by_type: dict[str, np.ndarray] = {}
    detected_by_type: dict[str, np.ndarray] = {}
    audit_rows: list[dict[str, object]] = []
    for cell_type, requested_genes in marker_sets["cell_types"].items():
        present = [gene for gene in requested_genes if gene.upper() in gene_lookup]
        indices = [gene_lookup[gene.upper()] for gene in present]
        if indices:
            matrix = raw.X[:, indices]
            scores_by_type[cell_type] = np.asarray(matrix.mean(axis=1)).ravel()
            detected_by_type[cell_type] = np.asarray((matrix > 0).sum(axis=1)).ravel()
        else:
            scores_by_type[cell_type] = np.zeros(adata.n_obs, dtype=float)
            detected_by_type[cell_type] = np.zeros(adata.n_obs, dtype=int)
        audit_rows.append({
            "cell_type": cell_type,
            "marker_genes": ";".join(requested_genes),
            "genes_present_in_dataset": ";".join(present),
            "genes_present_count": len(present),
            "source": marker_sets["source"],
            "source_url": marker_sets["source_url"],
        })
    annotation_rows = []
    for i in range(adata.n_obs):
        scores = {name: float(values[i]) for name, values in scores_by_type.items()}
        detected = {name: int(values[i]) for name, values in detected_by_type.items()}
        result = classify_marker_scores(scores, detected)
        label = result["cell_type"]
        annotation_rows.append({
            "cell_type": label,
            "confidence": result["confidence"],
            "top_marker_score": result["top_score"],
            "runner_up_cell_type": result["runner_up"],
            "score_margin": result["margin"],
            "top_scoring_cell_type": result["top_scoring_cell_type"],
            "top_marker_genes_detected": detected.get(result["top_scoring_cell_type"], 0),
            **{f"score_{name}": float(values[i]) for name, values in scores_by_type.items()},
            **{f"detected_{name}": int(values[i]) for name, values in detected_by_type.items()},
        })
    return {key: [row[key] for row in annotation_rows] for key in annotation_rows[0]}, audit_rows

def cluster_marker_score_summary(obs):
    """Summarize mean marker scores and marker detection by cluster and set."""
    import pandas as pd

    if "leiden" not in obs.columns:
        raise ValueError("obs must contain Leiden cluster labels in 'leiden'")
    score_columns = [column for column in obs.columns if str(column).startswith("score_")]
    if not score_columns:
        raise ValueError("obs must contain per-cell marker score columns")

    rows = []
    for cluster, group in obs.groupby("leiden", observed=True, sort=False):
        for score_column in score_columns:
            cell_type = str(score_column)[len("score_"):]
            detected_column = f"detected_{cell_type}"
            if detected_column not in group.columns:
                raise ValueError(f"Missing per-cell marker evidence column: {detected_column}")
            values = group[score_column].astype(float)
            detected = group[detected_column].astype(int)
            rows.append({
                "cluster": str(cluster),
                "cell_type": cell_type,
                "mean_marker_score": float(values.mean()),
                "pct_cells_with_2plus_markers": float((detected >= MIN_MARKERS_DETECTED).mean() * 100),
                "cells_in_cluster": int(len(group)),
            })
    return pd.DataFrame(rows)


def plot_cluster_marker_scores(summary, output_path: Path) -> None:
    """Draw a dot plot: color is mean score; size is fraction with >=2 markers."""
    import matplotlib.pyplot as plt
    import numpy as np

    clusters = list(dict.fromkeys(summary["cluster"].astype(str)))
    cell_types = list(dict.fromkeys(summary["cell_type"].astype(str)))
    x_lookup = {name: index for index, name in enumerate(cell_types)}
    y_lookup = {name: index for index, name in enumerate(clusters)}
    x = summary["cell_type"].astype(str).map(x_lookup).to_numpy()
    y = summary["cluster"].astype(str).map(y_lookup).to_numpy()
    scores = summary["mean_marker_score"].astype(float).to_numpy()
    fractions = summary["pct_cells_with_2plus_markers"].astype(float).to_numpy()

    fig_width = max(8, len(cell_types) * 1.25)
    fig_height = max(4, len(clusters) * 0.55 + 1.5)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    scatter = ax.scatter(
        x, y, c=scores, s=20 + fractions * 2.2, cmap="viridis",
        edgecolors="#333333", linewidths=0.35, alpha=0.9,
    )
    ax.set_xticks(range(len(cell_types)), cell_types, rotation=35, ha="right")
    ax.set_yticks(range(len(clusters)), clusters)
    ax.set_xlabel("Marker set / candidate cell type")
    ax.set_ylabel("Leiden cluster")
    ax.set_title("Per-cluster marker-score evidence")
    ax.grid(axis="both", color="#dddddd", linewidth=0.5, zorder=0)
    colorbar = fig.colorbar(scatter, ax=ax, pad=0.02)
    colorbar.set_label("Mean log1p-normalized marker score")
    for fraction in (25, 50, 75):
        ax.scatter([], [], s=20 + fraction * 2.2, color="#aaaaaa", edgecolors="#333333", label=f"{fraction}%")
    ax.legend(title="Cells with ≥2 markers", bbox_to_anchor=(1.18, 1), loc="upper left", frameon=False)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def run(args: argparse.Namespace) -> None:
    # Heavy dependencies are imported only when executing the full analysis.
    import matplotlib.pyplot as plt
    import scanpy as sc

    validate_parameters(args.min_genes, args.max_mito_pct, args.n_top_genes, args.resolution)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    matrix_dir = Path(args.input_dir) if args.input_dir else download_data(Path(args.data_dir))
    adata = sc.read_10x_mtx(matrix_dir, var_names="gene_symbols", cache=True)
    adata.var_names_make_unique()
    adata.obs["sample"] = "PBMC3K_healthy_donor"
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    axes[0].hist(adata.obs["n_genes_by_counts"], bins=40, color="#4776a8")
    axes[0].set(title="Genes detected per cell", xlabel="Genes", ylabel="Cells")
    axes[1].hist(adata.obs["pct_counts_mt"], bins=40, color="#b26b4d")
    axes[1].set(title="Mitochondrial counts", xlabel="Percent counts", ylabel="Cells")
    fig.tight_layout()
    fig.savefig(out / "qc_before_filtering.png", dpi=160)
    plt.close(fig)

    adata = adata[adata.obs["n_genes_by_counts"] >= args.min_genes].copy()
    adata = adata[adata.obs["pct_counts_mt"] <= args.max_mito_pct].copy()
    sc.pp.filter_genes(adata, min_cells=3)
    if adata.n_obs < 3 or adata.n_vars < 3:
        raise ValueError("Too few cells or genes remain after QC filtering")
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata.copy()
    n_hvg = min(args.n_top_genes, adata.n_vars)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor="seurat")
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=min(10, adata.n_obs - 1), n_pcs=min(40, adata.obsm["X_pca"].shape[1]))
    sc.tl.umap(adata, random_state=args.seed)
    sc.tl.leiden(adata, resolution=args.resolution, random_state=args.seed, key_added="leiden", flavor="igraph", n_iterations=2, directed=False)
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", use_raw=True)
    marker_sets = marker_set_definitions()
    annotations, marker_audit = score_marker_sets(adata, marker_sets)
    for column, values in annotations.items():
        adata.obs[column] = values
    with (out / "marker_sets.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=marker_audit[0].keys())
        writer.writeheader()
        writer.writerows(marker_audit)
    annotation_columns = ["cell_type", "confidence", "top_scoring_cell_type", "top_marker_score", "runner_up_cell_type", "score_margin", "top_marker_genes_detected"]
    annotation_columns.extend(column for column in adata.obs.columns if str(column).startswith("score_"))
    annotation_table = adata.obs[annotation_columns].copy()
    annotation_table.insert(0, "barcode", annotation_table.index.astype(str))
    annotation_table.to_csv(out / "cell_annotations.csv", index=False)
    marker_score_summary = cluster_marker_score_summary(adata.obs)
    marker_score_summary.to_csv(out / "cluster_marker_scores.csv", index=False)
    plot_cluster_marker_scores(marker_score_summary, out / "marker_score_dotplot.png")
    cluster_annotation = adata.obs.groupby(["leiden", "cell_type", "confidence"], observed=True).size().rename("cells").reset_index()
    cluster_annotation.to_csv(out / "cluster_celltype_summary.csv", index=False)

    sc.settings.figdir = out
    sc.pl.umap(adata, color=["leiden", "cell_type", "confidence"], show=False, save="_annotation.png")
    annotation_figure = out / "umap_annotation.png"
    if not annotation_figure.exists():
        candidates = list(out.glob("umap_annotation*"))
        if candidates:
            candidates[0].replace(annotation_figure)
    sc.pl.umap(adata, color=["leiden", "n_genes_by_counts", "pct_counts_mt"], show=False, save="_overview.png")
    generated_figure = out / "umap_overview.png"
    if not generated_figure.exists():
        candidates = list(out.glob("umap_overview*"))
        if candidates:
            candidates[0].replace(generated_figure)
    marker = sc.get.rank_genes_groups_df(adata, group=None)
    marker.to_csv(out / "cluster_markers.csv", index=False)
    adata.write_h5ad(out / "pbmc3k_processed.h5ad", compression="gzip")

    counts = adata.obs["leiden"].value_counts().sort_index()
    with (out / "cluster_cell_counts.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["cluster", "cells"])
        writer.writerows((str(cluster), int(count)) for cluster, count in counts.items())
    summary = {
        "dataset": "10x Genomics 3k PBMCs from a Healthy Donor",
        "source_url": DATA_URL,
        "cells_after_qc": int(adata.n_obs),
        "genes_after_qc_and_hvg": int(adata.n_vars),
        "clusters": int(adata.obs["leiden"].nunique()),
        "annotation_confidence_counts": {str(key): int(value) for key, value in adata.obs["confidence"].value_counts().items()},
        "annotation_score_method": marker_sets["score_method"],
        "annotation_confidence_thresholds": {"high_score": HIGH_SCORE, "high_margin": HIGH_MARGIN, "moderate_score": MODERATE_SCORE, "moderate_margin": MODERATE_MARGIN, "minimum_detected_markers": MIN_MARKERS_DETECTED},
        "marker_set_source": marker_sets["source_url"],
        "parameters": {"min_genes": args.min_genes, "max_mito_pct": args.max_mito_pct, "n_top_genes": args.n_top_genes, "resolution": args.resolution, "seed": args.seed},
        "interpretation_note": "One healthy donor; clusters and marker genes are exploratory and are not independent biological replicates.",
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (out / "analysis_summary.md").write_text(
        "# Analysis summary\n\n"
        f"- Cells retained after QC: {adata.n_obs}\n"
        f"- Highly variable genes used: {adata.n_vars}\n"
        f"- Leiden clusters: {adata.obs['leiden'].nunique()}\n"
        f"- High-confidence heuristic labels: {(adata.obs['confidence'] == 'High').sum()} cells\n"
        f"- Moderate-confidence heuristic labels: {(adata.obs['confidence'] == 'Moderate').sum()} cells\n"
        f"- Low-confidence or ambiguous labels: {adata.obs['confidence'].isin(['Low']).sum()} cells\n"
        f"- Unassigned cells: {(adata.obs['confidence'] == 'Unassigned').sum()}\n\n"
        f"Marker sets and citations: {marker_sets['source_url']}\n\n"
        "Confidence labels are transparent heuristics based on marker score, separation from the runner-up, and detected markers; they are not probabilities and have not been calibrated against an independent validation set. This is a single healthy donor dataset. Cell types and markers are exploratory and do not support donor-level differential expression or batch-integration claims.\n"
    )
    print(f"Analysis complete. Results: {out.resolve()}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", help="Existing 10x matrix directory; otherwise download the public dataset")
    parser.add_argument("--data-dir", default="data", help="Directory for downloaded input data")
    parser.add_argument("--output-dir", default="results", help="Output directory")
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--max-mito-pct", type=float, default=5.0)
    parser.add_argument("--n-top-genes", type=int, default=2000)
    parser.add_argument("--resolution", type=float, default=0.5)
    parser.add_argument("--seed", type=int, default=0)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    run(args)


if __name__ == "__main__":
    main()
