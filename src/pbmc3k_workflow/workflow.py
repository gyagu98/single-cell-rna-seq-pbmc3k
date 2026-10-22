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

    sc.settings.figdir = out
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
        "parameters": {"min_genes": args.min_genes, "max_mito_pct": args.max_mito_pct, "n_top_genes": args.n_top_genes, "resolution": args.resolution, "seed": args.seed},
        "interpretation_note": "One healthy donor; clusters and marker genes are exploratory and are not independent biological replicates.",
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (out / "analysis_summary.md").write_text(
        "# Analysis summary\n\n"
        f"- Cells retained after QC: {adata.n_obs}\n"
        f"- Highly variable genes used: {adata.n_vars}\n"
        f"- Leiden clusters: {adata.obs['leiden'].nunique()}\n\n"
        "This is a single healthy donor dataset. Clusters and marker genes are exploratory; this dataset does not support donor-level differential expression or batch-integration claims.\n"
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
