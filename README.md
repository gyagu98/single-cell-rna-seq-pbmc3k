# PBMC 3k: reproducible single-cell RNA-seq analysis

A small, reproducible Scanpy workflow for quality control, normalization, highly variable gene selection, PCA, neighborhood graph construction, UMAP, Leiden clustering, and cluster-marker ranking on the public 10x Genomics 3k PBMC healthy-donor dataset.

> **Scope:** this dataset contains 2,700 detected cells from one healthy donor. It is a useful workflow demonstration, not a cohort study. It cannot support donor-level inference, multi-donor batch correction, or claims about disease biology. Cluster labels are computational groups and should be interpreted with marker genes and biological context.

## Dataset and references

The workflow downloads the filtered gene-barcode matrix from the 10x Genomics public PBMC 3k sample archive. The dataset page describes PBMCs from one healthy donor and 2,700 detected cells. [10x Genomics dataset](https://www.10xgenomics.com/datasets/3-k-pbm-cs-from-a-healthy-donor-1-standard-1-1-0) · [Download archive](https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz)

The analysis uses Scanpy's AnnData-based single-cell workflow and follows standard Scanpy preprocessing concepts. The official Scanpy documentation also describes this same public dataset. [Scanpy PBMC 3k dataset](https://scanpy.readthedocs.io/en/stable/generated/scanpy.datasets.pbmc3k.html) · [Scanpy clustering tutorial](https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering-2017.html)

## Workflow

```text
10x filtered matrix → QC metrics and filtering → library-size normalization + log1p
                    → highly variable genes → scaling + PCA → neighbors + UMAP
                    → Leiden clusters + Wilcoxon marker ranking → tables and figures
```

## Install and run

Python 3.10 or newer is required. From a clone of this repository:

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\\Scripts\\activate
python -m pip install --upgrade pip
python -m pip install -e '.[test]'
pytest -q
```

Run the complete analysis. It downloads the public matrix on first use:

```bash
pbmc3k --data-dir data --output-dir results
```

Equivalent module command:

```bash
python -m pbmc3k_workflow --min-genes 200 --max-mito-pct 5 \\
  --n-top-genes 2000 --resolution 0.5 --seed 0 \\
  --data-dir data --output-dir results
```

To use a matrix already downloaded and extracted by 10x Cell Ranger, supply its directory (the folder containing `matrix.mtx`, `genes.tsv`, and `barcodes.tsv`):

```bash
pbmc3k --input-dir path/to/filtered_gene_bc_matrices/hg19 --output-dir results
```

## Outputs

- `pbmc3k_processed.h5ad`: processed AnnData object, with filtered raw counts in `layers['counts']`.
- `qc_before_filtering.png`: pre-filter QC overview.
- `umap_overview.png`: UMAP colored by cluster, detected genes, and mitochondrial percentage.
- `cluster_markers.csv`: ranked marker statistics by cluster.
- `cluster_cell_counts.csv`: number of cells assigned to each cluster.
- `summary.json` and `analysis_summary.md`: run parameters, dimensions, and interpretation caveats.

Run artifacts and downloaded data are excluded from version control. The script uses a fixed random seed by default, and records run parameters in the summary.

## Reproducibility and interpretation

The unit tests check parameter validation and command-line defaults without downloading data or running the full analysis. GitHub Actions runs these tests on pushes and pull requests. A complete analysis requires the Scanpy stack and network access for the initial dataset download.

QC thresholds and clustering resolution affect the resulting cells and clusters; inspect the QC figure and test sensitivity to reasonable parameter changes. Marker ranking is exploratory, especially because cells from this single donor are not independent biological replicates. Do not interpret cell counts as population prevalence or present cluster markers as validated biomarkers.

## License and data

The workflow code is released under the MIT License. The public matrix remains subject to its source dataset's terms; see the linked 10x Genomics dataset page. No donor-identifying information is included in this repository.
