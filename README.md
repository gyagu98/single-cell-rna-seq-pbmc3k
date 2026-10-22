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
                    → Leiden clusters + marker ranking → heuristic marker-based cell labels
                    → confidence audit, cited marker sets, tables and figures
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
- `cell_annotations.csv`: per-cell broad cell-type annotation, confidence label, top-scoring candidate, top score, runner-up, score margin, and marker detection count.
- `cluster_celltype_summary.csv`: counts by Leiden cluster, cell type, and confidence label, with cluster size, explicit ambiguous- and unassigned-cell totals, ambiguous/unassigned percentage, and high-ambiguity warning flag.
- `marker_sets.csv`: marker genes, genes present in the dataset, and citation URL for each marker set.
- `umap_annotation.png`: UMAP colored by cluster, heuristic annotation, and confidence.
- `marker_score_dotplot.png`: per-cluster marker-score dot plot; each row includes the assigned cell-type/confidence pair, with visual warnings (amber-shaded row highlighting, `[⚠️ % ambig/unassigned]` label tag, and legend status) for clusters where ≥20% of cells are ambiguous or unassigned. Color encodes mean log1p-normalized marker score and dot size encodes the fraction of cells with at least two detected genes in that marker set.
- `cluster_marker_scores.csv`: numeric values plotted, including the cluster’s assigned label and confidence, mean marker score, marker-detection fraction, cluster size, and ambiguity warning metadata. If cells within a cluster have mixed labels, the displayed label is the most frequent cell-type/confidence pair; ties are resolved alphabetically. The complete breakdown remains available in `cluster_celltype_summary.csv`.
- `summary.json` and `analysis_summary.md`: run parameters, dimensions, confidence counts, and interpretation caveats.

Run artifacts and downloaded data are excluded from version control. The script uses a fixed random seed by default, and records run parameters in the summary.

## Marker-based cell-type annotation

After clustering, the workflow scores seven broad PBMC marker sets: monocytes, CD4 T cells, CD8 T cells, natural killer cells, B cells, other T cells, and dendritic cells. The marker lists are taken from the Azimuth Human PBMC reference at the broad `celltype.l1` level: monocytes [CTSS, FCN1, LYZ, S100A9, TYROBP, etc.](https://azimuth.hubmapconsortium.org/references/human_pbmc/), CD4 T cells [IL7R, CD4, TRAC, CD3D, etc.](https://azimuth.hubmapconsortium.org/references/human_pbmc/), CD8 T cells [CD8A, CD8B, CD3D, CTSW, etc.](https://azimuth.hubmapconsortium.org/references/human_pbmc/), NK cells [NKG7, KLRD1, GNLY, PRF1, etc.](https://azimuth.hubmapconsortium.org/references/human_pbmc/), B cells [CD79A, MS4A1, CD74, IGHM, etc.](https://azimuth.hubmapconsortium.org/references/human_pbmc/), other T cells [TRDC, GZMK, KLRB1, TRGC2, etc.](https://azimuth.hubmapconsortium.org/references/human_pbmc/), and dendritic cells [CD74, HLA-DRA, HLA-DPA1, CST3, etc.](https://azimuth.hubmapconsortium.org/references/human_pbmc/). These published reference marker signatures are a starting point for interpretation, not a standalone validation of identities. [Azimuth Human PBMC marker reference](https://azimuth.hubmapconsortium.org/references/human_pbmc/) · [Hao et al., Cell 2021](https://doi.org/10.1016/j.cell.2021.04.048)

For each cell, the workflow calculates the mean log1p-normalized expression across available genes in each set, ranks the scores, and records the top label, runner-up, score margin, and number of detected top-set markers. The `High`, `Moderate`, `Low`, `Ambiguous`, and `Unassigned` labels are heuristic quality flags based on documented score, margin, and marker-detection thresholds in `summary.json`; they are **not probabilities**, calibrated confidence, or an automated reference-mapping result. `marker_sets.csv` records the exact genes and citation for each set, and `cell_annotations.csv` retains per-cell scores and detected-marker counts. `cluster_marker_scores.csv` and `marker_score_dotplot.png` summarize that evidence by cluster so readers can compare candidate labels directly. Dot color shows mean marker score; dot size shows the percentage of cluster cells with at least two markers detected. Clusters with high ambiguous or unassigned cell counts (≥20% based on `cluster_celltype_summary.csv`) are visually flagged with an amber row highlight and a `[⚠️ % ambig/unassigned]` label tag, alerting readers to clusters where cells could not be definitively resolved by marker signatures alone. These are descriptive summaries, not calibrated probabilities. Treat broad labels as tentative, especially where related lymphocyte signatures overlap.

## Reproducibility and interpretation

The unit tests check parameter validation and command-line defaults without downloading data or running the full analysis. GitHub Actions runs these tests on pushes and pull requests. A complete analysis requires the Scanpy stack and network access for the initial dataset download.

QC thresholds and clustering resolution affect the resulting cells and clusters; inspect the QC figure and test sensitivity to reasonable parameter changes. Marker ranking and marker-based cell labels are exploratory, especially because cells from this single donor are not independent biological replicates. The heuristic scores are not probabilities and are not validated against independent labels. Do not interpret cell counts as population prevalence or present cluster markers as validated biomarkers.

## License and data

The workflow code is released under the MIT License. The public matrix remains subject to its source dataset's terms; see the linked 10x Genomics dataset page. No donor-identifying information is included in this repository.
