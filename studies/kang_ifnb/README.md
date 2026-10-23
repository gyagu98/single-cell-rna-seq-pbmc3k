# Donor-paired interferon response: analysis components

Question: how does the PBMC transcriptional response to interferon-beta vary
across cell types and individual donors?

Dataset selected: Kang et al. (2018), GSE96583, available as
`muscData::Kang18_8vs8()` / ExperimentHub EH2259. The documented object has
eight lupus donors measured before and after stimulation. This is not a
lupus-versus-healthy comparison.

Sources: [dataset documentation](https://bioconductor.org/packages/release/data/experiment/vignettes/muscData/inst/doc/muscData.html)
and [original paper](https://doi.org/10.1038/nbt.4042).

## Status

These are tested aggregation components and a synthetic statistical test,
not a completed real-data study. The original PBMC3k workflow is unchanged.
No real-data findings are asserted. Input retrieval, metadata-value inspection,
QC, annotation review, diagnostics, pathway analysis, sensitivity runs and a
rendered real-data report remain to be implemented and executed.

## Components

`paired_pseudobulk()` takes a genes-by-cells raw-count matrix and exactly
aligned cell metadata with donor, condition and cell_type columns. It sums
counts per donor/condition within one cell type. It retains only donors with
both conditions meeting the minimum cell threshold (default 20) and requires
at least three complete pairs. The audit records cell counts and excluded
donor strata for the selected cell type. Missing metadata, fractional counts,
misalignment, insufficient pairs and invalid designs cause errors.

`fit_paired_deseq2()` uses `~ donor + condition` and an explicit comparison
versus reference. Genes must have count >=10 in at least three pseudobulk
samples. Results include unshrunken effect sizes and BH adjusted p-values.
The minimum sample and filtering rules are demonstration choices, not a
guarantee of power or appropriate study design for arbitrary datasets.

## Run tests from the repository root

```bash
Rscript studies/kang_ifnb/test_pseudobulk.R
# Requires DESeq2:
Rscript studies/kang_ifnb/test_deseq2.R
```

## Real-data acceptance criteria

Inspect `ind`, `stim`, `cell` and `multiplets` values in the original object;
do not guess label encodings or invent donors. Record source provenance,
checksum, versions and applicable reuse terms. Verify raw counts, retain
documented singlets, inspect per-sample QC, review inherited cell annotations,
and report attrition. Repeat aggregation at minimum cell thresholds 10, 20
and 50. Show donor-pair completeness, sample diagnostics, donor-level effect
consistency and leave-one-donor-out sensitivity. Pathway analysis must use a
versioned reference and the actual tested-gene background. Different
significance levels across cell types do not prove different treatment effects.
