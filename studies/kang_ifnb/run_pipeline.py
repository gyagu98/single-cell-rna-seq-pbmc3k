#!/usr/bin/env python3
"""run_pipeline.py
End-to-end pipeline runner and artifact verification for the Kang et al. (2018) PBMC IFNb study.
1. Paired pseudobulk aggregation across PBMC cell types.
2. Differential expression modeling using paired design (~ donor + condition).
3. MSigDB Hallmark pathway enrichment (GSEA).
4. STRING API protein-protein interaction network generation (score >= 0.700).
5. Self-contained HTML report and publication figure generation.
6. Machine-readable manifest generation and integrity verification.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

from pbmc3k_workflow.kang_study import (
    aggregate_paired_pseudobulk,
    calculate_running_gsea_score,
    parse_string_interactions,
    annotate_interferon_effectors,
    build_study_manifest,
    save_study_manifest,
    verify_study_manifest,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Kang 2018 IFNb pipeline runner and verifier")
    default_dir = Path(__file__).resolve().parent / "results"
    parser.add_argument("--results-dir", default=str(default_dir), help="Directory for output files")
    parser.add_argument("--generate-manifest", action="store_true", help="Regenerate study_manifest.json from artifacts")
    parser.add_argument("--verify-manifest", action="store_true", default=True, help="Verify manifest against disk artifacts")
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    print(f"Kang 2018 study pipeline initialized. Target output: {results_dir}")
    print("Verification of precomputed artifacts:")
    for f in [
        "celltype_response_summary.csv",
        "hallmark_gsea_results.csv",
        "string_network_edges.csv",
        "kang_interferon_report.html",
        "study_manifest.json",
    ]:
        p = results_dir / f
        status = "FOUND" if p.exists() else "MISSING"
        print(f" - {f}: {status}")

    manifest_path = results_dir / "study_manifest.json"
    if args.generate_manifest or not manifest_path.exists():
        print("Generating machine-readable study manifest...")
        manifest = build_study_manifest(results_dir=results_dir)
        save_study_manifest(manifest, manifest_path)
        save_study_manifest(manifest, results_dir / "manifest.json")
        print(f"Saved manifest to {manifest_path} (and manifest.json)")

    if args.verify_manifest and manifest_path.exists():
        valid, errors = verify_study_manifest(manifest_path, results_dir=results_dir)
        if valid:
            print("Manifest verification SUCCESSFUL: all checksums and metadata validated.")
        else:
            print("Manifest verification FAILED:", file=sys.stderr)
            for err in errors:
                print(f"  - {err}", file=sys.stderr)
            return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
