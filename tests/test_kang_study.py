"""Unit tests for the Kang 2018 paired pseudobulk, GSEA, and STRING network workflow."""

import pytest
import pandas as pd
import numpy as np

from pbmc3k_workflow.kang_study import (
    aggregate_paired_pseudobulk,
    calculate_running_gsea_score,
    parse_string_interactions,
    annotate_interferon_effectors,
)


def test_aggregate_paired_pseudobulk_valid():
    """Verify donor pairing, count aggregation, and complete pair filtering."""
    np.random.seed(42)
    genes = [f"gene_{i}" for i in range(10)]
    n_cells = 120
    counts = pd.DataFrame(
        np.random.poisson(lam=5, size=(n_cells, len(genes))),
        columns=genes,
        index=[f"cell_{i}" for i in range(n_cells)],
    )

    # 3 donors with 20 cells per condition (3 * 2 * 20 = 120 cells)
    donors = ["d1"] * 40 + ["d2"] * 40 + ["d3"] * 40
    conds = (["ctrl"] * 20 + ["stim"] * 20) * 3
    meta = pd.DataFrame({
        "ind": donors,
        "stim": conds,
        "cell": ["CD14+ Monocytes"] * n_cells,
        "multiplets": ["singlet"] * n_cells,
    }, index=counts.index)

    pb_counts, pb_meta = aggregate_paired_pseudobulk(
        counts, meta, target_cell_type="CD14+ Monocytes", min_cells=20, min_pairs=3
    )

    assert pb_counts.shape == (6, 10)
    assert pb_meta.shape == (6, 4)
    assert set(pb_meta["condition"]) == {"ctrl", "stim"}
    assert set(pb_meta["donor"]) == {"d1", "d2", "d3"}
    # Cell count per sample should be 20
    assert (pb_meta["n_cells"] == 20).all()


def test_aggregate_paired_pseudobulk_incomplete_pairs_filtered():
    """Donors missing a condition or with fewer than min_cells are excluded."""
    genes = ["gene_1", "gene_2"]
    counts = pd.DataFrame(np.ones((50, 2)), columns=genes)
    # donor 1 has both (20, 20), donor 2 has only ctrl (10)
    meta = pd.DataFrame({
        "ind": ["d1"] * 40 + ["d2"] * 10,
        "stim": (["ctrl"] * 20 + ["stim"] * 20) + ["ctrl"] * 10,
        "cell": ["B cells"] * 50,
        "multiplets": ["singlet"] * 50,
    })

    with pytest.raises(ValueError, match="Fewer than 3 complete donor pairs"):
        aggregate_paired_pseudobulk(
            counts, meta, target_cell_type="B cells", min_cells=20, min_pairs=3
        )


def test_calculate_running_gsea_score():
    """Verify running enrichment score calculation and leading edge extraction."""
    genes = [f"gene_{i}" for i in range(50)]
    stats = pd.Series(np.linspace(5, -5, 50), index=genes)

    top_set = set(genes[:10])
    res = calculate_running_gsea_score(stats, top_set)
    assert res["es"] > 0.8
    assert len(res["lead_genes"]) == 10

    bot_set = set(genes[-10:])
    res_bot = calculate_running_gsea_score(stats, bot_set)
    assert res_bot["es"] < -0.8


def test_parse_string_interactions_and_annotation():
    """Verify deduplication, minimum score thresholding, and functional cluster assignment."""
    raw_interactions = [
        {"preferredName_A": "ISG15", "preferredName_B": "USP18", "score": 0.999},
        {"preferredName_A": "USP18", "preferredName_B": "ISG15", "score": 0.999},
        {"preferredName_A": "IFIT1", "preferredName_B": "IFIT2", "score": 0.998},
        {"preferredName_A": "ISG15", "preferredName_B": "IFIT1", "score": 0.850},
        {"preferredName_A": "LOW_CONF", "preferredName_B": "ISG15", "score": 0.400},
    ]

    edges_df, nodes_df = parse_string_interactions(raw_interactions, min_score=0.700)
    assert len(edges_df) == 3
    assert set(nodes_df["gene"]) == {"ISG15", "USP18", "IFIT1", "IFIT2"}

    annotated = annotate_interferon_effectors(nodes_df)
    assert "functional_cluster" in annotated.columns
    isg_cluster = annotated.loc[annotated["gene"] == "ISG15", "functional_cluster"].iloc[0]
    assert isg_cluster == "ISGylation cascade"
    ifit_cluster = annotated.loc[annotated["gene"] == "IFIT1", "functional_cluster"].iloc[0]
    assert ifit_cluster == "Viral RNA recognition (IFIT)"
