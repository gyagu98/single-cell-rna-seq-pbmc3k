import pytest

from pbmc3k_workflow.workflow import build_parser, validate_parameters


def test_defaults_are_explicit_and_parse():
    args = build_parser().parse_args([])
    assert args.min_genes == 200
    assert args.max_mito_pct == 5.0
    assert args.seed == 0


@pytest.mark.parametrize(
    "kwargs",
    [
        {"min_genes": -1, "max_mito_pct": 5, "n_top_genes": 2000, "resolution": 0.5},
        {"min_genes": 200, "max_mito_pct": 101, "n_top_genes": 2000, "resolution": 0.5},
        {"min_genes": 200, "max_mito_pct": 5, "n_top_genes": 0, "resolution": 0.5},
        {"min_genes": 200, "max_mito_pct": 5, "n_top_genes": 2000, "resolution": 0},
    ],
)
def test_invalid_parameters_are_rejected(kwargs):
    with pytest.raises(ValueError):
        validate_parameters(**kwargs)


def test_valid_parameters_pass():
    validate_parameters(min_genes=200, max_mito_pct=5, n_top_genes=2000, resolution=0.5)


def test_marker_confidence_high_requires_strong_separated_signature():
    from pbmc3k_workflow.workflow import classify_marker_scores

    result = classify_marker_scores({"B cell": 0.9, "Monocyte": 0.3}, {"B cell": 4, "Monocyte": 2})
    assert result["cell_type"] == "B cell"
    assert result["confidence"] == "High"


def test_marker_confidence_flags_ambiguous_top_scores():
    from pbmc3k_workflow.workflow import classify_marker_scores

    result = classify_marker_scores({"CD4 T cell": 0.4, "Other T cell": 0.35}, {"CD4 T cell": 3, "Other T cell": 3})
    assert result["cell_type"] == "Ambiguous"
    assert result["confidence"] == "Low"


def test_marker_confidence_leaves_weak_signatures_unassigned():
    from pbmc3k_workflow.workflow import classify_marker_scores

    result = classify_marker_scores({"Natural killer cell": 0.05, "B cell": 0.01}, {"Natural killer cell": 1, "B cell": 0})
    assert result["cell_type"] == "Unassigned"


def test_cluster_marker_score_summary_aggregates_score_and_detection_fraction():
    import pandas as pd
    from pbmc3k_workflow.workflow import cluster_marker_score_summary

    obs = pd.DataFrame({
        "leiden": ["0", "0", "1"],
        "score_B cell": [0.8, 0.4, 0.1],
        "detected_B cell": [3, 1, 0],
        "score_Monocyte": [0.1, 0.2, 0.9],
        "detected_Monocyte": [0, 2, 4],
    })
    result = cluster_marker_score_summary(obs)
    b_cluster_0 = result[(result["cluster"] == "0") & (result["cell_type"] == "B cell")].iloc[0]
    mono_cluster_1 = result[(result["cluster"] == "1") & (result["cell_type"] == "Monocyte")].iloc[0]
    assert b_cluster_0["mean_marker_score"] == pytest.approx(0.6)
    assert b_cluster_0["pct_cells_with_2plus_markers"] == pytest.approx(50.0)
    assert b_cluster_0["cells_in_cluster"] == 2
    assert mono_cluster_1["pct_cells_with_2plus_markers"] == pytest.approx(100.0)


def test_cluster_marker_score_summary_requires_marker_evidence_columns():
    import pandas as pd
    from pbmc3k_workflow.workflow import cluster_marker_score_summary

    with pytest.raises(ValueError, match="Missing per-cell marker evidence"):
        cluster_marker_score_summary(pd.DataFrame({"leiden": ["0"], "score_B cell": [0.5]}))
