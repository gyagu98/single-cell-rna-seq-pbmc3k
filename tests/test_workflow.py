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


def test_cluster_annotation_summary_includes_labels_and_ambiguous_counts():
    import pandas as pd
    from pbmc3k_workflow.workflow import cluster_annotation_summary

    obs = pd.DataFrame({
        "leiden": ["0", "0", "0", "0", "1", "1"],
        "cell_type": ["B cell", "Ambiguous", "Ambiguous", "Unassigned", "Monocyte", "Monocyte"],
        "confidence": ["High", "Low", "Low", "Unassigned", "Moderate", "High"],
    })
    result = cluster_annotation_summary(obs)

    cluster0 = result[result["cluster"] == "0"]
    ambiguous0 = cluster0[cluster0["cell_type"] == "Ambiguous"]
    b_cell0 = cluster0[cluster0["cell_type"] == "B cell"].iloc[0]
    assert len(cluster0) == 3
    assert ambiguous0["cells"].tolist() == [2]
    assert set(cluster0["ambiguous_cells"]) == {2}
    assert set(cluster0["unassigned_cells"]) == {1}
    assert set(cluster0["cluster_cells"]) == {4}
    assert b_cell0["confidence"] == "High"
    assert b_cell0["cells"] == 1

    cluster1 = result[result["cluster"] == "1"]
    assert len(cluster1) == 2  # Same cell type with different confidence labels stays distinguishable.
    assert set(cluster1["ambiguous_cells"]) == {0}
    assert set(cluster1["cluster_cells"]) == {2}


def test_cluster_annotation_summary_requires_annotation_columns():
    import pandas as pd
    from pbmc3k_workflow.workflow import cluster_annotation_summary

    with pytest.raises(ValueError, match="missing required annotation columns"):
        cluster_annotation_summary(pd.DataFrame({"leiden": ["0"]}))


def test_cluster_marker_summary_includes_most_common_label_and_confidence():
    import pandas as pd
    from pbmc3k_workflow.workflow import cluster_marker_score_summary

    obs = pd.DataFrame({
        "leiden": ["0", "0", "0", "1", "1"],
        "cell_type": ["B cell", "B cell", "Ambiguous", "Monocyte", "Monocyte"],
        "confidence": ["High", "High", "Low", "Moderate", "High"],
        "score_B cell": [0.8, 0.6, 0.1, 0.1, 0.2],
        "detected_B cell": [3, 2, 1, 0, 1],
        "score_Monocyte": [0.1, 0.2, 0.3, 0.9, 0.8],
        "detected_Monocyte": [0, 1, 2, 4, 3],
    })
    result = cluster_marker_score_summary(obs)

    cluster0 = result[result["cluster"] == "0"]
    assert set(cluster0["assigned_cell_type"]) == {"B cell"}
    assert set(cluster0["assigned_confidence"]) == {"High"}
    cluster1 = result[result["cluster"] == "1"]
    # The two labels tie by count; stable alphabetical tie-break chooses High.
    assert set(cluster1["assigned_cell_type"]) == {"Monocyte"}
    assert set(cluster1["assigned_confidence"]) == {"High"}


def test_cluster_annotation_summary_calculates_ambiguity_warning():
    import pandas as pd
    from pbmc3k_workflow.workflow import cluster_annotation_summary

    obs = pd.DataFrame({
        "leiden": ["0", "0", "0", "0", "1", "1", "1", "1", "1"],
        "cell_type": ["B cell", "Ambiguous", "Ambiguous", "Unassigned", "Monocyte", "Monocyte", "Monocyte", "Monocyte", "Ambiguous"],
        "confidence": ["High", "Low", "Low", "Unassigned", "High", "High", "High", "High", "Low"],
    })
    result = cluster_annotation_summary(obs, ambiguity_threshold=0.20)

    # Cluster 0: 3 of 4 cells ambiguous or unassigned (75.0%) -> True
    cluster0 = result[result["cluster"] == "0"]
    assert set(cluster0["pct_ambiguous_or_unassigned"]) == {75.0}
    assert set(cluster0["high_ambiguity_warning"]) == {True}

    # Cluster 1: 1 of 5 cells ambiguous (20.0%) -> True (threshold is >= 20%)
    cluster1 = result[result["cluster"] == "1"]
    assert set(cluster1["pct_ambiguous_or_unassigned"]) == {20.0}
    assert set(cluster1["high_ambiguity_warning"]) == {True}


def test_cluster_marker_score_summary_includes_ambiguity_metrics():
    import pandas as pd
    from pbmc3k_workflow.workflow import cluster_marker_score_summary, cluster_annotation_summary

    obs = pd.DataFrame({
        "leiden": ["0", "0", "0", "0", "1", "1", "1", "1", "1"],
        "cell_type": ["B cell", "Ambiguous", "Ambiguous", "Unassigned", "Monocyte", "Monocyte", "Monocyte", "Monocyte", "Monocyte"],
        "confidence": ["High", "Low", "Low", "Unassigned", "High", "High", "High", "High", "High"],
        "score_B cell": [0.8, 0.2, 0.3, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
        "detected_B cell": [3, 1, 1, 0, 0, 0, 0, 0, 0],
    })
    annot = cluster_annotation_summary(obs, ambiguity_threshold=0.20)
    result = cluster_marker_score_summary(obs, cluster_summary=annot)

    c0 = result[result["cluster"] == "0"].iloc[0]
    c1 = result[result["cluster"] == "1"].iloc[0]
    assert c0["high_ambiguity_warning"] is True or c0["high_ambiguity_warning"] == 1
    assert c0["pct_ambiguous_or_unassigned"] == pytest.approx(75.0)
    assert c1["high_ambiguity_warning"] is False or c1["high_ambiguity_warning"] == 0
    assert c1["pct_ambiguous_or_unassigned"] == pytest.approx(0.0)


def test_plot_cluster_marker_scores_renders_with_ambiguity_warnings(tmp_path):
    import pandas as pd
    from pbmc3k_workflow.workflow import plot_cluster_marker_scores

    summary = pd.DataFrame({
        "cluster": ["0", "0", "1", "1"],
        "cell_type": ["B cell", "Monocyte", "B cell", "Monocyte"],
        "mean_marker_score": [0.8, 0.1, 0.2, 0.9],
        "pct_cells_with_2plus_markers": [75.0, 10.0, 20.0, 85.0],
        "assigned_cell_type": ["B cell", "B cell", "Monocyte", "Monocyte"],
        "assigned_confidence": ["High", "High", "Low", "Low"],
    })
    cluster_summary = pd.DataFrame({
        "cluster": ["0", "1"],
        "cluster_cells": [10, 10],
        "ambiguous_cells": [1, 4],
        "unassigned_cells": [0, 1],
        "pct_ambiguous_or_unassigned": [10.0, 50.0],
        "high_ambiguity_warning": [False, True],
    })
    out_file = tmp_path / "marker_dotplot.png"
    plot_cluster_marker_scores(summary, out_file, cluster_summary=cluster_summary)
    assert out_file.exists()
    assert out_file.stat().st_size > 1000
