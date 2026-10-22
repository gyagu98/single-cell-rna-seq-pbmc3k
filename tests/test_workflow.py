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
