from pathlib import Path


def test_kang_full_rerun_command_is_documented_and_smoke_tested():
    root = Path(__file__).resolve().parents[1]
    readme = (root / "studies/kang_ifnb/README.md").read_text()
    root_readme = (root / "README.md").read_text()
    workflow = (root / ".github/workflows/paired-study.yml").read_text()
    runner = root / "studies/kang_ifnb/run_study.R"

    assert runner.is_file()
    assert "Rscript studies/kang_ifnb/run_study.R" in readme
    assert "Rscript studies/kang_ifnb/run_study.R" in root_readme
    assert "--refresh-data" in readme
    assert "Smoke-test documented full-run command" in workflow
    for output in (
        "celltype_response_summary.csv",
        "hallmark_gsea_results.csv",
        "string_network_edges.csv",
        "string_network.png",
        "kang_interferon_report.html",
        "study_manifest.json",
    ):
        assert output in runner.read_text()
