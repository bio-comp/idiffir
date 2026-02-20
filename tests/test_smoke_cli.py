import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = [
    "iDiffIR/SpliceGrapher/scripts/build_classifiers.py",
    "iDiffIR/SpliceGrapher/scripts/classify_sites.py",
    "iDiffIR/SpliceGrapher/scripts/convert_models.py",
    "iDiffIR/SpliceGrapher/scripts/ests_to_splicegraph.py",
    "iDiffIR/SpliceGrapher/scripts/find_splice_forms.py",
    "iDiffIR/SpliceGrapher/scripts/fix_unresolved.py",
    "iDiffIR/SpliceGrapher/scripts/gene_model_to_splicegraph.py",
    "iDiffIR/SpliceGrapher/scripts/generate_known_junctions.py",
    "iDiffIR/SpliceGrapher/scripts/generate_predicted_junctions.py",
    "iDiffIR/SpliceGrapher/scripts/generate_splice_site_data.py",
    "iDiffIR/SpliceGrapher/scripts/generate_putative_sequences.py",
    "iDiffIR/SpliceGrapher/scripts/generate_roc.py",
    "iDiffIR/SpliceGrapher/scripts/get_good_pairs.py",
    "iDiffIR/SpliceGrapher/scripts/gtf2gff.py",
    "iDiffIR/SpliceGrapher/scripts/isolasso_pipeline.py",
    "iDiffIR/SpliceGrapher/scripts/isolasso_update_graphs.py",
    "iDiffIR/SpliceGrapher/scripts/plotter.py",
    "iDiffIR/SpliceGrapher/scripts/predict_graphs.py",
    "iDiffIR/SpliceGrapher/scripts/predict_splicegraph.py",
    "iDiffIR/SpliceGrapher/scripts/psginfer_pipeline.py",
    "iDiffIR/SpliceGrapher/scripts/psginfer_update_graphs.py",
    "iDiffIR/SpliceGrapher/scripts/realignment_pipeline.py",
    "iDiffIR/SpliceGrapher/scripts/sam_collate.py",
    "iDiffIR/SpliceGrapher/scripts/sam_filter.py",
    "iDiffIR/SpliceGrapher/scripts/sam_split.py",
    "iDiffIR/SpliceGrapher/scripts/sam_to_depths.py",
    "iDiffIR/SpliceGrapher/scripts/select_model_parameters.py",
    "iDiffIR/SpliceGrapher/scripts/splice_junction_pipeline.py",
    "iDiffIR/SpliceGrapher/scripts/splicegraph_statistics.py",
    "iDiffIR/SpliceGrapher/scripts/view_splicegraph_multiplot.py",
    "iDiffIR/SpliceGrapher/scripts/view_splicegraphs.py",
    "scripts/convertSam.py",
    "scripts/getDepths.py",
    "scripts/get_gene_expression.py",
    "scripts/get_intron_expression.py",
    "scripts/idiffir.py",
    "scripts/idiffir_plotter.py",
    "scripts/make_MISO_AS_GFF.py",
    "scripts/make_MISO_IR_GFF.py",
    "scripts/make_MISO_SE_GFF.py",
    "scripts/run_miso_ir.py",
    "scripts/simulate_IR.py",
]


def run_help(script_path: str, cwd: Path | None = None):
    return subprocess.run(
        [sys.executable, str(ROOT / script_path), "--help"],
        cwd=cwd,
        capture_output=True,
        text=True,
        check=False,
    )


@pytest.mark.parametrize("script_path", SCRIPTS)
def test_script_help_runs(script_path):
    result = run_help(script_path)
    assert result.returncode == 0, result.stderr
    assert "usage" in result.stdout.lower()


def test_make_miso_as_help_has_no_side_effect_files(tmp_path: Path):
    result = run_help("scripts/make_MISO_AS_GFF.py", cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert not (tmp_path / "a3Maps.txt").exists()
    assert not (tmp_path / "a5Maps.txt").exists()
