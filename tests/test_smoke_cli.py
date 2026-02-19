import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = [
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
