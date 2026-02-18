import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def run_help(script_path: str):
    return subprocess.run(
        [sys.executable, str(ROOT / script_path), "--help"],
        capture_output=True,
        text=True,
        check=False,
    )


def test_idiffir_help_runs():
    result = run_help("scripts/idiffir.py")
    assert result.returncode == 0, result.stderr
    assert "usage: idiffir.py" in result.stdout


def test_convert_sam_help_runs():
    result = run_help("scripts/convertSam.py")
    assert result.returncode == 0, result.stderr
    assert "usage: convertSam.py" in result.stdout
