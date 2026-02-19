import re
import subprocess
import sys
from pathlib import Path

import iDiffIR


ROOT = Path(__file__).resolve().parents[1]


def _pyproject_version() -> str:
    pyproject_text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    match = re.search(r'^version\s*=\s*"([^"]+)"', pyproject_text, re.MULTILINE)
    assert match, "Could not find project version in pyproject.toml"
    return match.group(1)


def test_package_version_matches_pyproject():
    assert iDiffIR.__version__ == _pyproject_version()


def test_compileall_strict_syntaxwarning_gate():
    result = subprocess.run(
        [
            sys.executable,
            "-W",
            "error::SyntaxWarning",
            "-m",
            "compileall",
            "-f",
            "-q",
            "iDiffIR",
            "scripts",
            "tests",
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
