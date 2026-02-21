from pathlib import Path
import re

import pytest


ROOT = Path(__file__).resolve().parents[1]


def _ci_workflow_text() -> str:
    return (ROOT / ".github" / "workflows" / "ci.yml").read_text(encoding="utf-8")


def _pyproject_text() -> str:
    return (ROOT / "pyproject.toml").read_text(encoding="utf-8")


@pytest.mark.unit
def test_ci_workflow_does_not_use_python_m_pip():
    ci_workflow = _ci_workflow_text()
    assert "python -m pip" not in ci_workflow


@pytest.mark.unit
def test_ci_workflow_uses_dependency_groups_not_extras():
    ci_workflow = _ci_workflow_text()
    assert "--group" in ci_workflow
    assert "--extra test" not in ci_workflow


@pytest.mark.unit
def test_ci_workflow_matrix_drops_python39():
    ci_workflow = _ci_workflow_text()
    assert '"3.9"' not in ci_workflow
    assert '"3.10"' in ci_workflow


@pytest.mark.unit
def test_project_requires_python_310_or_newer():
    pyproject_text = _pyproject_text()
    assert re.search(
        r'^\s*requires-python\s*=\s*">=3\.10"\s*$',
        pyproject_text,
        re.MULTILINE,
    )
    assert '"Programming Language :: Python :: 3.9"' not in pyproject_text


@pytest.mark.unit
def test_ci_workflow_skips_performance_tests() -> None:
    ci_workflow = _ci_workflow_text()
    assert '-m "not performance"' in ci_workflow
