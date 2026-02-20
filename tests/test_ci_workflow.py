from pathlib import Path
import tomllib


ROOT = Path(__file__).resolve().parents[1]


def test_ci_workflow_does_not_use_python_m_pip():
    ci_workflow = (ROOT / ".github" / "workflows" / "ci.yml").read_text(
        encoding="utf-8"
    )
    assert "python -m pip" not in ci_workflow


def test_ci_workflow_uses_dependency_groups_not_extras():
    ci_workflow = (ROOT / ".github" / "workflows" / "ci.yml").read_text(
        encoding="utf-8"
    )
    assert "--group" in ci_workflow
    assert "--extra test" not in ci_workflow


def test_ci_workflow_matrix_drops_python39():
    ci_workflow = (ROOT / ".github" / "workflows" / "ci.yml").read_text(
        encoding="utf-8"
    )
    assert '"3.9"' not in ci_workflow
    assert '"3.10"' in ci_workflow


def test_project_requires_python_310_or_newer():
    pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    assert pyproject["project"]["requires-python"] == ">=3.10"
    assert (
        "Programming Language :: Python :: 3.9"
        not in pyproject["project"]["classifiers"]
    )
