from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_ci_workflow_does_not_use_python_m_pip():
    ci_workflow = (ROOT / ".github" / "workflows" / "ci.yml").read_text(
        encoding="utf-8"
    )
    assert "python -m pip" not in ci_workflow

