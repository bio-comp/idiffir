from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TARGET_SCRIPTS = [
    "scripts/make_MISO_AS_GFF.py",
    "scripts/make_MISO_IR_GFF.py",
    "scripts/make_MISO_SE_GFF.py",
    "scripts/simulate_IR.py",
]


def test_target_scripts_use_argparse_only():
    for script in TARGET_SCRIPTS:
        content = (ROOT / script).read_text(encoding="utf-8")
        assert "optparse" not in content
        assert "argparse" in content
