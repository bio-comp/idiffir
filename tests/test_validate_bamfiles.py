from pathlib import Path

import importlib.util


ROOT = Path(__file__).resolve().parents[1]


def load_script_module(module_name: str, relative_path: str):
    spec = importlib.util.spec_from_file_location(module_name, ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


def test_validate_bamfiles_missing_paths_return_false(tmp_path: Path):
    idiffir = load_script_module("script_idiffir", "scripts/idiffir.py")
    idiffir_plotter = load_script_module("script_idiffir_plotter", "scripts/idiffir_plotter.py")
    get_intron_expression = load_script_module(
        "script_get_intron_expression", "scripts/get_intron_expression.py"
    )
    get_gene_expression = load_script_module(
        "script_get_gene_expression", "scripts/get_gene_expression.py"
    )

    missing = str(tmp_path / "does-not-exist.bam")
    validators = [
        idiffir._validateBamfiles,
        idiffir_plotter._validateBamfiles,
        get_intron_expression._validateBamfiles,
        get_gene_expression._validateBamfiles,
    ]
    for validator in validators:
        assert validator([missing]) is False
