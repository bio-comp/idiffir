import importlib
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_script_module(module_name: str, relative_path: str):
    spec = importlib.util.spec_from_file_location(module_name, ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


def test_core_modules_import():
    modules = [
        "iDiffIR",
        "iDiffIR.IntronModel",
        "iDiffIR.Stat",
        "iDiffIR.Plot",
        "iDiffIR.SpliceGrapher.SpliceGraph",
    ]
    for module_name in modules:
        importlib.import_module(module_name)

    load_script_module("script_idiffir", "scripts/idiffir.py")
    load_script_module("script_idiffir_plotter", "scripts/idiffir_plotter.py")
    load_script_module("script_get_intron_expression", "scripts/get_intron_expression.py")
    load_script_module("script_get_gene_expression", "scripts/get_gene_expression.py")
    load_script_module("script_make_miso_as_gff", "scripts/make_MISO_AS_GFF.py")
    load_script_module("script_make_miso_ir_gff", "scripts/make_MISO_IR_GFF.py")
    load_script_module("script_make_miso_se_gff", "scripts/make_MISO_SE_GFF.py")
    load_script_module("script_simulate_ir", "scripts/simulate_IR.py")
