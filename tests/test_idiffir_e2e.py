import csv
import math
import subprocess
import sys
from pathlib import Path

from iDiffIR.SpliceGrapher.formats.GeneModel import GeneModel
from iDiffIR.SpliceGrapher.shared.GeneModelConverter import makeSpliceGraph

TESTS_DIR = Path(__file__).resolve().parent
if str(TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(TESTS_DIR))

from helpers.idiffir_fixture_builder import build_fixture


ROOT = Path(__file__).resolve().parents[1]
IDIFFIR_SCRIPT = ROOT / "scripts" / "idiffir.py"


def _run_idiffir(model: Path, factor1_bam: Path, factor2_bam: Path, output_dir: Path):
    return subprocess.run(
        [
            sys.executable,
            str(IDIFFIR_SCRIPT),
            str(model),
            str(factor1_bam),
            str(factor2_bam),
            "-n",
            "-p",
            "1",
            "-c",
            "0.1",
            "-o",
            str(output_dir),
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )


def _read_introns_table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def test_idiffir_gff3_end_to_end_outputs_directional_ir_metrics(tmp_path: Path):
    fixture = build_fixture(tmp_path)
    out_dir = tmp_path / "out"

    result = _run_idiffir(fixture.gff3, fixture.factor1_bam, fixture.factor2_bam, out_dir)
    assert result.returncode == 0, result.stdout + result.stderr

    lists_dir = out_dir / "lists"
    introns_path = lists_dir / "allIntrons.txt"
    assert introns_path.exists()
    assert (lists_dir / "allDIRs.txt").exists()
    assert (lists_dir / "allDIRGenes.txt").exists()

    rows = _read_introns_table(introns_path)
    assert rows
    assert {row["geneID"] for row in rows} == {"GENE1", "GENE2"}

    for row in rows:
        for field in ("pValue", "adjPValue"):
            value = float(row[field])
            assert math.isfinite(value)
            assert 0.0 <= value <= 1.0

    gene1_rows = [row for row in rows if row["geneID"] == "GENE1"]
    assert gene1_rows

    row = gene1_rows[0]
    assert float(row["logFoldChange"]) > 0
    assert float(row["IRR_ratio_diff"]) > 0
    assert float(row["IRRratio_1"]) > float(row["IRRratio_2"])


def test_gtf_loader_regression_path_processes_gene_and_builds_graph(tmp_path: Path):
    fixture = build_fixture(tmp_path)
    model = GeneModel(str(fixture.gtf))
    genes = model.getGeneRecords("chr1")

    assert genes
    graph = makeSpliceGraph(genes[0])
    assert len(graph.nodeDict) > 0
