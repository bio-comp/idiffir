import csv
import math
import subprocess
import sys
from pathlib import Path

import pytest

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


def _assert_bounded_probability_fields(rows: list[dict[str, str]]) -> None:
    for row in rows:
        for field in ("pValue", "adjPValue"):
            value = float(row[field])
            assert math.isfinite(value)
            assert 0.0 <= value <= 1.0


def _write_strict_gtf(path: Path) -> Path:
    path.write_text(
        'chr1\ttest\texon\t1\t100\t.\t+\t.\tgene_id "gene_gtf1"; transcript_id "gene_gtf1.1";\n'
        'chr1\ttest\texon\t201\t300\t.\t+\t.\tgene_id "gene_gtf1"; transcript_id "gene_gtf1.1";\n'
        'chr1\ttest\texon\t401\t500\t.\t+\t.\tgene_id "gene_gtf2"; transcript_id "gene_gtf2.1";\n'
        'chr1\ttest\texon\t601\t700\t.\t+\t.\tgene_id "gene_gtf2"; transcript_id "gene_gtf2.1";\n',
        encoding="utf-8",
    )
    return path


def _run_gff3_fixture(tmp_path: Path) -> list[dict[str, str]]:
    tmp_path.mkdir(parents=True, exist_ok=True)
    fixture = build_fixture(tmp_path)
    out_dir = tmp_path / "out"
    result = _run_idiffir(fixture.gff3, fixture.factor1_bam, fixture.factor2_bam, out_dir)
    assert result.returncode == 0, result.stdout + result.stderr
    return _read_introns_table(out_dir / "lists" / "allIntrons.txt")


@pytest.mark.integration
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
    _assert_bounded_probability_fields(rows)

    gene1_rows = [row for row in rows if row["geneID"] == "GENE1"]
    assert gene1_rows

    row = gene1_rows[0]
    assert float(row["logFoldChange"]) > 0
    assert float(row["IRR_ratio_diff"]) > 0
    assert float(row["IRRratio_1"]) > float(row["IRRratio_2"])


@pytest.mark.integration
def test_gtf_loader_regression_path_processes_gene_and_builds_graph(tmp_path: Path):
    fixture = build_fixture(tmp_path)
    model = GeneModel(str(fixture.gtf))
    genes = model.getGeneRecords("chr1")

    assert genes
    graph = makeSpliceGraph(genes[0])
    assert len(graph.nodeDict) > 0


@pytest.mark.integration
def test_idiffir_gtf_end_to_end_runs_and_writes_introns_table(tmp_path: Path):
    fixture = build_fixture(tmp_path)
    strict_gtf = _write_strict_gtf(tmp_path / "strict_model.gtf")
    out_dir = tmp_path / "out_gtf"

    result = _run_idiffir(strict_gtf, fixture.factor1_bam, fixture.factor2_bam, out_dir)
    assert result.returncode == 0, result.stdout + result.stderr

    lists_dir = out_dir / "lists"
    introns_path = lists_dir / "allIntrons.txt"
    assert introns_path.exists()
    assert (lists_dir / "allDIRs.txt").exists()
    assert (lists_dir / "allDIRGenes.txt").exists()

    rows = _read_introns_table(introns_path)
    assert rows
    assert {row["geneID"] for row in rows} == {"GENE_GTF1", "GENE_GTF2"}
    _assert_bounded_probability_fields(rows)


@pytest.mark.integration
def test_idiffir_gff3_end_to_end_is_deterministic_for_core_invariants(tmp_path: Path):
    run_one_rows = _run_gff3_fixture(tmp_path / "run_one")
    run_two_rows = _run_gff3_fixture(tmp_path / "run_two")

    assert len(run_one_rows) == len(run_two_rows)
    assert {row["geneID"] for row in run_one_rows} == {row["geneID"] for row in run_two_rows}

    run_one_gene_map = {row["geneID"]: row for row in run_one_rows}
    run_two_gene_map = {row["geneID"]: row for row in run_two_rows}
    for gene_id in run_one_gene_map:
        run_one = run_one_gene_map[gene_id]
        run_two = run_two_gene_map[gene_id]
        assert (float(run_one["logFoldChange"]) > 0) == (float(run_two["logFoldChange"]) > 0)
        assert (float(run_one["IRR_ratio_diff"]) > 0) == (float(run_two["IRR_ratio_diff"]) > 0)
