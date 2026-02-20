import csv
import math
import subprocess
import sys
from pathlib import Path

import pysam

from iDiffIR.SpliceGrapher.formats.GeneModel import GeneModel
from iDiffIR.SpliceGrapher.shared.GeneModelConverter import makeSpliceGraph


ROOT = Path(__file__).resolve().parents[1]
IDIFFIR_SCRIPT = ROOT / "scripts" / "idiffir.py"


def _write_models(tmp_path: Path) -> tuple[Path, Path]:
    gff3 = tmp_path / "model.gff3"
    gtf = tmp_path / "model.gtf"

    gff3.write_text(
        "##gff-version 3\n"
        "chr1\t.\tgene\t1\t300\t.\t+\t.\tID=gene1;Name=gene1\n"
        "chr1\t.\tmRNA\t1\t300\t.\t+\t.\tID=gene1.1;Parent=gene1\n"
        "chr1\t.\texon\t1\t100\t.\t+\t.\tParent=gene1.1\n"
        "chr1\t.\texon\t201\t300\t.\t+\t.\tParent=gene1.1\n"
        "chr1\t.\tgene\t401\t700\t.\t+\t.\tID=gene2;Name=gene2\n"
        "chr1\t.\tmRNA\t401\t700\t.\t+\t.\tID=gene2.1;Parent=gene2\n"
        "chr1\t.\texon\t401\t500\t.\t+\t.\tParent=gene2.1\n"
        "chr1\t.\texon\t601\t700\t.\t+\t.\tParent=gene2.1\n",
        encoding="utf-8",
    )

    gtf.write_text(
        "chr1\ttest\tgene\t1\t300\t.\t+\t.\tID=gene_gtf;Name=gene_gtf\n"
        "chr1\ttest\tmRNA\t1\t300\t.\t+\t.\tID=gene_gtf.1;Parent=gene_gtf\n"
        "chr1\ttest\texon\t1\t100\t.\t+\t.\tParent=gene_gtf.1\n"
        "chr1\ttest\texon\t201\t300\t.\t+\t.\tParent=gene_gtf.1\n",
        encoding="utf-8",
    )
    return gff3, gtf


def _segment(name: str, start: int, cigar: tuple[tuple[int, int], ...], query_len: int):
    seg = pysam.AlignedSegment()
    seg.query_name = name
    seg.query_sequence = "A" * query_len
    seg.flag = 0
    seg.reference_id = 0
    seg.reference_start = start
    seg.mapping_quality = 60
    seg.cigar = cigar
    seg.next_reference_id = -1
    seg.next_reference_start = -1
    seg.template_length = 0
    seg.query_qualities = pysam.qualitystring_to_array("I" * query_len)
    return seg


def _add_gene_reads(
    reads: list[tuple[int, tuple[tuple[int, int], ...], int]],
    exon1_start: int,
    intron_start: int,
    exon2_start: int,
    intron_reads: int,
    junction_reads: int,
    exon_reads: int,
) -> None:
    reads.extend([(exon1_start, ((0, 50),), 50)] * (exon_reads // 2))
    reads.extend([(exon2_start, ((0, 50),), 50)] * (exon_reads - (exon_reads // 2)))
    reads.extend([(intron_start, ((0, 50),), 50)] * intron_reads)
    reads.extend([(exon1_start + 40, ((0, 40), (3, 100), (0, 40)), 80)] * junction_reads)


def _write_bam(path: Path, gene1_counts: tuple[int, int, int], gene2_counts: tuple[int, int, int]) -> None:
    header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "chr1", "LN": 2000}]}
    reads: list[tuple[int, tuple[tuple[int, int], ...], int]] = []
    _add_gene_reads(reads, 20, 130, 220, *gene1_counts)
    _add_gene_reads(reads, 420, 530, 620, *gene2_counts)
    reads.sort(key=lambda record: record[0])

    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for idx, (start, cigar, query_len) in enumerate(reads):
            out.write(_segment(f"r{idx}", start, cigar, query_len))
    pysam.index(str(path))


def _build_fixture(tmp_path: Path) -> tuple[Path, Path, Path, Path]:
    gff3, gtf = _write_models(tmp_path)
    f1_bam = tmp_path / "f1.bam"
    f2_bam = tmp_path / "f2.bam"
    _write_bam(f1_bam, gene1_counts=(60, 20, 20), gene2_counts=(40, 25, 20))
    _write_bam(f2_bam, gene1_counts=(5, 40, 20), gene2_counts=(10, 35, 20))
    return gff3, gtf, f1_bam, f2_bam


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
    gff3, _, f1_bam, f2_bam = _build_fixture(tmp_path)
    out_dir = tmp_path / "out"

    result = _run_idiffir(gff3, f1_bam, f2_bam, out_dir)
    assert result.returncode == 0, result.stdout + result.stderr

    lists_dir = out_dir / "lists"
    introns_path = lists_dir / "allIntrons.txt"
    assert introns_path.exists()
    assert (lists_dir / "allDIRs.txt").exists()

    rows = _read_introns_table(introns_path)
    assert rows
    gene1_rows = [row for row in rows if row["geneID"] == "GENE1"]
    assert gene1_rows

    row = gene1_rows[0]
    assert float(row["logFoldChange"]) > 0
    assert float(row["IRR_ratio_diff"]) > 0
    for field in ("pValue", "adjPValue"):
        value = float(row[field])
        assert math.isfinite(value)
        assert 0.0 <= value <= 1.0


def test_gtf_loader_regression_path_processes_gene_and_builds_graph(tmp_path: Path):
    _, gtf, _, _ = _build_fixture(tmp_path)
    model = GeneModel(str(gtf))
    genes = model.getGeneRecords("chr1")

    assert genes
    graph = makeSpliceGraph(genes[0])
    assert len(graph.nodeDict) > 0
