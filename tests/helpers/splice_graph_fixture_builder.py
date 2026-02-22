from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pysam


@dataclass(frozen=True)
class SpliceGraphFixture:
    """Synthetic fixture paths used by splice-core graph and filtering tests."""

    reference_fasta: Path
    gtf_path: Path
    bam_path: Path
    sam_path: Path
    chromosome: str
    gene_id: str


def _build_segment(
    *,
    name: str,
    start_0based: int,
    cigar: tuple[tuple[int, int], ...],
    query_length: int,
) -> pysam.AlignedSegment:
    """Create one aligned segment with deterministic read metadata."""
    segment = pysam.AlignedSegment()
    segment.query_name = name
    segment.query_sequence = "A" * query_length
    segment.flag = 0
    segment.reference_id = 0
    segment.reference_start = start_0based
    segment.mapping_quality = 60
    segment.cigar = cigar
    segment.next_reference_id = -1
    segment.next_reference_start = -1
    segment.template_length = 0
    segment.query_qualities = pysam.qualitystring_to_array("I" * query_length)
    return segment


def build_splice_graph_fixture(tmp_path: Path) -> SpliceGraphFixture:
    """Build a compact BAM/SAM + GTF fixture with IR/ES/A3 evidence."""
    chromosome = "chr1"
    gene_id = "GENE_SYNTH"

    reference_fasta = tmp_path / "reference.fa"
    reference_fasta.write_text(f">{chromosome}\n" + ("A" * 2000) + "\n", encoding="utf-8")
    pysam.faidx(str(reference_fasta))

    # Exons:
    #   E1: 100-149
    #   E2: 200-249
    #   E3: 300-349
    # Transcript T1: E1-E2-E3 (inclusion)
    # Transcript T2: E1-E3    (exon skipping)
    gtf_path = tmp_path / "model.gtf"
    exon_records = [
        (100, 149, "T1"),
        (200, 249, "T1"),
        (300, 349, "T1"),
        (100, 149, "T2"),
        (300, 349, "T2"),
    ]
    gtf_path.write_text(
        "\n".join(
            [
                (
                    f'{chromosome}\tsynth\texon\t{start}\t{end}\t.\t+\t.\t'
                    f'gene_id "{gene_id}"; transcript_id "{transcript_id}";'
                )
                for start, end, transcript_id in exon_records
            ]
            + [""]
        ),
        encoding="utf-8",
    )

    header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": chromosome, "LN": 2000}]}

    segments = [
        # Junction E1->E2: 50M50N50M
        _build_segment(
            name="jct_e1_e2",
            start_0based=99,
            cigar=((0, 50), (3, 50), (0, 50)),
            query_length=100,
        ),
        # Junction E2->E3: 50M50N50M
        _build_segment(
            name="jct_e2_e3",
            start_0based=199,
            cigar=((0, 50), (3, 50), (0, 50)),
            query_length=100,
        ),
        # Skipping junction E1->E3: 50M150N50M
        _build_segment(
            name="jct_e1_e3",
            start_0based=99,
            cigar=((0, 50), (3, 150), (0, 50)),
            query_length=100,
        ),
        # Retained intron read spanning E1+I1+E2: 150M
        _build_segment(
            name="retained_i1",
            start_0based=99,
            cigar=((0, 150),),
            query_length=150,
        ),
    ]

    bam_path = tmp_path / "reads.bam"
    with pysam.AlignmentFile(bam_path, "wb", header=header) as bam_stream:
        for segment in segments:
            bam_stream.write(segment)

    sam_path = tmp_path / "reads.sam"
    with pysam.AlignmentFile(sam_path, "w", header=header) as sam_stream:
        for segment in segments:
            sam_stream.write(segment)

    return SpliceGraphFixture(
        reference_fasta=reference_fasta,
        gtf_path=gtf_path,
        bam_path=bam_path,
        sam_path=sam_path,
        chromosome=chromosome,
        gene_id=gene_id,
    )
