from pathlib import Path

import pysam

from iDiffIR.BamfileIO import CigarOp, getDepthsFromBam
from iDiffIR.SpliceGrapher.formats.alignment_io import getSamReadData


def _make_segment(
    name: str,
    start: int,
    cigar: tuple[tuple[int, int], ...],
    query_len: int,
) -> pysam.AlignedSegment:
    """Build a synthetic aligned segment for coordinate regression tests."""
    segment = pysam.AlignedSegment()
    segment.query_name = name
    segment.query_sequence = "A" * query_len
    segment.flag = 0
    segment.reference_id = 0
    segment.reference_start = start
    segment.mapping_quality = 60
    segment.cigar = cigar
    segment.next_reference_id = -1
    segment.next_reference_start = -1
    segment.template_length = 0
    segment.query_qualities = pysam.qualitystring_to_array("I" * query_len)
    return segment


def _write_bam(tmp_path: Path, segments: list[pysam.AlignedSegment]) -> tuple[Path, str]:
    """Write indexed BAM for the provided segments."""
    chrom = "chr1"
    bam = tmp_path / "reads.bam"
    header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": chrom, "LN": 200}]}
    with pysam.AlignmentFile(bam, "wb", header=header) as stream:
        for segment in segments:
            stream.write(segment)
    pysam.index(str(bam))
    return bam, chrom


def test_get_depths_from_bam_includes_left_boundary_base(tmp_path: Path) -> None:
    """A query covering genomic base 1 should be counted when requesting range [1, 1]."""
    segment = _make_segment("left_boundary", start=0, cigar=((0, 1),), query_len=1)
    bam_path, chrom = _write_bam(tmp_path, [segment])

    depths, junctions = getDepthsFromBam(str(bam_path), chrom, 1, 1)

    assert depths.tolist() == [1]
    assert dict(junctions) == {}


def test_bamfileio_uses_enum_cigar_ops() -> None:
    """CIGAR operation codes should be represented as an enum, not loose constants."""
    assert CigarOp.MATCH.value == 0
    assert CigarOp.INSERT.value == 1
    assert CigarOp.DELETE.value == 2
    assert CigarOp.GAP.value == 3


def test_get_sam_read_data_handles_complex_cigar_reference_positions(tmp_path: Path) -> None:
    """Soft clips/insertions/deletions must not shift depth/junction reference coordinates."""
    # 2S 3M 2I 2M 2D 1M 3N 4M 1H
    segment = _make_segment(
        "complex_cigar",
        start=0,
        cigar=((4, 2), (0, 3), (1, 2), (0, 2), (2, 2), (0, 1), (3, 3), (0, 4), (5, 1)),
        query_len=14,
    )
    bam_path, chrom = _write_bam(tmp_path, [segment])

    depths, junctions = getSamReadData(str(bam_path), chromosomes=[chrom], minjct=1)
    observed = [depths[chrom][position] for position in range(1, 16)]

    assert observed == [1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1]
    assert len(junctions[chrom]) == 1
    assert junctions[chrom][0].donor() == 8
    assert junctions[chrom][0].acceptor() == 12
