import sys
from pathlib import Path

import numpy

from iDiffIR.BamfileIO import getDepthsFromBam

TESTS_DIR = Path(__file__).resolve().parent
if str(TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(TESTS_DIR))

from helpers.alignment_fixture_builder import build_alignment_fixture
from helpers.legacy_depth_reference import compute_depths_and_junctions


def _assert_depths_and_junctions_equal(
    actual_depths: numpy.ndarray,
    actual_junctions: dict[tuple[int, int], int],
    expected_depths: numpy.ndarray,
    expected_junctions: dict[tuple[int, int], int],
) -> None:
    assert numpy.array_equal(actual_depths, expected_depths)
    assert actual_junctions == expected_junctions


def test_alignment_fixture_builds_sam_bam_cram_and_annotations(tmp_path: Path):
    fixture = build_alignment_fixture(tmp_path)
    assert fixture.gff3.exists()
    assert fixture.gtf.exists()
    assert fixture.bam.exists()
    assert (fixture.bam.parent / (fixture.bam.name + ".bai")).exists()
    assert fixture.sam.exists()
    assert fixture.cram.exists()
    assert (fixture.cram.parent / (fixture.cram.name + ".crai")).exists()
    assert fixture.reference_fasta.exists()


def test_bam_depth_and_junction_parity_against_legacy_reference(tmp_path: Path):
    fixture = build_alignment_fixture(tmp_path)

    actual_depths, actual_junctions = getDepthsFromBam(
        str(fixture.bam), fixture.chrom, fixture.region_start, fixture.region_end
    )
    expected_depths, expected_junctions = compute_depths_and_junctions(
        fixture.bam, fixture.chrom, fixture.region_start, fixture.region_end
    )
    _assert_depths_and_junctions_equal(actual_depths, actual_junctions, expected_depths, expected_junctions)


def test_cram_reference_depths_available_for_parity_harness(tmp_path: Path):
    fixture = build_alignment_fixture(tmp_path)

    expected_depths, expected_junctions = compute_depths_and_junctions(
        fixture.cram,
        fixture.chrom,
        fixture.region_start,
        fixture.region_end,
        reference_fasta=fixture.reference_fasta,
    )
    assert len(expected_depths) == fixture.region_end - fixture.region_start + 1
    assert expected_junctions


def test_sam_reference_depths_available_for_parity_harness(tmp_path: Path):
    fixture = build_alignment_fixture(tmp_path)

    expected_depths, expected_junctions = compute_depths_and_junctions(
        fixture.sam, fixture.chrom, fixture.region_start, fixture.region_end
    )
    assert len(expected_depths) == fixture.region_end - fixture.region_start + 1
    assert expected_junctions
