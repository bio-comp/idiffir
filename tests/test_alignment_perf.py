import statistics
import time
from collections.abc import Callable
from pathlib import Path

import numpy
import pytest

from iDiffIR.BamfileIO import getDepthsFromBam
from tests.helpers.alignment_fixture_builder import build_alignment_fixture
from tests.helpers.legacy_depth_reference import compute_depths_and_junctions


def _median_runtime_seconds(
    func: Callable[[], object],
    *,
    warmup: int = 1,
    repeats: int = 7,
) -> float:
    for _ in range(warmup):
        func()

    samples = []
    for _ in range(repeats):
        start = time.perf_counter()
        func()
        samples.append(time.perf_counter() - start)
    return statistics.median(samples)


@pytest.mark.integration
@pytest.mark.performance
def test_depth_counting_path_meets_relative_performance_gate(tmp_path: Path) -> None:
    # Keep workload large enough that compute dominates setup/IO jitter in CI.
    fixture = build_alignment_fixture(tmp_path, repeat_scale=6000)

    def legacy_runner() -> tuple[numpy.ndarray, dict[tuple[int, int], int]]:
        return compute_depths_and_junctions(
            fixture.bam, fixture.chrom, fixture.region_start, fixture.region_end
        )

    def current_runner() -> tuple[numpy.ndarray, dict[tuple[int, int], int]]:
        return getDepthsFromBam(
            str(fixture.bam), fixture.chrom, fixture.region_start, fixture.region_end
        )

    legacy_depths, legacy_junctions = legacy_runner()
    current_depths, current_junctions = current_runner()

    assert numpy.array_equal(current_depths, legacy_depths)
    assert current_junctions == legacy_junctions

    legacy_median = _median_runtime_seconds(legacy_runner)
    current_median = _median_runtime_seconds(current_runner)

    assert current_median <= legacy_median * 0.5, (
        f"Expected >=2x speedup over legacy baseline; "
        f"current={current_median:.6f}s legacy={legacy_median:.6f}s"
    )
