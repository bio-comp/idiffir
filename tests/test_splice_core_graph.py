from __future__ import annotations

from pathlib import Path

from iDiffIR.splice_core.events import (
    AlternativeSpliceSiteType,
    detect_alternative_splice_sites,
    detect_exon_skipping,
    detect_intron_retention,
)
from iDiffIR.splice_core.graph import build_splice_graph
from tests.helpers.splice_graph_fixture_builder import build_splice_graph_fixture


def test_build_splice_graph_detects_ir_es_and_alt_acceptor_events(tmp_path: Path) -> None:
    fixture = build_splice_graph_fixture(tmp_path)
    context = build_splice_graph(
        annotation_path=fixture.gtf_path,
        alignment_path=fixture.bam_path,
    )

    intron_retention_events = detect_intron_retention(context)
    assert intron_retention_events, "expected at least one intron-retention event"
    first_ir = intron_retention_events[0]
    assert first_ir.gene_id == fixture.gene_id
    assert first_ir.retained_support >= 1
    assert first_ir.spliced_support >= 1
    assert first_ir.retained_ratio > 0.0

    exon_skipping_events = detect_exon_skipping(context)
    assert exon_skipping_events, "expected at least one exon-skipping event"
    assert any(
        event.skipped_exon_start == 200 and event.skipped_exon_end == 249
        for event in exon_skipping_events
    )

    alt_events = detect_alternative_splice_sites(context)
    assert alt_events, "expected at least one alternative splice-site event"
    assert any(event.site_type == AlternativeSpliceSiteType.ALT3 for event in alt_events)


def test_build_splice_graph_supports_sam_input(tmp_path: Path) -> None:
    fixture = build_splice_graph_fixture(tmp_path)
    context = build_splice_graph(
        annotation_path=fixture.gtf_path,
        alignment_path=fixture.sam_path,
    )

    intron_retention_events = detect_intron_retention(context)
    assert intron_retention_events


def test_build_splice_graph_is_deterministic(tmp_path: Path) -> None:
    fixture = build_splice_graph_fixture(tmp_path)
    context_one = build_splice_graph(
        annotation_path=fixture.gtf_path,
        alignment_path=fixture.bam_path,
    )
    context_two = build_splice_graph(
        annotation_path=fixture.gtf_path,
        alignment_path=fixture.bam_path,
    )

    signature_one = [
        (
            event.gene_id,
            event.intron_start,
            event.intron_end,
            event.retained_support,
            event.spliced_support,
        )
        for event in detect_intron_retention(context_one)
    ]
    signature_two = [
        (
            event.gene_id,
            event.intron_start,
            event.intron_end,
            event.retained_support,
            event.spliced_support,
        )
        for event in detect_intron_retention(context_two)
    ]
    assert signature_one == signature_two
