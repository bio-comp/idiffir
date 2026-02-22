"""Deterministic traversal queries over splice-core NetworkX graphs."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from enum import Enum

from iDiffIR.splice_core.graph import IntronSupport, NodeKey, SpliceGraphContext


class AlternativeSpliceSiteType(str, Enum):
    """Alternative splice-site classes represented in traversal outputs."""

    ALT5 = "alt5"
    ALT3 = "alt3"


@dataclass(frozen=True)
class IntronRetentionEvent:
    """Intron-retention evidence summary for one intron span."""

    gene_id: str
    chrom: str
    strand: str
    intron_start: int
    intron_end: int
    retained_support: int
    spliced_support: int
    retained_ratio: float


@dataclass(frozen=True)
class ExonSkippingEvent:
    """Exon-skipping evidence summary around one candidate skipped exon."""

    gene_id: str
    chrom: str
    strand: str
    skipped_exon_start: int
    skipped_exon_end: int
    skipped_junction_support: int
    included_junction_support: int


@dataclass(frozen=True)
class AlternativeSpliceSiteEvent:
    """Alternative donor/acceptor event summary."""

    gene_id: str
    chrom: str
    strand: str
    site_type: AlternativeSpliceSiteType
    donor_position: int
    acceptor_position: int
    support: int


def detect_intron_retention(
    context: SpliceGraphContext,
    *,
    min_retained_support: int = 1,
    min_total_support: int = 1,
) -> list[IntronRetentionEvent]:
    """Detect intron-retention events from intron support summaries."""
    events: list[IntronRetentionEvent] = []
    for intron in context.introns:
        if intron.retained_support < min_retained_support:
            continue
        total = intron.retained_support + intron.junction_support
        if total < min_total_support:
            continue
        retained_ratio = float(intron.retained_support) / float(total)
        events.append(
            IntronRetentionEvent(
                gene_id=intron.gene_id,
                chrom=intron.chrom,
                strand=intron.strand,
                intron_start=intron.intron_start,
                intron_end=intron.intron_end,
                retained_support=intron.retained_support,
                spliced_support=intron.junction_support,
                retained_ratio=retained_ratio,
            )
        )
    events.sort(
        key=lambda event: (
            event.chrom,
            event.gene_id,
            event.intron_start,
            event.intron_end,
        )
    )
    return events


def detect_exon_skipping(
    context: SpliceGraphContext,
    *,
    min_skip_support: int = 1,
    min_include_support: int = 1,
) -> list[ExonSkippingEvent]:
    """Detect exon-skipping events using local predecessor/successor topology."""
    graph = context.graph
    signatures: set[tuple[str, str, int, int, NodeKey, NodeKey]] = set()
    events: list[ExonSkippingEvent] = []

    for node_key in graph.nodes:
        node_data = graph.nodes[node_key]
        predecessors = sorted(graph.predecessors(node_key))
        successors = sorted(graph.successors(node_key))
        for predecessor_key in predecessors:
            for successor_key in successors:
                if not graph.has_edge(predecessor_key, successor_key):
                    continue
                skipped_support = int(
                    graph[predecessor_key][successor_key].get("junction_support", 0)
                )
                if skipped_support < min_skip_support:
                    continue
                included_support = int(
                    graph[predecessor_key][node_key].get("junction_support", 0)
                ) + int(
                    graph[node_key][successor_key].get("junction_support", 0)
                )
                if included_support < min_include_support:
                    continue
                signature = (
                    str(node_data["gene_id"]),
                    str(node_data["chrom"]),
                    int(node_data["start"]),
                    int(node_data["end"]),
                    predecessor_key,
                    successor_key,
                )
                if signature in signatures:
                    continue
                signatures.add(signature)
                events.append(
                    ExonSkippingEvent(
                        gene_id=str(node_data["gene_id"]),
                        chrom=str(node_data["chrom"]),
                        strand=str(node_data["strand"]),
                        skipped_exon_start=int(node_data["start"]),
                        skipped_exon_end=int(node_data["end"]),
                        skipped_junction_support=skipped_support,
                        included_junction_support=included_support,
                    )
                )

    events.sort(
        key=lambda event: (
            event.chrom,
            event.gene_id,
            event.skipped_exon_start,
            event.skipped_exon_end,
        )
    )
    return events


def detect_alternative_splice_sites(
    context: SpliceGraphContext,
    *,
    min_support: int = 1,
) -> list[AlternativeSpliceSiteEvent]:
    """Detect ALT5/ALT3 events from per-intron donor/acceptor support groups."""
    introns = [intron for intron in context.introns if intron.junction_support >= min_support]
    events: list[AlternativeSpliceSiteEvent] = []
    event_signatures: set[tuple[str, str, str, str, int, int]] = set()

    donor_groups: dict[tuple[str, str, str, int], list[IntronSupport]] = defaultdict(list)
    acceptor_groups: dict[tuple[str, str, str, int], list[IntronSupport]] = defaultdict(list)
    for intron in introns:
        donor_groups[
            (intron.gene_id, intron.chrom, intron.strand, intron.donor_position)
        ].append(intron)
        acceptor_groups[
            (intron.gene_id, intron.chrom, intron.strand, intron.acceptor_position)
        ].append(intron)

    for (gene_id, chrom, strand, donor_position), grouped_introns in donor_groups.items():
        acceptor_positions = sorted({intron.acceptor_position for intron in grouped_introns})
        if len(acceptor_positions) < 2:
            continue
        for intron in grouped_introns:
            signature = (
                gene_id,
                chrom,
                strand,
                AlternativeSpliceSiteType.ALT3.value,
                donor_position,
                intron.acceptor_position,
            )
            if signature in event_signatures:
                continue
            event_signatures.add(signature)
            events.append(
                AlternativeSpliceSiteEvent(
                    gene_id=gene_id,
                    chrom=chrom,
                    strand=strand,
                    site_type=AlternativeSpliceSiteType.ALT3,
                    donor_position=donor_position,
                    acceptor_position=intron.acceptor_position,
                    support=intron.junction_support,
                )
            )

    for (gene_id, chrom, strand, acceptor_position), grouped_introns in acceptor_groups.items():
        donor_positions = sorted({intron.donor_position for intron in grouped_introns})
        if len(donor_positions) < 2:
            continue
        for intron in grouped_introns:
            signature = (
                gene_id,
                chrom,
                strand,
                AlternativeSpliceSiteType.ALT5.value,
                intron.donor_position,
                acceptor_position,
            )
            if signature in event_signatures:
                continue
            event_signatures.add(signature)
            events.append(
                AlternativeSpliceSiteEvent(
                    gene_id=gene_id,
                    chrom=chrom,
                    strand=strand,
                    site_type=AlternativeSpliceSiteType.ALT5,
                    donor_position=intron.donor_position,
                    acceptor_position=acceptor_position,
                    support=intron.junction_support,
                )
            )

    events.sort(
        key=lambda event: (
            event.chrom,
            event.gene_id,
            event.site_type.value,
            event.donor_position,
            event.acceptor_position,
        )
    )
    return events
