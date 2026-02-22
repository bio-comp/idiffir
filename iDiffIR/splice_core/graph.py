"""Build deterministic splice DAGs from annotations and alignment support."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from enum import IntEnum
from pathlib import Path

import networkx as nx
import pysam

from iDiffIR.SpliceGrapher.formats.annotation_io import load_gene_models


class CigarOperation(IntEnum):
    """SAM/BAM CIGAR operations encoded by pysam."""

    MATCH = 0
    INSERTION = 1
    DELETION = 2
    REFERENCE_SKIP = 3
    SOFT_CLIP = 4
    HARD_CLIP = 5
    PADDING = 6
    SEQUENCE_MATCH = 7
    SEQUENCE_MISMATCH = 8


REFERENCE_ADVANCING_OPS = {
    CigarOperation.MATCH,
    CigarOperation.DELETION,
    CigarOperation.REFERENCE_SKIP,
    CigarOperation.SEQUENCE_MATCH,
    CigarOperation.SEQUENCE_MISMATCH,
}

NodeKey = tuple[str, int, int, str]
IntronKey = tuple[str, str, str, int, int]


@dataclass
class IntronSupport:
    """Support metrics for one annotated intron span."""

    chrom: str
    gene_id: str
    strand: str
    intron_start: int
    intron_end: int
    donor_position: int
    acceptor_position: int
    annotation_support: int = 0
    junction_support: int = 0
    retained_support: int = 0


@dataclass(frozen=True)
class SpliceGraphContext:
    """Container for a splice DAG and per-intron support metrics."""

    graph: nx.DiGraph
    introns: tuple[IntronSupport, ...]


def build_splice_graph(
    *,
    annotation_path: str | Path,
    alignment_path: str | Path | None = None,
) -> SpliceGraphContext:
    """Construct a deterministic splice DAG from annotation and optional alignment evidence."""
    annotation_model = load_gene_models(str(annotation_path))
    graph = nx.DiGraph()
    intron_by_key: dict[IntronKey, IntronSupport] = {}
    intron_to_edges: dict[IntronKey, list[tuple[NodeKey, NodeKey]]] = defaultdict(list)

    sorted_genes = sorted(
        annotation_model.getAllGenes(),
        key=lambda item: (item.chromosome, item.minpos, item.id),
    )
    for gene in sorted_genes:
        transcript_items = sorted(gene.isoforms.items(), key=lambda item: item[0])
        for transcript_id, isoform in transcript_items:
            exons = sorted(
                isoform.exons,
                key=lambda exon: (exon.minpos, exon.maxpos),
                reverse=gene.strand == "-",
            )
            if len(exons) < 2:
                continue
            for exon in exons:
                node_key = _node_key(gene.id, exon.minpos, exon.maxpos, gene.strand)
                if node_key not in graph:
                    graph.add_node(
                        node_key,
                        chrom=gene.chromosome,
                        gene_id=gene.id,
                        strand=gene.strand,
                        start=exon.minpos,
                        end=exon.maxpos,
                        transcript_ids={transcript_id},
                    )
                else:
                    graph.nodes[node_key]["transcript_ids"].add(transcript_id)

            for upstream, downstream in zip(exons[:-1], exons[1:]):
                upstream_key = _node_key(gene.id, upstream.minpos, upstream.maxpos, gene.strand)
                downstream_key = _node_key(
                    gene.id,
                    downstream.minpos,
                    downstream.maxpos,
                    gene.strand,
                )
                intron_start = min(upstream.maxpos, downstream.maxpos) + 1
                intron_end = max(upstream.minpos, downstream.minpos) - 1
                if intron_start > intron_end:
                    continue

                donor_position, acceptor_position = _splice_boundaries(
                    strand=gene.strand,
                    intron_start=intron_start,
                    intron_end=intron_end,
                )
                intron_key = _intron_key(
                    chrom=gene.chromosome,
                    gene_id=gene.id,
                    strand=gene.strand,
                    intron_start=intron_start,
                    intron_end=intron_end,
                )
                intron = intron_by_key.get(intron_key)
                if intron is None:
                    intron = IntronSupport(
                        chrom=gene.chromosome,
                        gene_id=gene.id,
                        strand=gene.strand,
                        intron_start=intron_start,
                        intron_end=intron_end,
                        donor_position=donor_position,
                        acceptor_position=acceptor_position,
                    )
                    intron_by_key[intron_key] = intron
                intron.annotation_support += 1
                intron_to_edges[intron_key].append((upstream_key, downstream_key))

                if graph.has_edge(upstream_key, downstream_key):
                    edge = graph[upstream_key][downstream_key]
                    edge["annotation_support"] += 1
                    edge["transcript_ids"].add(transcript_id)
                else:
                    graph.add_edge(
                        upstream_key,
                        downstream_key,
                        intron_start=intron_start,
                        intron_end=intron_end,
                        intron_key=intron_key,
                        annotation_support=1,
                        junction_support=0,
                        retained_support=0,
                        transcript_ids={transcript_id},
                    )

    if alignment_path is not None and intron_by_key:
        _apply_alignment_support(
            alignment_path=Path(alignment_path),
            intron_by_key=intron_by_key,
            intron_to_edges=intron_to_edges,
            graph=graph,
        )

    if not nx.is_directed_acyclic_graph(graph):
        raise ValueError("Splice graph must be acyclic after annotation/alignment ingest.")

    introns = tuple(
        sorted(
            intron_by_key.values(),
            key=lambda intron: (
                intron.chrom,
                intron.gene_id,
                intron.intron_start,
                intron.intron_end,
            ),
        )
    )
    return SpliceGraphContext(graph=graph, introns=introns)


def _apply_alignment_support(
    *,
    alignment_path: Path,
    intron_by_key: dict[IntronKey, IntronSupport],
    intron_to_edges: dict[IntronKey, list[tuple[NodeKey, NodeKey]]],
    graph: nx.DiGraph,
) -> None:
    """Populate junction and retained-read support from one SAM/BAM/CRAM file."""
    introns_by_chrom: dict[str, list[IntronSupport]] = defaultdict(list)
    introns_by_span: dict[tuple[str, int, int], list[IntronKey]] = defaultdict(list)
    for intron_key, intron in intron_by_key.items():
        introns_by_chrom[intron.chrom].append(intron)
        introns_by_span[(intron.chrom, intron.intron_start, intron.intron_end)].append(intron_key)

    for chrom in introns_by_chrom:
        introns_by_chrom[chrom].sort(key=lambda intron: (intron.intron_start, intron.intron_end))

    with pysam.AlignmentFile(
        str(alignment_path),
        _alignment_mode(alignment_path),
    ) as alignment_stream:
        for read in alignment_stream.fetch(until_eof=True):
            if read.is_unmapped or read.reference_id < 0:
                continue
            chrom = alignment_stream.get_reference_name(read.reference_id)
            introns = introns_by_chrom.get(chrom)
            if not introns:
                continue
            _update_junction_support(
                read=read,
                chrom=chrom,
                introns_by_span=introns_by_span,
                intron_by_key=intron_by_key,
            )
            _update_retained_support(read=read, introns=introns)

    for intron_key, intron in intron_by_key.items():
        for upstream_key, downstream_key in intron_to_edges[intron_key]:
            edge = graph[upstream_key][downstream_key]
            edge["junction_support"] = intron.junction_support
            edge["retained_support"] = intron.retained_support


def _update_junction_support(
    *,
    read: pysam.AlignedSegment,
    chrom: str,
    introns_by_span: dict[tuple[str, int, int], list[IntronKey]],
    intron_by_key: dict[IntronKey, IntronSupport],
) -> None:
    """Increment per-intron junction support for reference-skip CIGAR operations."""
    if not read.cigartuples:
        return
    ref_pos = read.reference_start + 1
    for op_code, op_length in read.cigartuples:
        cigar_op = CigarOperation(op_code)
        if cigar_op == CigarOperation.REFERENCE_SKIP:
            intron_start = ref_pos
            intron_end = ref_pos + op_length - 1
            intron_keys = introns_by_span.get((chrom, intron_start, intron_end), ())
            for intron_key in intron_keys:
                intron_by_key[intron_key].junction_support += 1
            ref_pos += op_length
            continue
        if cigar_op in REFERENCE_ADVANCING_OPS:
            ref_pos += op_length


def _update_retained_support(*, read: pysam.AlignedSegment, introns: list[IntronSupport]) -> None:
    """Increment retained support when one aligned block spans an intron interior."""
    blocks = read.get_blocks()
    if not blocks:
        return
    block_ranges = [(block_start + 1, block_end) for block_start, block_end in blocks]
    read_start = block_ranges[0][0]
    read_end = block_ranges[-1][1]

    for intron in introns:
        if intron.intron_start < read_start:
            continue
        if intron.intron_start > read_end:
            break
        if intron.intron_end > read_end:
            continue
        if _spans_intron(
            block_ranges=block_ranges,
            intron_start=intron.intron_start,
            intron_end=intron.intron_end,
        ):
            intron.retained_support += 1


def _spans_intron(
    *,
    block_ranges: list[tuple[int, int]],
    intron_start: int,
    intron_end: int,
) -> bool:
    """Return True when any contiguous alignment block covers the full intron span."""
    for block_start, block_end in block_ranges:
        if block_start <= intron_start and block_end >= intron_end:
            return True
    return False


def _alignment_mode(path: Path) -> str:
    """Select pysam alignment mode from file extension."""
    suffix = path.suffix.lower()
    if suffix == ".sam":
        return "r"
    if suffix == ".cram":
        return "rc"
    return "rb"


def _node_key(gene_id: str, start: int, end: int, strand: str) -> NodeKey:
    """Return deterministic exon node identity."""
    return gene_id, start, end, strand


def _intron_key(
    chrom: str,
    gene_id: str,
    strand: str,
    intron_start: int,
    intron_end: int,
) -> IntronKey:
    """Return deterministic intron identity."""
    return chrom, gene_id, strand, intron_start, intron_end


def _splice_boundaries(*, strand: str, intron_start: int, intron_end: int) -> tuple[int, int]:
    """Return donor/acceptor positions from genomic intron boundaries."""
    if strand == "-":
        return intron_end + 1, intron_start - 1
    return intron_start - 1, intron_end + 1
