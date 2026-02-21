
"""Bamfile fetching and writing functionality

Bamfile I/O operations for **iDiffIR**
"""
import sys
from collections import defaultdict
from enum import IntEnum
from typing import Protocol

import numpy
import pysam


class CigarOp(IntEnum):
    """Pysam CIGAR operation codes with symbol lookup."""

    MATCH = 0
    INSERT = 1
    DELETE = 2
    GAP = 3
    SOFT_CLIP = 4
    HARD_CLIP = 5
    PAD = 6
    EQUAL = 7
    DIFF = 8

    @property
    def symbol(self) -> str:
        """Return the canonical SAM CIGAR symbol for this op."""
        return {
            CigarOp.MATCH: "M",
            CigarOp.INSERT: "I",
            CigarOp.DELETE: "D",
            CigarOp.GAP: "N",
            CigarOp.SOFT_CLIP: "S",
            CigarOp.HARD_CLIP: "H",
            CigarOp.PAD: "P",
            CigarOp.EQUAL: "=",
            CigarOp.DIFF: "X",
        }[self]

    @property
    def consumes_query(self) -> bool:
        """Return whether this CIGAR operation consumes query sequence."""
        return self in {
            CigarOp.MATCH,
            CigarOp.INSERT,
            CigarOp.SOFT_CLIP,
            CigarOp.EQUAL,
            CigarOp.DIFF,
        }

    @property
    def consumes_reference(self) -> bool:
        """Return whether this CIGAR operation consumes reference coordinates."""
        return self in {
            CigarOp.MATCH,
            CigarOp.DELETE,
            CigarOp.GAP,
            CigarOp.EQUAL,
            CigarOp.DIFF,
        }

    @classmethod
    def from_code(cls, code: int) -> "CigarOp | None":
        """Resolve a raw pysam CIGAR opcode to an enum value."""
        try:
            return cls(code)
        except ValueError:
            return None


REFERENCE_MATCH_OPS = {CigarOp.MATCH, CigarOp.EQUAL, CigarOp.DIFF}
JunctionCounts = dict[tuple[int, int], int]
DepthsAndJunctions = tuple[numpy.ndarray, JunctionCounts]
ReplicateDepthsAndJunctions = tuple[
    numpy.ndarray,
    numpy.ndarray,
    list[JunctionCounts],
    list[JunctionCounts],
]


class GeneBounds(Protocol):
    """Minimal shape required for depth extraction from gene-like records."""

    chrom: str | bytes
    minpos: int
    maxpos: int


def _normalize_reference_name(name: str | bytes) -> str:
    """Return a lowercase chromosome name as text."""
    if isinstance(name, bytes):
        return name.decode("utf-8").lower()
    return str(name).lower()


def processRead(
    read: pysam.AlignedSegment,
    depths: numpy.ndarray,
    junctions: JunctionCounts,
    min_idx: int,
    max_idx: int,
) -> None:
    """Get read depths and junctions from read

    Get read depths and junctions for aligned read
    and update current depths and junctions.

    Parameters
    ----------
    read : pysam.AlignedRead
           Read to process
    depths : numpy.ndarray
             Read depth vector to update, 0-based indexing
    junctions : dict
                Splice junction library to update
    minIdx, maxIdx : int
                     The minimum and maximum genomic positions
                     for the read depth vector.  len(depths)
                     must equal maxIdx - minIdx
    """
    del max_idx
    if read.cigartuples is None:
        return

    # current genomic position in read
    pos = read.reference_start
    # iterate through alignment CIGAR fields
    for op_code, op_length in read.cigartuples:
        cigar_op = CigarOp.from_code(op_code)
        if cigar_op is None:
            continue

        # contiguous match in reference,
        # increment all mathed positions in
        # gene depths
        if cigar_op in REFERENCE_MATCH_OPS:
            start = max(0, pos - min_idx)
            end = min(len(depths), max(0, (pos + op_length) - min_idx))
            depths[start:end] += 1
            pos = pos + op_length

        #insertion
        elif cigar_op in {CigarOp.INSERT, CigarOp.SOFT_CLIP, CigarOp.HARD_CLIP, CigarOp.PAD}:
            pass
        #deletion
        elif cigar_op == CigarOp.DELETE:
            adj_position = pos - min_idx
            if adj_position >= 0 and adj_position < len(depths):
                depths[adj_position] += 1
            pos += op_length
        # junction
        elif cigar_op == CigarOp.GAP:
            jct = (pos - min_idx, (pos + op_length) - min_idx)
            junctions[jct] += 1
            pos += op_length


def getDepthsFromBam(
    bamfile: str,
    chrom: str | bytes,
    start: int,
    end: int,
) -> DepthsAndJunctions:
    """Get read depths in bamfile for specific location

    Get read depths from a single bamfile.  Helper function
    for fetching multiple replicates/conditions.

    See Also
    --------
    getDepthsFromBamfiles : Wrapper function

    Parameters
    ----------
    bamfile : pysam.Samfile
              Bamfile to fetch read depths
    chrom : str
            Name of chromosome (region)
    start : int
            Minimum location of range
    end : int
          Maximum location of range

    Returns
    -------
    d : numpy.ndarray
        Read depth vector for given location

    """
    if start > end:
        raise ValueError("start must be <= end")

    with pysam.AlignmentFile(bamfile, "rb") as bam_stream:
        chrom_map = {}
        for ref in bam_stream.references:
            normalized = _normalize_reference_name(ref)
            chrom_map[normalized] = ref
            chrom_map[str(ref)] = ref

        normalized_chrom = _normalize_reference_name(chrom)
        if normalized_chrom not in chrom_map:
            sys.stderr.write(f"Chromosome {chrom} not found in bamfile\n")
            return numpy.empty(0), {}

        depths = numpy.zeros(end - start + 1, int)
        junctions = defaultdict(int)
        fetch_start = max(0, start - 1)
        fetch_end = end
        for read in bam_stream.fetch(chrom_map[normalized_chrom], fetch_start, fetch_end):
            processRead(read, depths, junctions, start - 1, end)
    return depths, junctions

def getDepthsFromBamfiles(
    gene: GeneBounds,
    f1files: list[str],
    f2files: list[str],
) -> ReplicateDepthsAndJunctions:
    """Get read depths from bamfiles

    Wrapper function for getting read depths from bamfiles
    from 2 exprerimental conditions.

    See Also
    --------
    getDepthsFromBam : Helper function

    Parameters
    ----------
    gene : iDiffIR.IntronModel.IntronModel
           Gene for which to compute read depths
    f1files, f2files : list
                       Lists of file paths to bamfiles for
                       factor\\ :math:`_1`,
                       factor\\ :math:`_2`, respectively.
    f2files : list
              List of file paths to bamfiles for
              factor\\ :math:`_1`.

    Returns
    -------
    factor1Depths, factor2Depths : numpy.ndarry
                                   :math:`r \\cross N` read depth vectors where
                                   :math:`r` is the number of replicates for the
                                   associate factor and :math:`N` is the length
                                   of the gene.
    factor1juncs, factor2juncs : dict
                                 Splice junctions counts mapping,
                                 (minpos, maxpos) -> count

    """
    factor1djs = [
        getDepthsFromBam(bamfile, gene.chrom, gene.minpos, gene.maxpos)
        for bamfile in f1files
    ]
    factor2djs = [
        getDepthsFromBam(bamfile, gene.chrom, gene.minpos, gene.maxpos)
        for bamfile in f2files
    ]
    factor1depths = numpy.array([entry[0] for entry in factor1djs])
    factor2depths = numpy.array([entry[0] for entry in factor2djs])

    factor1juncs = [entry[1] for entry in factor1djs]
    factor2juncs = [entry[1] for entry in factor2djs]

    return numpy.array(factor1depths), numpy.array(factor2depths), factor1juncs, factor2juncs
