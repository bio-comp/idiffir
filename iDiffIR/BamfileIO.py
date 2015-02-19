
"""Bamfile fetching and writing functionality

Bamfile I/O operations for **iDiffIR**
"""
import pysam, sys, numpy
from collections import defaultdict

MATCH    = 0   #M
INSERT   = 1   #I
DELETE   = 2   #D
GAP      = 3   #N
SOFT_CLIP= 4   #4
HARD_CLIP= 5   #H
PAD      = 6   #P
EQUAL    = 7   #=
DIFF     = 8   #X

def processRead( read, depths, junctions, minIdx, maxIdx ):
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
    # current genomic position in read
    pos = read.pos 
    # iterate through alignment CIGAR fields
    for i,rec in enumerate(read.cigar):
        t,l = rec
        # contiguous match in reference,
        # increment all mathed positions in 
        # gene depths
        if t == MATCH:
            start = max(0, pos-minIdx)
            end   = min(len(depths), max(0, (pos+l)-minIdx))
            depths[start:end] += 1
            pos = pos+ l

        #insertion
        elif t == INSERT:
            pass
        #deletion
        elif t == DELETE:
            adjPosition = pos-minIdx
            if adjPosition >= 0:
                depths[adjPosition] += 1
            pos += l
        # junction
        elif t == GAP:
            jct = (pos-minIdx, (pos+l)-minIdx)
            junctions[jct] += 1
            pos += l

def getDepthsFromBam( bamfile, chrom, start, end ):
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
    bamfile = pysam.Samfile(bamfile, 'rb')
    chromMap = { }
    for ref in bamfile.references:
        chromMap[ref.lower()] = ref
        chromMap[ref] = ref
    if chrom not in chromMap:
        sys.stderr.write('Chromosome %s not found in bamfile\n')
        return numpy.empty(0), {}

    depths = numpy.zeros( end-start+1, int)
    junctions = defaultdict(int)
    readItr = bamfile.fetch(chromMap[chrom], start, end)
    for read in readItr:
        processRead( read, depths, junctions, start-1, end )
    return depths, junctions

def getDepthsFromBamfiles( gene, f1files, f2files ):
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
                       factor\ :math:`_1`, 
                       factor\ :math:`_2`, respectively.
    f2files : list 
              List of file paths to bamfiles for 
              factor\ :math:`_1`.
    
    Returns
    -------
    factor1Depths, factor2Depths : numpy.ndarry
                                   :math:`r \cross N` read depth vectors where
                                   :math:`r` is the number of replicates for the
                                   associate factor and :math:`N` is the length
                                   of the gene.
    factor1juncs, factor2juncs : dict
                                 Splice junctions counts mapping, 
                                 (minpos, maxpos) -> count           
    
    """
    factor1djs = [getDepthsFromBam(bamfile, gene.chrom, 
                                   gene.minpos, gene.maxpos) \
                     for bamfile in f1files]
    factor2djs = [getDepthsFromBam(bamfile, gene.chrom, 
                                   gene.minpos, gene.maxpos) \
                     for bamfile in f2files]
    factor1depths = numpy.array( [ t[0] for t in factor1djs])
    factor2depths = numpy.array( [ t[0] for t in factor2djs])

    factor1juncs = [ t[1] for t in factor1djs]
    factor2juncs = [ t[1] for t in factor2djs]

    return numpy.array(factor1depths), numpy.array(factor2depths), factor1juncs, factor2juncs


    
