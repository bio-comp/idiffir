
"""Bamfile fetching and writing functionality

Bamfile I/O operations for **iDiffIR**
"""
import pysam, sys, numpy

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
    if chrom not in bamfile.references:
        sys.stderr.write('Chromosome %s not found in bamfile\n')
        return numpy.empty(0)

    readItr = bamfile.fetch(chrom, start, end)
    
def getDepthsFromBamfiles( start, end, f1files, f2files ):
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
    d1, d2 : numpy.ndarry
             :math:`r \cross N` read depth vectors where
             :math:`r` is the number of replicates for the
             associate factor and :math:`N` is the length
             of the gene.
    
    """
    factor1depths = [getDepthsFromBam(bamfile, start, end) \
                     for bamfile in f1files]
    factor2depths = [getDepthsFromBam(bamfile, start, end) \
                     for bamfile in f2files]

    return numpy.array(factor1depths, int), numpy.array(factor2depths, int)


