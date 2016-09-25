#!/usr/bin/env python

#
#    iDiffIR- identifying differential intron retention (IR) from RNA-seq
#    Copyright (C) 2015  Michael Hamilton
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""
iDiffIR plotting wrapper
"""
import matplotlib
matplotlib.use('agg')
from iDiffIR.Plot import *
import os, sys, numpy, pysam
from SpliceGrapher.formats.fasta import *
from argparse import ArgumentParser, ArgumentTypeError
from SpliceGrapher.shared.utils      import *
from SpliceGrapher.formats.loader import loadGeneModels
from iDiffIR.IntronModel import *

def fileList( raw ):
    """Converts a **:**-separated string of file paths into a list
    
    Parameters
    ----------
    raw : string containing paths separated by a **:**

    Returns
    -------
    p : list
        List of file paths
        
    """
    return [ r.strip() for r in raw.strip().split(':')]

def parseArgs():
    parser = ArgumentParser(description='Plot highlighted regions between two factors.')
    parser.add_argument('-v', '--verbose',dest='verbose', action='store_true', 
                        default=False, help="verbose output [default is quiet running]")

    parser.add_argument('-l', '--factorlabel', dest="factorlabels",action='store', nargs=2,
                        type=str, default=['factor1','factor2'], 
                        help="factor labels, example:  -f Mutant Wildtype", metavar='FACTORLABEL')
    parser.add_argument('-o', '--output-dir',dest='outdir', action='store', 
                        default='iDiffIR_plots', help="output file directory name")
    parser.add_argument('-s', '--shrink_introns', dest='shrink_introns', action='store_true', 
                        default=False, help='shrink introns for depth plots [default is no shrinking]')
    parser.add_argument('-g', '--graph-dirs', dest='graphDirs', type=fileList, 
                        help='colon-separated list of directories to recursively search for SpliceGrapher predictions')
    parser.add_argument('-p', '--procs', dest='procs', action='store', default=1, 
                        type=int, help='Number of processing cores to use, [default = 1]')

    parser.add_argument('genemodel', type=str,
                        help="gene model file: NAME.gtf[.gz] | NAME.gff[.gz]")
    parser.add_argument('genelist', type=str,
                        help="File containing gene records to plot.  Format for lines is geneID\tstart_1,end_1;...start_n,end_n.")

    parser.add_argument('factor1bamfiles', type=fileList,
                        help="colon-separated list of bamfiles: PATH-TO-REPLICATE_1[:PATH-TO-REPLICATE_2,...]")
    parser.add_argument('factor2bamfiles', type=fileList,
                        help="colon-separated list of bamfiles: PATH-TO-REPLICATE_1[:PATH-TO-REPLICATE_2,...]")

    args = parser.parse_args()
    if not validateArgs( args ):
        raise Exception("Argument Errors: check arguments and usage!")
    return args

def _validateBamfiles(bamfileList):
    """Check if bamfiles exist

    Parameters
    ----------
    bamfileList : List of paths to check

    Returns
    -------
    b : bool
        True if all files exist at given paths

    """
    bamfilesOK = True
    for f in bamfileList:
        if not os.path.exists(f):
            sys.stderr.write('**bamfile %s not found\n' % f )
            countFilesOK = False
    return bamfilesOK
    
def validateArgs( nspace ):
    """Verify program arguments

    Checks all arguments to **idiffir.py** for consistency
    and feasibilty.

    Parameters
    ----------
    nspace : argparse.Namespace object containing **idiffir.py** arguments
    
    .. todo:: add rest of argument checks
 
    Returns
    -------
    b : bool
        True if parameters are feasible

    """
    geneModelOK  = True
        
    # bamfiles
    cfOK1 = _validateBamfiles(nspace.factor1bamfiles)
    cfOK2 = _validateBamfiles(nspace.factor2bamfiles)
    countFilesOK = cfOK1 and cfOK2

    # gene model
    if not os.path.isfile(nspace.genemodel):
        sys.stderr.write('**Genene model file %s not found\n' % nspace.genemodel )
        geneModelOK = False

    
    return countFilesOK and geneModelOK

def makeOutputDir(nspace):
    """Make output directories

    Make output directories and subdirectories for **iDiffIR** output

    Parameters
    ----------
    nspace : argparse.Namespace object with **idiffir.py** parameters

    Returns
    -------
    status : bool
             True if directories were able to be created or already exist.
    """
    try:
        # build main directory
        if not os.path.exists(nspace.outdir):
            os.makedirs(nspace.outdir)

        # build figures subdirectory
        if not os.path.exists(os.path.join(nspace.outdir, 'figures')) and not nspace.noplot:
            os.makedirs(os.path.join(nspace.outdir, 'figures'))

        # build figuresLog subdirectory
        if not os.path.exists(os.path.join(nspace.outdir, 'figuresLog')) and not nspace.noplot:
            os.makedirs(os.path.join(nspace.outdir, 'figuresLog'))
        
        # HTML output not implemented yet
        # # build HTML directory
        # if not os.path.exists(os.path.join(nspace.outdir, 'HTML')):
        #     os.makedirs(os.path.join(nspace.outdir, 'HTML'))

        # build lists subdirectory
        if not os.path.exists(os.path.join(nspace.outdir, 'lists')):
            os.makedirs(os.path.join(nspace.outdir, 'lists'))

    # couldn't make a directory...
    except OSError:
        return False
    # all ok
    return True

def writeStatus( status, verbose=False ):
    """Print status message
    
    Pretty-print status message to stdout

    Parameters
    ----------
    status  : string to write
    verbose : whether to display status

    """
    # do nothing if we're keeping quiet
    if not verbose: return

    # print otherwise
    n = len(status)+2
    sys.stderr.write( '%s\n' % ( '-' * n ) )
    sys.stderr.write( ' %s \n' % ( status ) )
    sys.stderr.write( '%s\n' % ( '-' * n ) )

def runIntron(geneRecords, geneModel, nspace, validChroms, f1LNorm, f2LNorm):
    """Run differential IR analysis

    Run **iDiffIR** differential intron retention analysis.

    Parameters
    ----------
    geneRecords : list of IntronModel.IntronModel objects
    geneModel   : SpliceGrapher.formats.GeneModel.GeneModel
                  Gene model object for the provided genome
    nspace      : argparse.ArgumentParser 
                  Command line arguments for **idiffir.py**
    
    """
    # depricated bias adjustment
    #    if nspace.biasAdj:
    #        f1Codes, f2Codes, f1Cnt, f2Cnt = removePositionalBias(
    #            geneRecords, f1Dict, f1Dict, nspace.numClusts)
    #        rmi = rmsip(geneRecords, f1Dict, f2Dict )
    #        plotGBCs( f1Codes, f2Codes, f1Cnt, f2Cnt, rmi, nspace.outdir)
    aVals = range(nspace.krange[0], nspace.krange[1]+1)
    # compute differential IR statistics
    writeStatus('Computing statistics', nspace.verbose)
    computeIRStatistics( geneRecords, nspace, validChroms, f1LNorm, f2LNorm )

    # create a summary dictionary of results
    summaryDict = summary( geneRecords, aVals, nspace.fdrlevel)
    # write full latex table
    fullTexTable(summaryDict,os.path.join(nspace.outdir, 'lists')) 
    # write gene lists
    writeLists( summaryDict, os.path.join(nspace.outdir, 'lists'))
    # write all IR events
    writeAll( geneRecords, aVals, os.path.join(nspace.outdir, 'lists'))
    #writeGeneExpression(geneRecords, os.path.join(nspace.outdir, 'lists'))

    # plot figures for significant events
    writeStatus('Plotting Depths', nspace.verbose)
    f1labs = [ '%s Rep %d' % (nspace.factorlabels[0], i+1) for i in xrange( len(nspace.factor1bamfiles))]
    f2labs = [ '%s Rep %d' % (nspace.factorlabels[1], i+1) for i in xrange( len(nspace.factor2bamfiles))]

    # finish up if we're not plotting
    if nspace.noplot: return

    # plot diagnostic figures (p-value distribution and MvA )
    plotPDist(geneRecords, os.path.join(nspace.outdir, 'figures'))
    plotMVA(geneRecords, aVals, os.path.join(nspace.outdir, 'figures'))

    plotResults( geneRecords, f1labs+f2labs, 
                 nspace, geneModel, False,
                 os.path.join(nspace.outdir, 'figures'))
    plotResults( geneRecords, f1labs+f2labs, 
                 nspace, geneModel, True,
                 os.path.join(nspace.outdir, 'figuresLog'))

    
def runExon(geneRecords, geneModel, nspace, validChroms, f1LNorm, f2LNorm):
    """Run differential exon skipping analysis

    Run iDiffIR's differential exon skipping analysis

    Parameters
    ----------
    geneRecords : list
                  List of IntronModel.IntronModel objects
    geneModel   : SpliceGrapher.formats.GeneModel.GeneModel
                  Gene model object for the provided genome
    nspace      : argparse.ArgumentParser 
                  Command line arguments for **idiffir.py**

    """
    aVals = range(nspace.krange[0], nspace.krange[1]+1)
    writeStatus("Computing SE statistics", nspace.verbose)
    computeSEStatistics( geneRecords, nspace, validChroms, f1LNorm, f2LNorm )
    summaryDictSE = summarySE( geneRecords, aVals, nspace.fdrlevel)
    fullTexTableSE(summaryDictSE,os.path.join(nspace.outdir, 'lists')) 
    #writeListsSE( summaryDict, os.path.join(nspace.outdir, 'lists'))
    writeAllSE( geneRecords, aVals, os.path.join(nspace.outdir, 'lists'))
    writeStatus("Plotting Depths", nspace.verbose)
    f1labs = [ '%s Rep %d' % (nspace.factorlabels[0], i+1) \
               for i in xrange( len(nspace.factor1bamfiles))]
    f2labs = [ '%s Rep %d' % (nspace.factorlabels[1], i+1) \
               for i in xrange( len(nspace.factor2bamfiles))]

    if nspace.noplot: return
    plotPDistSE(geneRecords, os.path.join(nspace.outdir, 'figures'))
    plotMVASE(geneRecords, aVals, os.path.join(nspace.outdir, 'figures'))
    plotResultsSE( geneRecords, f1labs+f2labs, 
                 nspace, geneModel, False,
                 os.path.join(nspace.outdir, 'figures'))
    plotResultsSE( geneRecords, f1labs+f2labs, 
                 nspace, geneModel, True,
                 os.path.join(nspace.outdir, 'figuresLog'))


def _dispatch(event):
    """Run differential event analysis for given event

    Dispatch differential analysis to corresponding function

    Parameters
    ----------
    event : str
            Either IR or SE indicating intron retention or exon skipping
    """
    # run differential IR analysis
    if event == 'IR':
        return runIntron

    # run differential SE analysis
    elif event == 'SE':
        return runExon

    # invalalid event ID
    else:
        raise ValueError

def getValidChromosomes(nspace):
    """Find chromosomes that are covered in 
       each bamfile

    """
    f1Chroms = [ ]
    for path in nspace.factor1bamfiles:
        bamfile = pysam.Samfile(path, 'rb')
        f1Chroms.append(set([chrom.lower() for chrom in bamfile.references]))
    f2Chroms = [ ]
    for path in nspace.factor2bamfiles:
        bamfile = pysam.Samfile(path, 'rb')
        f2Chroms.append(set([chrom.lower() for chrom in bamfile.references]))
    f1Chroms = set.intersection( *f1Chroms )
    f2Chroms = set.intersection( *f2Chroms )
    return set.intersection(f1Chroms, f2Chroms)

def getNormFactors(nspace):
    f1Cnts = [ ]
    try:
        for path in nspace.factor1bamfiles:
            bamfile = pysam.Samfile(path, 'rb')
            f1Cnts.append(bamfile.mapped)
        f2Cnts = [ ]
        for path in nspace.factor2bamfiles:
            bamfile = pysam.Samfile(path, 'rb')
            f2Cnts.append(bamfile.mapped)
        mVal = max(f1Cnts+f2Cnts)
        f1LNorm = [mVal/float(x) for x in f1Cnts]
        f2LNorm = [mVal/float(x) for x in f2Cnts]
    except ValueError:
        seps = '*'*41
        sys.exit('%s\n Make sure bamfiles are indexed--exiting\n%s' %
                 (seps,seps))
    return f1LNorm, f2LNorm

def writeStatus( status, verbose=False ):
    """Print status message
    
    Pretty-print status message to stdout

    Parameters
    ----------
    status  : string to write
    verbose : whether to display status

    """
    # do nothing if we're keeping quiet
    if not verbose: return

    # print otherwise
    n = len(status)+2
    sys.stderr.write( '%s\n' % ( '-' * n ) )
    sys.stderr.write( ' %s \n' % ( status ) )
    sys.stderr.write( '%s\n' % ( '-' * n ) )


def main():
    """Main program execution

    """
    nspace = parseArgs()
    validChroms = getValidChromosomes(nspace)
    #f1LNorm, f2LNorm, = getNormFactors(nspace)
    if not os.path.exists(nspace.outdir):
        os.makedirs(nspace.outdir)

    writeStatus('Loading models', nspace.verbose)
    geneModel = loadGeneModels( nspace.genemodel, verbose=nspace.verbose )
    writeStatus('Making reduced models', nspace.verbose)
    geneRecords = makeModels( geneModel, None, verbose=nspace.verbose, 
                              graphDirs=nspace.graphDirs,
                              graphDirsOnly=False,
                              exonic=False, procs=nspace.procs )
    f1labs = [ '%s Rep %d' % (nspace.factorlabels[0], i+1) \
               for i in xrange( len(nspace.factor1bamfiles))]
    f2labs = [ '%s Rep %d' % (nspace.factorlabels[1], i+1) \
               for i in xrange( len(nspace.factor2bamfiles))]
    labs = f1labs + f2labs
    # get genes to be plotted and regions to highlight
    geneDict = {}
    with open(nspace.genelist, 'r') as fin:
        for line in fin:
            geneID, rawlist = line.strip().split()
            highlights = []
            for pair in rawlist.split(';'):
                start,end = pair.split(',')
                highlights.append(sorted((int(start), int(end))))
            geneDict[geneID.lower()] = highlights
    print geneDict
    plotList(geneRecords, geneDict, labs, nspace, geneModel, useLog=True, odir=nspace.outdir)
if __name__ == '__main__':
    main()
