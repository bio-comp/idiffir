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
Main iDiffIR script

:author: Mike Hamilton
"""
import os, sys, numpy
from iDiffIR.IntronModel import *
from iDiffIR.Plot import *
from iDiffIR.Stat import *
import matplotlib
matplotlib.use('agg')
from SpliceGrapher.formats.fasta import *
from argparse import ArgumentParser, ArgumentTypeError
from SpliceGrapher.shared.utils      import *
from SpliceGrapher.formats.loader import loadGeneModels


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
    """Parse command line arguments

    Returns
    -------
    a : argparse.ArgumentParser

    """
    parser = ArgumentParser(description='Identify differentially expressed introns.')
    parser.add_argument('-v', '--verbose',dest='verbose', action='store_true', 
                        default=False, help="verbose output [default is quiet running]")
    parser.add_argument('-n', '--noplot', dest='noplot', action='store_true', 
                        default=False, help="Do not plot figures [default is to make figures]")

    parser.add_argument('-l', '--factorlabel', dest="factorlabels",action='store', nargs=2,
                        type=str, default=['factor1','factor2'], 
                        help="factor labels, example:  -f Mutant Wildtype", metavar='FACTORLABEL')
    parser.add_argument('-o', '--output-dir',dest='outdir', action='store', 
                        default='iDiffIR_output', help="output file directory name")
    parser.add_argument('-s', '--shrink_introns', dest='shrink_introns', action='store_true', 
                        default=False, help='shrink introns for depth plots [default is no shrinking]')
    parser.add_argument('-k', '--krange', dest='krange', action='store', 
                        default=[None,None], type=int, help='kmin kmax; [default is to search for kmax]', nargs=2)
    parser.add_argument('-c', '--coverage', dest='coverage', action='store', default=0.99, 
                        type=float, help='coverage cutoff, default = 0.99')
    parser.add_argument('-d', '--dexpThresh', dest='dexpThresh', action='store', default=10, 
                        type=float, help='differential gene expression threshold, [default = 10]')
    parser.add_argument('-p', '--procs', dest='procs', action='store', default=1, 
                        type=int, help='Number of processing cores to use, [default = 1]')

    parser.add_argument('-f', '--fdrlevel', dest='fdrlevel', action='store', default=0.05, 
                        type=float, help='FDR test level, [default = 0.05]')
#    parser.add_argument('-n', '--numClusts', dest='numClusts', action='store', default=5, 
#                        help='number of clusters for generating global bias curves (for addressing read depth bias')
#    parser.add_argument('-d', '--diagnositcs', dest='plotDiagnostics', action='store_true', default=True, 
#                        help='Plot diagnostic figures, default = True')
#    parser.add_argument('-b', '--biasAdj', dest='biasAdj', action='store_true', 
#                        default=False, help='Plot diagnostic figures, default = True')
    parser.add_argument('-g', '--graph-dirs', dest='graphDirs', type=fileList, 
                        help='colon-separated list of directories to recursively search for SpliceGrapher predictions')
    parser.add_argument('-m', '--multTest', dest='multTest', choices=['BF', 'BH', 'QV'], type=str, default='QV',
                        help='Multiple testing adjustment method BF: Bonferroni, BH: Benjamini-Hochberg, QV: q-values [default = QV]')

    parser.add_argument('-e', '--event', choices=['IR', 'SE'], type=str, default='IR',
                        help='AS event to test, IR: Intron Retention, SE: Exon Skipping [default = IR]')
    parser.add_argument('genemodel', type=str,
                        help="gene model file: NAME.gtf[.gz] | NAME.gff[.gz]")
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
            sys.stderr.write('**Counts directory %s not found\n' % f )
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
    cParamOK     = True
        
    # bamfiles
    cfOK1 = _validateBamfiles(nspace.factor1bamfiles)
    cfOK2 = _validateBamfiles(nspace.factor2bamfiles)
    countFilesOK = cfOK1 and cfOK2

    # gene model
    if not os.path.isfile(nspace.genemodel):
        sys.stderr.write('**Genene model file %s not found\n' % nspace.genemodel )
        geneModelOK = False

    # check event coverage parameter
    if nspace.coverage < 0 or nspace.coverage > 1:
        sys.stderr.write( 'Feasible values of coverage filter `-c`, `--coverage`: 0 <= c <= 1\n')
        cParamOK = False

    
    return countFilesOK and cParamOK and geneModelOK

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

def runIntron(geneRecords, geneModel, nspace):
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

    # compute differential IR statistics
    writeStatus('Computing statistics', nspace.verbose)
    testedGenes, aVals = computeStatistics( geneRecords, nspace)

    # create a summary dictionary of results
    summaryDict = summary( testedGenes, aVals, nspace.fdrlevel)
    # write full latex table
    fullTexTable(summaryDict,os.path.join(nspace.outdir, 'lists')) 
    # write gene lists
    writeLists( summaryDict, os.path.join(nspace.outdir, 'lists'))
    # write all IR events
    writeAll( testedGenes, aVals, os.path.join(nspace.outdir, 'lists'))
    #writeGeneExpression(geneRecords, os.path.join(nspace.outdir, 'lists'))

    # plot figures for significant events
    writeStatus('Plotting Depths', nspace.verbose)
    f1labs = [ '%s Rep %d' % (nspace.factorlabels[0], i+1) for i in xrange( len(nspace.factor1Dirs))]
    f2labs = [ '%s Rep %d' % (nspace.factorlabels[1], i+1) for i in xrange( len(nspace.factor2Dirs))]

    # finish up if we're not plotting
    if nspace.noplot: return

    # plot diagnostic figures (p-value distribution and MvA )
    plotPDist(testedGenes, os.path.join(nspace.outdir, 'figures'))
    plotMVA(testedGenes, aVals, os.path.join(nspace.outdir, 'figures'))

    # plot figures in standard scale
    plotResults( testedGenes, aVals, f1Dict, f2Dict, f1labs+f2labs, 
                 nspace, geneModel, False,
                 os.path.join(nspace.outdir, 'figures'))

    # plot figures in log scale
    plotResults( testedGenes, aVals, f1Dict, f2Dict, f1labs+f2labs, 
                 nspace, geneModel, True,
                 os.path.join(nspace.outdir, 'figuresLog'))
    
def runExon(geneRecords, geneModel, nspace):
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
    writeStatus("Computing SE statistics", nspace.verbose)
    SErecords, aVals = computeSEStatistics( geneRecords, nspace )
    plotPDistSE(SErecords, os.path.join(nspace.outdir, 'figures'))
    plotMVASE(SErecords, aVals, os.path.join(nspace.outdir, 'figures'))
    summaryDictSE = summarySE( SErecords, aVals, nspace.fdrlevel)
    fullTexTableSE(summaryDictSE,os.path.join(nspace.outdir, 'lists')) 
    #writeListsSE( summaryDict, os.path.join(nspace.outdir, 'lists'))
    writeAllSE( SErecords, aVals, os.path.join(nspace.outdir, 'lists'))
    writeStatus("Plotting Depths", nspace.verbose)
    f1labs = [ '%s Rep %d' % (nspace.factorlabels[0], i+1) \
               for i in xrange( len(nspace.factor1Dirs))]
    f2labs = [ '%s Rep %d' % (nspace.factorlabels[1], i+1) \
               for i in xrange( len(nspace.factor2Dirs))]

    if nspace.noplot: return
    plotPDistSE(SErecords, os.path.join(nspace.outdir, 'figures'))
    plotMVASE(SErecords, aVals, os.path.join(nspace.outdir, 'figures'))
    plotResultsSE( SErecords, aVals, f1Dict, f2Dict, f1labs+f2labs, 
                 nspace, geneModel, False,
                 os.path.join(nspace.outdir, 'figures'))
    plotResultsSE( SErecords, aVals, f1Dict, f2Dict, f1labs+f2labs, 
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

def main():
    """Main program execution

    """
    nspace = parseArgs()
    if not makeOutputDir(nspace): 
        sys.exit('Could not create directory: %s' % nspace.outdir)

    writeStatus('Loading models', nspace.verbose)
    geneModel = loadGeneModels( nspace.genemodel, verbose=nspace.verbose )

    writeStatus('Making reduced models', nspace.verbose)
    geneRecords = makeModels( geneModel, verbose=nspace.verbose, 
                              graphDirs=nspace.graphDirs, 
                              exonic=nspace.event=='SE', procs=nspace.procs )
    
    # Depricated 
    # writeStatus( 'Loading Depths', nspace.verbose )
    # f1Dict, f2Dict = loadData( nspace, geneModel, geneRecords )
    # genered = [ ]
    # for gene in geneRecords:
    #     if gene.gid in f1Dict and gene.gid in f2Dict and not gene.chrom.isalpha(): genered.append(gene)
    # geneRecords = genered

    _dispatch(nspace.event)(geneRecords, geneModel, f1Dict, f2Dict, nspace)

if __name__ == "__main__":
    main()

