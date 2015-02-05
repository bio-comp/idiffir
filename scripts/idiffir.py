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
    return [ r.strip() for r in raw.strip().split(':')]

def parseArgs():
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

def validateArgs( nspace ):
    """
    Verify arguments
    """
    countFilesOK = True
    geneModelOK  = True
    cParamOK     = True
        
    # count files
    for f in nspace.factor1Dirs:
        if not os.path.exists(f):
            sys.stderr.write('**Counts directory %s not found\n' % f )
            countFilesOK = False

    for f in nspace.factor2Dirs:
        if not os.path.exists(f):
            sys.stderr.write('**Counts directory %s not found\n' % f )
            countFilesOK = False

    # gene model
    if not os.path.isfile(nspace.genemodel):
        sys.stderr.write('**Genene model file %s not found\n' % nspace.genemodel )
        geneModelOK = False

    if nspace.coverage < 0 or nspace.coverage > 1:
        sys.stderr.write( 'Feasible values of coverage filter `-c`, `--coverage`: 0 <= c <= 1\n')
        cParamOK = False

    return countFilesOK and cParamOK and geneModelOK

def makeOutputDir(nspace):
    """
    Make directories and subdirectories for output.
    """
    try:
        if not os.path.exists(nspace.outdir):
            os.makedirs(nspace.outdir)
        
        if not os.path.exists(os.path.join(nspace.outdir, 'figures')) and not nspace.noplot:
            os.makedirs(os.path.join(nspace.outdir, 'figures'))
        if not os.path.exists(os.path.join(nspace.outdir, 'figuresLog')) and not nspace.noplot:
            os.makedirs(os.path.join(nspace.outdir, 'figuresLog'))

#        if not os.path.exists(os.path.join(nspace.outdir, 'HTML')):
#            os.makedirs(os.path.join(nspace.outdir, 'HTML'))
        if not os.path.exists(os.path.join(nspace.outdir, 'lists')):
            os.makedirs(os.path.join(nspace.outdir, 'lists'))

    except OSError:
        return False
    return True

def writeStatus( status, verbose=False ):
    """
    Pretty print status
    """
    if not verbose: return
    n = len(status)+2
    sys.stderr.write( '%s\n' % ( '-' * n ) )
    sys.stderr.write( ' %s \n' % ( status ) )
    sys.stderr.write( '%s\n' % ( '-' * n ) )

def loadData( nspace, geneModel, geneRecords ):
    """
    Load gene depths
    """
    f1Dict = { }
    f2Dict = { }
    grdict = { }
    for gene in geneRecords:
        grdict[gene.gid] = gene
    def load( factorfiles, fDict ):
        for chrm in sorted(geneModel.getChromosomes()):
            genes     = geneModel.getGeneRecords(chrm)
            genes.sort()

            for i in xrange(len(factorfiles)):
                fname = os.path.join(factorfiles[i], 
                                     '%s.cnt.npy' % chrm) 
                if not os.path.exists(fname):
                    sys.stderr.write('Depths file %s not found\n' % fname )
                    continue
                depths = numpy.load(fname)
                if nspace.verbose:
                    sys.stderr.write("Reading depths from %s\n" % ( fname ) )
                depths = csr_matrix( depths )
                counter = 0
                for gene in genes:
                    grec = grdict[gene.id]
                    start = grec.minpos - 1
                    end = grec.maxpos
                    l = fDict.get(gene.id, [ ] )
                    l.append( depths[0,start:end])
                    fDict[gene.id] = l
                    counter += 1
                if nspace.verbose:
                    sys.stderr.write("Read %d gene depths\n" % (counter))
    if nspace.verbose:
        sys.stderr.write( '-' * 30 + '\n' )
        sys.stderr.write("|Reading factor 1 gene depths|\n")
        sys.stderr.write( '-' * 30 + '\n' )

    load( nspace.factor1Dirs, f1Dict)
    if nspace.verbose:
        sys.stderr.write( '-' * 30 + '\n' )
        sys.stderr.write("|Reading factor 2 gene depths|\n")
        sys.stderr.write( '-' * 30 + '\n' )
    load( nspace.factor2Dirs, f2Dict)
    return f1Dict, f2Dict 

def runIntron(geneRecords, geneModel, f1Dict, f2Dict, nspace):
    """
    Run differential IR analysis
    """
    writeStatus('Computing statistics', nspace.verbose)
#    if nspace.biasAdj:
#        f1Codes, f2Codes, f1Cnt, f2Cnt = removePositionalBias(
#            geneRecords, f1Dict, f1Dict, nspace.numClusts)
#        rmi = rmsip(geneRecords, f1Dict, f2Dict )
#        plotGBCs( f1Codes, f2Codes, f1Cnt, f2Cnt, rmi, nspace.outdir)
    testedGenes, aVals = computeStatistics( geneRecords, f1Dict, f2Dict, 
                                      nspace )
    summaryDict = summary( testedGenes, aVals, nspace.fdrlevel)
    fullTexTable(summaryDict,os.path.join(nspace.outdir, 'lists')) 
    writeLists( summaryDict, os.path.join(nspace.outdir, 'lists'))
    writeAll( testedGenes, aVals, os.path.join(nspace.outdir, 'lists'))
    #writeGeneExpression(geneRecords, os.path.join(nspace.outdir, 'lists'))
    writeStatus('Plotting Depths', nspace.verbose)
    f1labs = [ '%s Rep %d' % (nspace.factorlabels[0], i+1) for i in xrange( len(nspace.factor1Dirs))]
    f2labs = [ '%s Rep %d' % (nspace.factorlabels[1], i+1) for i in xrange( len(nspace.factor2Dirs))]
    if nspace.noplot: return
    plotPDist(testedGenes, os.path.join(nspace.outdir, 'figures'))
    plotMVA(testedGenes, aVals, os.path.join(nspace.outdir, 'figures'))

    plotResults( testedGenes, aVals, f1Dict, f2Dict, f1labs+f2labs, 
                 nspace, geneModel, False,
                 os.path.join(nspace.outdir, 'figures'))
    plotResults( testedGenes, aVals, f1Dict, f2Dict, f1labs+f2labs, 
                 nspace, geneModel, True,
                 os.path.join(nspace.outdir, 'figuresLog'))
    
def runExon(geneRecords, geneModel, f1Dict, f2Dict, nspace):
    """
    Run differential exon skipping analysis
    """
    writeStatus("Computing SE statistics", nspace.verbose)
    SErecords, aVals = computeSEStatistics( geneRecords, f1Dict, f2Dict, 
                                            nspace )
    plotPDistSE(SErecords, os.path.join(nspace.outdir, 'figures'))
    plotMVASE(SErecords, aVals, os.path.join(nspace.outdir, 'figures'))
    summaryDictSE = summarySE( SErecords, aVals, nspace.fdrlevel)
    fullTexTableSE(summaryDictSE,os.path.join(nspace.outdir, 'lists')) 
    #writeListsSE( summaryDict, os.path.join(nspace.outdir, 'lists'))
    writeAllSE( SErecords, aVals, os.path.join(nspace.outdir, 'lists'))
    writeStatus("Plotting Depths", nspace.verbose)
    f1labs = [ '%s Rep %d' % (nspace.factorlabels[0], i+1) for i in xrange( len(nspace.factor1Dirs))]
    f2labs = [ '%s Rep %d' % (nspace.factorlabels[1], i+1) for i in xrange( len(nspace.factor2Dirs))]

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
    """
    Run differential event analysis for given event
    """
    if event == 'IR':
        return runIntron

    elif event == 'SE':
        return runExon
        
    else:
        raise ValueError

def main():
    nspace = parseArgs()
    if not makeOutputDir(nspace): 
        sys.exit('Could not create directory: %s' % nspace.outdir)

    writeStatus('Loading models', nspace.verbose)
    geneModel = loadGeneModels( nspace.genemodel, verbose=nspace.verbose )
    writeStatus('Making reduced models', nspace.verbose)
    geneRecords = makeModels( geneModel, verbose=nspace.verbose, 
                              graphDirs=nspace.graphDirs, 
                              exonic=nspace.event=='SE', procs=nspace.procs )
    writeStatus( 'Loading Depths', nspace.verbose )
    f1Dict, f2Dict = loadData( nspace, geneModel, geneRecords )
    genered = [ ]
    for gene in geneRecords:
        if gene.gid in f1Dict and gene.gid in f2Dict and not gene.chrom.isalpha(): genered.append(gene)
    geneRecords = genered

    _dispatch(nspace.event)(geneRecords, geneModel, f1Dict, f2Dict, nspace)

if __name__ == "__main__":
    main()

