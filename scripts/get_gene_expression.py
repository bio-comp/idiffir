#!/bin/python


#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
#    Author: Mike Hamilton, Colorado State University, 2013
#
from iDiffIR.IntronModel import *

from iDiffIR.BamfileIO import *
from iDiffIR.SpliceGrapher.formats.fasta import *
from argparse import ArgumentParser, ArgumentTypeError
import os, sys, numpy
from iDiffIR.SpliceGrapher.shared.utils      import *
from iDiffIR.SpliceGrapher.formats.loader import loadGeneModels



def fileList( raw ):
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
    parser.add_argument('-p', '--procs', dest='procs', action='store', default=1,
                        type=int, help='Number of processing cores to use, [default = 1]')
    parser.add_argument('-o', '--outfile', dest='outfile', action='store', default='expression.txt',
                        type=str, help='output file name')

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
            sys.stderr.write('**bamfile %s not found\n' % f )
            bamfilesOK = False
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

def writeStatus( status ):
    """
    Pretty print status
    """
    n = len(status)+2
    sys.stderr.write( '%s\n' % ( '-' * n ) )
    sys.stderr.write( ' %s \n' % ( status ) )
    sys.stderr.write( '%s\n' % ( '-' * n ) )

def loadData( nspace, geneModel ):
    """
    Load gene depths
    """
    f1Dict = { }
    f2Dict = { }
    def load( factorfiles, fDict ):
        for chrm in sorted(geneModel.getChromosomes()):
            genes     = geneModel.getGeneRecords(chrm)
            genes.sort()

            for i in range(len(factorfiles)):
                fname = os.path.join(factorfiles[i],
                                     '%s.cnt.gz' % chrm)
                if not os.path.exists(fname):
                    sys.stderr.write('Depths file %s not found\n' % fname )
                    continue
                itr = fasta_itr( fname )
                if nspace.verbose:
                    sys.stderr.write("Reading depths from %s\n" % ( fname ) )
                rec = next(itr)
                key = rec.header
                depths = numpy.array( rec.sequence.strip().split(), int )
                counter = 0
                for gene in genes:
                    start = min( gene.start(), gene.end() ) - 1
                    end = max( gene.start(), gene.end() )
                    l = fDict.get(gene.id, [ ] )
                    l.append( depths[start:end] )
                    fDict[gene.id] = l
                    counter += 1
                if nspace.verbose:
                    sys.stderr.write("Read %d gene depths\n" % (counter))
    if nspace.verbose:
        sys.stderr.write( '-' * 30 + '\n' )
        sys.stderr.write("|Reading factor 1 gene depths|\n")
        sys.stderr.write( '-' * 30 + '\n' )

    load( nspace.factor1files, f1Dict)
    if nspace.verbose:
        sys.stderr.write( '-' * 30 + '\n' )
        sys.stderr.write("|Reading factor 2 gene depths|\n")
        sys.stderr.write( '-' * 30 + '\n' )
    load( nspace.factor2files, f2Dict)
    return f1Dict, f2Dict

def main():
    nspace = parseArgs()
    writeStatus('Loading models')
    geneModel = loadGeneModels( nspace.genemodel, verbose=nspace.verbose )
    writeStatus('Making reduced models')
    geneRecords = makeModels( geneModel, None,
                              verbose=nspace.verbose,
                              graphDirs=None,
                              exonic=False,
                              procs=nspace.procs )

    writeStatus('Computing depths')
    ofile = open(nspace.outfile, 'w')
    ofile.write('geneID\tf1Exp_gene\tf2Exp_gene\tf1Exp_IR\tf2Exp_IR\n')
    for gene in geneRecords:
        f1EV, f2EV, f1Juncs, f2Juncs =  getDepthsFromBamfiles( gene,
                                                               nspace.factor1bamfiles,
                                                               nspace.factor2bamfiles
                                                           )
        F1C = numpy.array([ [(f1EV[i][s:(e+1)]).mean() \
                             for s,e in gene.exonsI] for i in range(len( f1EV))]).mean(0)
        F2C = numpy.array([ [(f2EV[i][s:(e+1)]).mean() \
                             for s,e in gene.exonsI] for i in range(len( f2EV))]).mean(0)
        f1depth = numpy.mean(F1C)
        f2depth = numpy.mean(F2C)
        F1I = numpy.array([ [(f1EV[i][s:(e+1)]).mean() \
                             for s,e in gene.introns] for i in range(len( f1EV))]).mean(0)
        ires1 = [ ]
        ires2 = [ ]
        for i,r in enumerate(gene.introns):
            s,e = r
            if gene.retained[i]:
                ires1.append( numpy.array([(f1EV[i][s:(e)]).mean() \
                                          for i in range(len( f1EV))]).mean(0) )
                ires2.append( numpy.array([(f2EV[i][s:(e)]).mean() \
                                          for i in range(len( f2EV))]).mean(0) )


        if len(ires1) > 0:
            ofile.write('%s\t%f\t%f\t%f\t%f\n' % ( gene.gid, f1depth, f2depth, numpy.mean(ires1), numpy.mean(ires2)))
        else:
            ofile.write('%s\t%f\t%f\t-\t-\n' % ( gene.gid, f1depth, f2depth) )
    ofile.close()

if __name__ == "__main__":
    main()

def nspaceS():
    class w: pass
    nspace = w()
    nspace.outdir = '.'
    nspace.krange=[2,6]
    #nspace.genemodel = os.getenv('HUMANGTFPROT')
    nspace.genemodel = os.getenv('TAIRGTF')
    nspace.factor1files = fileList( '/s/waffles/b/tmp/hamiltom/simulation/mutant_idiffir')
    nspace.factor2files = fileList( '/s/waffles/b/tmp/hamiltom/simulation/wildtype_idiffir')
    nspace.verbose = True
    nspace.biasAdj=False
    nspace.numClusts=5
    nspace.coverage=0.90
    nspace.fdrlevel = 0.05
    nspace.shrink_introns=False
