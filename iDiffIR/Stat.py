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
#    Author: Michael Hamilton, Colorado State University, 2013
#    Contact: <hamiltom@cs.colostate.edu>
"""Statistical functions for **iDiffIR**

"""

from SpliceGrapher.SpliceGraph       import *
from scipy.stats import t,gmean, hmean, sem, ttest_1samp
from scipy.stats import norm as sNorm
#from scipy.stats import variance as coeffVar
from scipy.cluster.vq import kmeans,vq
from multi_test import bh,qvalues,bonferroni
from itertools import chain
import numpy, os, sys, multi_test
sys.setrecursionlimit(10000)
from multiprocessing import Pool, freeze_support, Queue, Process
from iDiffIR.BamfileIO import *
from numpy.linalg import norm
EPS=10**-5

class Results(object):
    """
    Container for results
    """
    
    def __init__(self, nspace):
        self.nspace = nspace

    def addLibrarySizes(f1LibSizes, f2LibSizes):
        self.f1LibSizes = f1LibSizes
        self.f2LibSizes = f2LibSizes


def rmsip( geneRecords, f1Dict, f2Dict ):
    cosL = []
    for gene in geneRecords:
        v1 = f1Dict[gene.gid][0]
        v2 = f2Dict[gene.gid][0]
        if numpy.sum(v1) == 0 or numpy.sum(v2) == 0: continue
        cosL.append( numpy.dot( v1, v2) / ( numpy.linalg.norm( v1 ) *  numpy.linalg.norm( v2) ))
    return numpy.sqrt( numpy.sum( [ x**2 for x in cosL ]) / float(len(cosL)) )
        
def removePositionalBias( geneRecords, f1Dict, f2Dict, numClusts, verbose=True ):
    """Removes primer bias by smoothing read depths.

    Depricated
    """
    THREE_PRIME_BIAS = 0
    FIVE_PRIME_BIAS  = 1
    CENTRAL_BIAS     = 2
    TERMINAL_BIAS    = 3
    UNIFORM_BIAS     = 4

    bMap = {0:'3', 1:'5', 2:'C', 3:'T', 4:'U'}
    sys.stderr.write('\tAdjusting read depths for positional bias...\n')
    f1Reps = len(f1Dict[f1Dict.keys()[0]])
    f2Reps = len(f2Dict[f2Dict.keys()[0]])
    f1Codes    = [ list() for _ in xrange( f1Reps ) ]
    f2Codes    = [ list() for _ in xrange( f2Reps ) ]
    f1Cnt      = {}
    f2Cnt      = {}
    f1GeneCodes    = [ list() for _ in xrange( f1Reps ) ]
    f2GeneCodes    = [ list() for _ in xrange( f2Reps ) ]

    f1Weights  = [ list() for _ in xrange( f1Reps ) ]
    f2Weights  = [ list() for _ in xrange( f2Reps ) ]

    soffset = 0
    eoffset = 1

    for i in xrange( f1Reps ):
        allBins = []
        tested  = []
        if verbose:
            sys.stderr.write('\tAdjusting read depths for Factor 1 Rep %d...\n'%(i+1))
        for gene in geneRecords:
            expR = []
            for s,e in gene.exons:
                expR.extend( 
                    f1Dict[gene.gid][i][(s+soffset):(e+eoffset)].tolist())

            tot = float(sum( expR))
            mu = float(numpy.mean(expR))
            # no expression
            if tot == 0: 
                tested.append(False)
            else:
                tested.append(True)
                bins = numpy.array_split(expR, 10)
                allBins.append(numpy.array([ (sum(x)+mu)/tot for x in bins ]))
        allBins = numpy.array(allBins)
        codes, distortion = kmeans(allBins, 
                                   k_or_guess=numClusts, 
                                   thresh=10**(-100))

        f1Codes[i] = codes
        rlabels, dist = vq( allBins, codes)
        for c in xrange( numClusts ):
            f1Cnt[c] = rlabels.tolist().count(c)

        wts = [ [ max(X)/xi for xi in X ] for X in codes ]
        gWts = [ ]
        idx = 0
        for test in tested:
            if test:
                #lWts = numpy.array([ max(allBins[idx]) / xi for xi in allBins[idx]])
                #gWts.append( (wts[rlabels[idx]]+lWts)/2.0 )
                gc = rlabels[idx]
                gWts.append( wts[rlabels[idx]]  )
                f1GeneCodes[i].append(gc)
                idx += 1
            else:
                f1GeneCodes[i].append(None)
                gWts.append(numpy.array([1.0]*10))

        f1Weights[i] = gWts

    for i in xrange( f2Reps ):
        allBins = []
        tested  = []
        if verbose:
            sys.stderr.write('\tAdjusting read depths for Factor 2 Rep %d...\n'%(i+1))
        for gene in geneRecords:
            expR = []
            for s,e in gene.exons:
                expR.extend( 
                    f2Dict[gene.gid][i][(s+soffset):(e+eoffset)].tolist())
            
            tot = float(sum( expR))
            mu = float(numpy.mean(expR))
            # no expression
            if tot == 0: 
                tested.append(False)
            else:
                tested.append(True)
                bins = numpy.array_split(expR, 10)
                allBins.append(numpy.array([ (sum(x)+mu)/tot for x in bins ]))
        allBins = numpy.array(allBins)
        codes, distortion = kmeans(allBins, 
                                   k_or_guess=numClusts, 
                                   thresh=10**(-100))

        f2Codes[i] = codes
        rlabels, dist = vq( allBins, codes)
        for c in xrange( numClusts ):
            f2Cnt[c] = rlabels.tolist().count(c)

        wts = [ [ max(X)/xi for xi in X ] for X in codes ]
        gWts = [ ]
        idx = 0
        for test in tested:
            if test:
                #lWts = numpy.array([ max(allBins[idx]) / xi for xi in allBins[idx]])
                #gWts.append( (wts[rlabels[idx]]+lWts)/2.0 )
                gc = rlabels[idx]
                gWts.append(wts[rlabels[idx]])
                f2GeneCodes[i].append(gc)
                idx += 1
            else:
                f2GeneCodes[i].append(None)
                gWts.append(numpy.array([1.0]*10))

        f2Weights[i] = gWts

    for gidx, gene in enumerate( geneRecords):
        gene.f1BiasWts = [ f1Weights[idx][gidx] for idx in range(f1Reps) ]
        tot = sum( [e-s for s,e in gene.exons])
        fracs = [ (e-s)/float(tot) for s,e in gene.exons]
        gene.f1ExonWts = [ ]
        for i in xrange( f1Reps ):

            wts = [ ]
            offset = 0
            for idx in range(len(gene.exons)):
                s,e = (int(numpy.floor(offset*10)), 
                       int(numpy.ceil( (offset+fracs[idx])*10)))
                wts.append( numpy.mean(gene.f1BiasWts[i][s:e]) )
                offset += fracs[idx]
            wts = [ xi/min(wts) for xi in wts]
            gene.f1ExonWts.append(wts)
            
        gene.f2BiasWts = [ f2Weights[idx][gidx] for idx in range(f2Reps) ]
        gene.f2ExonWts = [ ]

        for i in xrange( f2Reps ):
            wts = [ ]
            offset = 0
            for idx in range(len(gene.exons)):
                s,e = (int(numpy.floor(offset*10)), 
                       int(numpy.ceil( (offset+fracs[idx])*10)))
                wts.append( numpy.mean(gene.f2BiasWts[i][s:e]) )
                offset += fracs[idx]
            wts = [ xi/min(wts) for xi in wts]
            gene.f2ExonWts.append(wts)

    return f1Codes, f2Codes, f1Cnt, f2Cnt

def computeNormFactors( geneRecords, f1Dict, f2Dict, verbose=False ):
    """Compute library norm factors for each replicate and
    separately for both introns and exons

    Depricated
    """
    if verbose:
        sys.stderr.write('\tComputing library size factors...\n')
    f1Reps = len(f1Dict[f1Dict.keys()[0]])
    f2Reps = len(f2Dict[f2Dict.keys()[0]])
    f1ExonExp   = [ list() for _ in xrange( f1Reps ) ]
    f2ExonExp   = [ list() for _ in xrange( f2Reps ) ]
    f1IntronExp = [ list() for _ in xrange( f1Reps ) ]
    f2IntronExp = [ list() for _ in xrange( f2Reps ) ]

    for gene in geneRecords:
        soffset = 0
        eoffset = 1
        for i in xrange( f1Reps ):
            f1ExonExp[i].extend( [ numpy.mean( 
                        f1Dict[gene.gid][i][(s+soffset):(e+eoffset)]) for s,e in gene.exons] )

            f1IntronExp[i].extend( [ numpy.mean( 
                        f1Dict[gene.gid][i][(s+soffset):(e+eoffset)]) for s,e in gene.introns ] )

        for i in xrange( f2Reps ):
            f2ExonExp[i].extend( [ numpy.mean( 
                        f2Dict[gene.gid][i][(s+soffset):(e+eoffset)]) for s,e in gene.exons] )
            f2IntronExp[i].extend( [ numpy.mean( 
                        f2Dict[gene.gid][i][(s+soffset):(e+eoffset)]) for s,e in gene.introns ] )

    k = int(numpy.log2(max([numpy.mean( 
                        [ x for x in b if x > 0 ]) for b in f1ExonExp+f2ExonExp])))
    
    f1ExonExp   = [ sum(x) for x in f1ExonExp ]
    f2ExonExp   = [ sum(x) for x in f2ExonExp ]
    f1IntronExp = [ sum(x) for x in f1IntronExp ]
    f2IntronExp = [ sum(x) for x in f2IntronExp ]

    exonNorm = max(f1ExonExp+f2ExonExp)
    f1ExonExp   = exonNorm / numpy.array( f1ExonExp )
    f2ExonExp   = exonNorm / numpy.array( f2ExonExp )

    intronNorm = max(f1IntronExp+f2IntronExp)

    f1IntronExp = intronNorm / numpy.array( f1IntronExp )
    f2IntronExp = intronNorm / numpy.array( f2IntronExp )

    return f1ExonExp, f2ExonExp, f1IntronExp, f2IntronExp, k

def computeNormFactorsAdj( geneRecords, f1Dict, f2Dict, verbose=False ):
    """Compute library norm factors for each replicate and
    separately for both introns and exons

    Depricated
    """
    if verbose:
        sys.stderr.write('\tComputing library size factors...\n')
    f1Reps = len(f1Dict[f1Dict.keys()[0]])
    f2Reps = len(f2Dict[f2Dict.keys()[0]])
    f1ExonExp   = [ list() for _ in xrange( f1Reps ) ]
    f2ExonExp   = [ list() for _ in xrange( f2Reps ) ]
    f1IntronExp = [ list() for _ in xrange( f1Reps ) ]
    f2IntronExp = [ list() for _ in xrange( f2Reps ) ]

    for gene in geneRecords:
        soffset = 0
        eoffset = 1
        for i in xrange( f1Reps ):
            f1ExonExp[i].extend( [ numpy.mean( 
                        f1Dict[gene.gid][i][(gene.exons[e][0]+soffset):(gene.exons[e][1]+eoffset)]*gene.f1ExonWts[i][e]) for e in range(len(gene.exons)) ] )

            f1IntronExp[i].extend( [ numpy.mean( 
                        f1Dict[gene.gid][i][(gene.introns[idx][0]+soffset):(gene.introns[idx][1]+eoffset)]*(gene.f1ExonWts[i][idx]+gene.f1ExonWts[i][idx+1])/2.0) for idx in range(len(gene.introns)) ] )

        for i in xrange( f2Reps ):
            f2ExonExp[i].extend( [ numpy.mean( 
                        f2Dict[gene.gid][i][(gene.exons[e][0]+soffset):(gene.exons[e][1]+eoffset)]*gene.f2ExonWts[i][e]) for e in range(len(gene.exons)) ] )
            f2IntronExp[i].extend( [ numpy.mean( 
                        f2Dict[gene.gid][i][(gene.introns[idx][0]+soffset):(gene.introns[idx][1]+eoffset)]*(gene.f2ExonWts[i][idx]+gene.f2ExonWts[i][idx+1])/2.0) for idx in range(len(gene.introns)) ] )

    k = int(numpy.log2(max([numpy.mean( 
                        [ x for x in b if x > 0 ]) for b in f1ExonExp+f2ExonExp])))
    
    f1ExonExp   = [ sum(x) for x in f1ExonExp ]
    f2ExonExp   = [ sum(x) for x in f2ExonExp ]
    f1IntronExp = [ sum(x) for x in f1IntronExp ]
    f2IntronExp = [ sum(x) for x in f2IntronExp ]

    exonNorm = max(f1ExonExp+f2ExonExp)
    f1ExonExp   = exonNorm / numpy.array( f1ExonExp )
    f2ExonExp   = exonNorm / numpy.array( f2ExonExp )

    intronNorm = max(f1IntronExp+f2IntronExp)

    f1IntronExp = intronNorm / numpy.array( f1IntronExp )
    f2IntronExp = intronNorm / numpy.array( f2IntronExp )

    return f1ExonExp, f2ExonExp, f1IntronExp, f2IntronExp, k

def procGeneStatsIR_new( tasks, test_status):
    """Wrapper to compute IR statistics for a gene
    """
    for gene, nspace, f1LNorm, f2LNorm in iter(tasks.get, 'STOP'):
        f1Depths, f2Depths, f1Juncs, f2Juncs =  getDepthsFromBamfiles( gene, 
                                                                       nspace.factor1bamfiles, 
                                                                       nspace.factor2bamfiles 
                                                                   )
        aVals = range(nspace.krange[0], nspace.krange[1]+1)
        # nothing to test 

        # .. todo:: change to search for undetected splice juntions
        #           if requested
        if len(gene.introns) == 0:
            test_status.put((gene.gid, False))
            continue
        # filter genes with no expression 
        if numpy.sum( numpy.array(f1Depths).mean(0) ) == 0 or \
                numpy.sum( numpy.array(f2Depths).mean(0) ) == 0:
            test_status.put((gene.gid, False))
            continue

        # create stats arrays for each 
        # IR event
        nIRs = len(gene.introns)
        IRFC = [ ]
        IRfc = [ ]
        IRexp = [ ]
        IRserrs = []
        IRstat = []
        IRTested = [ ]
        f1Depths += 1
        f2Depths += 1
        serr, f1Norm,f2Norm, f1exp, f2exp, fc = computeSE( gene, f1Depths, 
                                                           f2Depths, f1LNorm, 
                                                           f2LNorm)
        gene.serr = serr
        gene.f1Norm = f1Norm
        gene.f2Norm = f2Norm
        w = 5
        for k,r  in enumerate(gene.introns):
            s,e = r
            sr, er = gene.intronsR[k]
            ls,le = gene.exons[k]
            rs,re = gene.exons[k+1]

            # make sure intron and flanking exons exist
            if s+1 >= e or ls+1 >= le or rs+1 >= re:
                IRTested.append(False)
                IRFC.append(None)
                IRfc.append(None)
                IRexp.append(None)
                IRstat.append(None)
                continue

            if False:#biasAdj:
                f1Int = numpy.array([f1ExpV[i][(s+1):(e)] * f1IntronNorm[i] * ((gene.f1ExonWts[i][k] + gene.f1ExonWts[i][k+1]) / 2.0)\
                                         for i in xrange(len(f1IntronNorm))]).mean(0)
                f2Int = numpy.array([f2ExpV[i][(s+1):(e)] * f2LNorm[i] * ((gene.f2ExonWts[i][k] + gene.f2ExonWts[i][k+1]) / 2.0)\
                                         for i in xrange(len(f2LNorm))]).mean(0)

            else:
                f1Int = numpy.array([f1Depths[i][(s+1):(e)] * f1LNorm[i]\
                                         for i in xrange(len(f1LNorm))]).mean(0)
                f2Int = numpy.array([f2Depths[i][(s+1):(e)] * f2LNorm[i]\
                                         for i in xrange(len(f2LNorm))]).mean(0)

                f1IntR = numpy.array([(f1Depths[i][(s+1):(e)]-1) * f1LNorm[i]\
                                         for i in xrange(len(f1LNorm))]).mean(0)
                f2IntR = numpy.array([(f2Depths[i][(s+1):(e)]-1) * f2LNorm[i]\
                                         for i in xrange(len(f2LNorm))]).mean(0)

                f1Irv =  numpy.array([f1Depths[i][(s-w):(e+w)] * f1LNorm[i]\
                                         for i in xrange(len(f1LNorm))]).mean(0)
                f2Irv =  numpy.array([f2Depths[i][(s-w):(e+w)] * f2LNorm[i]\
                                         for i in xrange(len(f2LNorm))]).mean(0)
                f1Irv /= norm(f1Irv)
                f2Irv /= norm(f2Irv)


            if False:#biasAdj:
                f1ExonL = numpy.array([f1ExpV[i][ls:le]*f1LNorm[i]*gene.f1ExonWts[i][k] \
                                           for i in xrange(len(f1LNorm))]).mean(0)
                f1ExonR = numpy.array([f1ExpV[i][rs:re]*f1LNorm[i]*gene.f1ExonWts[i][k+1] \
                                           for i in xrange(len(f1LNorm))]).mean(0)
                f2ExonL = numpy.array([f2ExpV[i][ls:le]*f2LNorm[i]*gene.f2ExonWts[i][k] \
                                           for i in xrange(len(f2LNorm))]).mean(0)
                f2ExonR = numpy.array([f2ExpV[i][rs:re]*f2LNorm[i]*gene.f2ExonWts[i][k+1] \
                                           for i in xrange(len(f2LNorm))]).mean(0)
            else:
                f1ExonL = numpy.array([(f1Depths[i][ls:le]-1)*f1LNorm[i] \
                                           for i in xrange(len(f1LNorm))]).mean(0)
                f1ExonR = numpy.array([(f1Depths[i][rs:re]-1)*f1LNorm[i] \
                                           for i in xrange(len(f1LNorm))]).mean(0)
                f2ExonL = numpy.array([(f2Depths[i][ls:le]-1)*f2LNorm[i] \
                                           for i in xrange(len(f2LNorm))]).mean(0)
                f2ExonR = numpy.array([(f2Depths[i][rs:re]-1)*f2LNorm[i] \
                                           for i in xrange(len(f2LNorm))]).mean(0)

            #intserrs.append(serr)
            f1Intm = numpy.mean(f1Int)
            f2Intm = numpy.mean(f2Int)
            f1m = numpy.mean( f1Irv[(w):(-w)] )
            f2m = numpy.mean( f2Irv[(w):(-w)] )
            f1IntmR = numpy.mean(f1IntR)
            f2IntmR = numpy.mean(f2IntR)

            f1LSS = numpy.array(list(f1ExonL[-10:]) + list(f1IntR[:11])) + 1
            f1RSS = numpy.array(list(f1IntR[-10:]) + list(f1ExonR[:11])) + 1

            f2LSS = numpy.array(list(f2ExonL[-10:]) + list(f2IntR[:11])) + 1
            f2RSS = numpy.array(list(f2IntR[-10:]) + list(f2ExonR[:11])) + 1

            cosL = numpy.dot( f1LSS, f2LSS) / (EPS+numpy.linalg.norm(f1LSS) * numpy.linalg.norm(f2LSS))
            cosR = numpy.dot( f1RSS, f2RSS) / (EPS+numpy.linalg.norm(f1RSS) * numpy.linalg.norm(f2RSS))

            f1exon = 0.5 * ( numpy.mean(f1ExonL) + numpy.mean(f1ExonR) )
            f2exon = 0.5 * ( numpy.mean(f2ExonL) + numpy.mean(f2ExonR) )

            tested = True
            # if ambiguity exists from antisense gene, don't test!
            for os, oe in gene.overlap:
                if min(oe, er) - max(sr, os) > 0.25*(er-sr+1):
                    tested = False
                    break
            f1Coverage = numpy.sum(f1IntR >= 1 )/float(len(f1IntR))
            f2Coverage = numpy.sum(f2IntR >= 1 )/float(len(f2IntR))
            #if numpy.sum(f1LSS > 0 ) / float(len(f1LSS)) < 1 or numpy.sum(f1RSS > 0 ) / float(len(f1RSS)) < 1 or \
            #   numpy.sum(f2LSS > 0 ) / float(len(f2LSS)) < 1 or numpy.sum(f2RSS > 0 ) / float(len(f2RSS)) < 1:
            #    tested = False
            #if numpy.sum(f1LSS) == 0 or numpy.sum(f1RSS) == 0 or numpy.sum(f2LSS) == 0 or numpy.sum(f2RSS) == 0:
            #    tested = False

            if f1Coverage < nspace.coverage and f2Coverage < nspace.coverage: 
                tested = False

            # check if gene is DE and if sufficient read 
            # depth exists in the lower-expressed condition
            if tested:
                if gene.f1Norm < gene.f2Norm and f1IntmR > f2IntmR * f2Norm:
                    if  max(f1Norm, f2Norm) > nspace.dexpThresh:
                        tested = False
                    if f2Coverage < f1Coverage:
                        tested = False


                elif gene.f2Norm < gene.f1Norm and f2IntmR > f1IntmR * f1Norm:
                    if  max(f1Norm, f2Norm) > nspace.dexpThresh:
                        tested = False
                    if f1Coverage < f2Coverage:
                        tested = False

            if not tested:
                IRTested.append(True)
                IRFC.append(0)
                IRfc.append(0)
                IRexp.append(1)
                #IRstat.append(0)


            else:
                IRTested.append(True)
                # compute numerators for each value of a
                numer = [numpy.log2( 0 + f1m ) \
                             for a in aVals ]

                # compute denominators for each value of a
                denom = [numpy.log2( 0 + f2m) \
                             for a in aVals ]

                IRfc.append( numpy.log2(f1m) -  \
                             numpy.log2(f2m))

                IRexp.append(0.5 * (numpy.log2(f1Intm * f1Norm+EPS)
                                     + numpy.log2(f2Intm *f2Norm+EPS)))


                IRFC.append( [(numer[i] - denom[i]) / (cosL + cosR) \
                                   for i in xrange(len(numer))])

                IRstat.append( [ (numer[i] - denom[i] ) / gene.serr \
                                  for i in xrange(len(numer))])
        if IRTested.count(True) > 0:
            test_status.put((gene.gid, True, IRexp, IRfc, IRTested, IRstat))
        else:
            test_status.put((gene.gid, False))

def procGeneStatsIR( tasks, test_status):
    """Wrapper to compute IR statistics for a gene
    """
    for gene, nspace, f1LNorm, f2LNorm in iter(tasks.get, 'STOP'):
        f1Depths, f2Depths, f1Juncs, f2Juncs =  getDepthsFromBamfiles( gene, 
                                                                       nspace.factor1bamfiles, 
                                                                       nspace.factor2bamfiles 
                                                                   )
        aVals = range(nspace.krange[0], nspace.krange[1]+1)
        # nothing to test 

        # .. todo:: change to search for undetected splice juntions
        #           if requested
        if len(gene.introns) == 0:
            test_status.put((gene.gid, False))
            continue
        # filter genes with no expression or too much differential gene expression
        if numpy.sum( numpy.array(f1Depths).mean(0) ) == 0 or \
                numpy.sum( numpy.array(f2Depths).mean(0) ) == 0:
            test_status.put((gene.gid, False))
            continue

        # create stats arrays for each 
        # IR event
        nIRs = len(gene.introns)
        IRFC = [ ]
        IRfc = [ ]
        IRexp = [ ]
        IRserrs = []
        IRstat = []
        IRTested = [ ]
        IRrat1 = [ ]
        IRrat2 = [ ]
        serr, f1Norm,f2Norm, f1exp, f2exp, f1expr, f2expr, fc = computeSE( gene, f1Depths, 
                                                                           f2Depths, f1LNorm, 
                                                                           f2LNorm)
        if max(f1Norm, f2Norm) > nspace.dexpThresh:
            test_status.put((gene.gid, False))
            continue
            
        f1Depths += 1
        f2Depths += 1
        
        gene.serr = serr
        gene.f1Norm = f1Norm
        gene.f2Norm = f2Norm
        gene.f1exp = f1exp
        gene.f2exp = f2exp

        for k,r  in enumerate(gene.introns):
            s,e = r
            sr, er = gene.intronsR[k]
            ls,le = gene.exons[k]
            rs,re = gene.exons[k+1]

            # make sure intron and flanking exons exist
            if s+1 >= e or ls+1 >= le or rs+1 >= re:
                IRTested.append(False)
                IRFC.append(None)
                IRfc.append(None)
                IRexp.append(None)
                IRstat.append(None)
                IRrat1.append(None)
                IRrat2.append(None)

                continue

            if False:#biasAdj:
                f1Int = numpy.array([f1ExpV[i][(s+1):(e)] * f1IntronNorm[i] * ((gene.f1ExonWts[i][k] + gene.f1ExonWts[i][k+1]) / 2.0)\
                                         for i in xrange(len(f1IntronNorm))]).mean(0)
                f2Int = numpy.array([f2ExpV[i][(s+1):(e)] * f2LNorm[i] * ((gene.f2ExonWts[i][k] + gene.f2ExonWts[i][k+1]) / 2.0)\
                                         for i in xrange(len(f2LNorm))]).mean(0)

            else:
                f1Int = numpy.array([f1Depths[i][(s+1):(e)] * f1LNorm[i]\
                                     for i in xrange(len(f1LNorm))]).mean(0)
                f2Int = numpy.array([f2Depths[i][(s+1):(e)] * f2LNorm[i]\
                                     for i in xrange(len(f2LNorm))]).mean(0)

                f1IntR = numpy.array([(f1Depths[i][(s+1):(e)]-1) * f1LNorm[i]\
                                      for i in xrange(len(f1LNorm))]).mean(0)
                f2IntR = numpy.array([(f2Depths[i][(s+1):(e)]-1) * f2LNorm[i]\
                                      for i in xrange(len(f2LNorm))]).mean(0)
                f1Intexp = numpy.array([(f1Depths[i][(s+1):(e)]-1) \
                                      for i in xrange(len(f1LNorm))]).mean(0)
                f2Intexp = numpy.array([(f2Depths[i][(s+1):(e)]-1) \
                                      for i in xrange(len(f2LNorm))]).mean(0)

                IRrat1.append(f1Intexp/f1expr)
                IRrat2.append(f2Intexp/f2expr)

            if False:#biasAdj:
                f1ExonL = numpy.array([f1ExpV[i][ls:le]*f1LNorm[i]*gene.f1ExonWts[i][k] \
                                           for i in xrange(len(f1LNorm))]).mean(0)
                f1ExonR = numpy.array([f1ExpV[i][rs:re]*f1LNorm[i]*gene.f1ExonWts[i][k+1] \
                                           for i in xrange(len(f1LNorm))]).mean(0)
                f2ExonL = numpy.array([f2ExpV[i][ls:le]*f2LNorm[i]*gene.f2ExonWts[i][k] \
                                           for i in xrange(len(f2LNorm))]).mean(0)
                f2ExonR = numpy.array([f2ExpV[i][rs:re]*f2LNorm[i]*gene.f2ExonWts[i][k+1] \
                                           for i in xrange(len(f2LNorm))]).mean(0)
            else:
                f1ExonL = numpy.array([(f1Depths[i][ls:le]-1)*f1LNorm[i] \
                                           for i in xrange(len(f1LNorm))]).mean(0)
                f1ExonR = numpy.array([(f1Depths[i][rs:re]-1)*f1LNorm[i] \
                                           for i in xrange(len(f1LNorm))]).mean(0)
                f2ExonL = numpy.array([(f2Depths[i][ls:le]-1)*f2LNorm[i] \
                                           for i in xrange(len(f2LNorm))]).mean(0)
                f2ExonR = numpy.array([(f2Depths[i][rs:re]-1)*f2LNorm[i] \
                                           for i in xrange(len(f2LNorm))]).mean(0)

            #intserrs.append(serr)
            f1Intm = numpy.mean(f1Int)
            f2Intm = numpy.mean(f2Int)

            f1IntmR = numpy.mean(f1IntR)
            f2IntmR = numpy.mean(f2IntR)

            f1LSS = numpy.array(list(f1ExonL[-10:]) + list(f1IntR[:11])) + 1
            f1RSS = numpy.array(list(f1IntR[-10:]) + list(f1ExonR[:11])) + 1

            f2LSS = numpy.array(list(f2ExonL[-10:]) + list(f2IntR[:11])) + 1
            f2RSS = numpy.array(list(f2IntR[-10:]) + list(f2ExonR[:11])) + 1

            cosL = numpy.dot( f1LSS, f2LSS) / (EPS+numpy.linalg.norm(f1LSS) * numpy.linalg.norm(f2LSS))
            cosR = numpy.dot( f1RSS, f2RSS) / (EPS+numpy.linalg.norm(f1RSS) * numpy.linalg.norm(f2RSS))

            f1exon = 0.5 * ( numpy.mean(f1ExonL) + numpy.mean(f1ExonR) )
            f2exon = 0.5 * ( numpy.mean(f2ExonL) + numpy.mean(f2ExonR) )

            tested = True
            # if ambiguity exists from antisense gene, don't test!
            for os, oe in gene.overlap:
                if min(oe, er) - max(sr, os) > 0.25*(er-sr+1):
                    tested = False
                    break
            f1Coverage = numpy.sum(f1IntR >= 1 )/float(len(f1IntR))
            f2Coverage = numpy.sum(f2IntR >= 1 )/float(len(f2IntR))
            #if numpy.sum(f1LSS > 0 ) / float(len(f1LSS)) < 1 or numpy.sum(f1RSS > 0 ) / float(len(f1RSS)) < 1 or \
            #   numpy.sum(f2LSS > 0 ) / float(len(f2LSS)) < 1 or numpy.sum(f2RSS > 0 ) / float(len(f2RSS)) < 1:
            #    tested = False
            #if numpy.sum(f1LSS) == 0 or numpy.sum(f1RSS) == 0 or numpy.sum(f2LSS) == 0 or numpy.sum(f2RSS) == 0:
            #    tested = False

            if max(f1Coverage, f2Coverage) < nspace.coverage:
                tested = False
            # check if gene is DE and if sufficient read 
            # depth exists in the lower-expressed condition
            if tested:
                if gene.f1Norm < gene.f2Norm:
#                    if f2IntmR == 0 and (f2exon+EPS) * min(1.0, f1IntmR / (f1exon+EPS) ) < 1.0/len(f1Int):
#                        tested = False
                    if f2Coverage < f1Coverage:# and (f2exon+EPS) * min(1.0, f1IntmR / (f1exon+EPS) ) < 1.0/len(f1Int):
                        tested = False

                elif gene.f2Norm < gene.f1Norm:
#                    if f1IntmR == 0 and (f1exon+EPS)* min(1.0, f2IntmR / (f2exon+EPS) ) < 1.0/len(f2Int):
#                        tested = False
                    if f1Coverage < f2Coverage: # and (f1exon+EPS)* min(1.0, f2IntmR / (f2exon+EPS) ) < 1.0/len(f2Int):
                        tested = False

            if not tested:
                IRTested.append(False)
                IRFC.append(None)
                IRfc.append(None)
                IRexp.append(None)
                IRstat.append(None)


            else:
                IRTested.append(True)
                # compute numerators for each value of a
                numer = [numpy.log2( 2**a + f1Intm * f1Norm  ) \
                             for a in aVals ]

                # compute denominators for each value of a
                denom = [numpy.log2( 2**a + f2Intm * f2Norm ) \
                             for a in aVals ]

                IRfc.append( (numpy.log2(f1Intm * f1Norm+EPS)) -  \
                                  (numpy.log2(f2Intm * f2Norm+EPS)))

                IRexp.append(0.5 * (numpy.log2(f1Intm * f1Norm+EPS)
                                     + numpy.log2(f2Intm *f2Norm+EPS)))


                IRFC.append( [(numer[i] - denom[i]) / (cosL + cosR) \
                                   for i in xrange(len(numer))])

                IRstat.append( [ ((numer[i] - denom[i] ) / (cosL + cosR))/gene.serr \
                                  for i in xrange(len(numer))])


               
        if IRTested.count(True) > 0:
            test_status.put((gene.gid, True, IRexp, IRfc, IRTested, IRstat, IRrat1, IRrat2))
        else:
            test_status.put((gene.gid, False))


def overlap(start, stop, start2, stop2):
    if start==start2 or stop==stop2:
        return True
    if start2 >= start and start2 <= stop:
        return True
    if start >= start2 and start <= stop2:
        return True
    if start <= start2 and stop2 <= stop:
        return True
    if start2 <= start and stop <= stop2:
        return True
    return False

def get_weight( start, stop, gene, wts):
    retn = [ ]
    for i, exon in enumerate(gene.exons):
        start2, stop2 = exon
        if overlap( start, stop, start2, stop2 ):
            retn.append(wts[i])
    return numpy.mean(retn)

def testIR(geneRecords, aVals, nspace):
    """Hypothesis testing using Z-scores of 
    modified $\log FC$ statisitc
    """
    #compute Z-score distribution parameters
    for aidx in xrange(len(aVals)):
        X = list(chain.from_iterable([ numpy.array([geneRecords[i].IRstat[t][aidx]\
                                                    for t in xrange(len(geneRecords[i].IRstat)) \
                                                    if geneRecords[i].IRTested[t] and geneRecords[i].IRGTested]) \
                                       for i in xrange(len(geneRecords)) if geneRecords[i].IRGTested]) )
        mu = numpy.mean(X)
        #mu = numpy.median(X)

        sigma = numpy.std(X)

        pvals = [ ]
        N = len(geneRecords)
        # Assign z-scores and pvalues 
        for gene in geneRecords:
            if not gene.IRGTested: continue
            for i in xrange( len( gene.IRTested)):
                if gene.IRTested[i]:
                    z = zscore(gene.IRstat[i][aidx], mu, sigma)
                    gene.IRZ[i].append( z )
                    gene.IRPvals[i].append(min(1.0,2*sNorm.cdf( -abs(z), loc=0, scale=1.0)))
                    pvals.append(gene.IRPvals[i][-1])
                else:
                    gene.IRZ[i].append(0)
                    gene.IRPvals[i].append(2)
        assert len(pvals) == len(X)
        if nspace.multTest == 'BH':
            qvals = bh(pvals)
        elif nspace.multTest == 'BF':
            qvals = bonferroni(pvals)
        elif nspace.multTest == 'QV':
            qvals = qvalues(pvals)[0]

        i = 0
        for gene in geneRecords:
            if not gene.IRGTested: continue
            ci = 0
            tested = gene.IRTested.count(True)
            for idx in xrange(len(gene.introns)):
                if gene.IRTested[idx]:
                    gene.IRQvals[idx].append(qvals[i+ci])
                    ci += 1
                else:
                    gene.IRQvals[idx].append(2)
            i = i + tested


def computeIRStatistics(geneRecords, nspace, validChroms, f1LNorm, f2LNorm):
    """Compute SE statistics

    Parameters
    ----------
    geneRecords : list
                  List of iDiffIR.IntronModel.IntronModel objects to 
                  process.
    nspace : argparse.ArgumentParser
            Command line arguments for **idiffir.py**
    validChroms : set
                  Set of valid chromosomes to test. All depth files
                  contain them.
    f1Norm : list
             List of effective library size norm scaling values for 
             factor 1
    f2Norm : list
             List of effective library size norm scaling values for 
             factor 2

    """
    indicator = ProgressIndicator(10000, verbose=nspace.verbose)
    exonExp = [ ] # list of exonic expression across all genes
    intronExp = [ ] # list of exonic expression across all genes
    if nspace.verbose:
        sys.stderr.write('Computing differential splicing scores...\n')

    aVals = range(nspace.krange[0], nspace.krange[1]+1)
    geneStatus = { }

    # parallel call
    if 1:
        task_queue = Queue()
        status_queue = Queue()
        nTasks = 0
        for gene in geneRecords:
            if gene.chrom in validChroms: 
                task_queue.put( (gene, nspace, f1LNorm, f2LNorm) )
                nTasks += 1

        for _ in xrange(nspace.procs):
            task_queue.put('STOP')

        for _ in xrange(nspace.procs):
            Process(target=procGeneStatsIR,
                    args=(task_queue, status_queue)).start() 

        for _ in xrange(nTasks):
            result = status_queue.get()
            if len(result) > 2:
                geneStatus[result[0]] = result[2:]
            indicator.update()

        for gene in geneRecords:
            if gene.gid in geneStatus:
                gene.IRexp, gene.IRfc, gene.IRTested,\
                    gene.IRstat, gene.IRRats, gene.IRexp_1, \
                    gene.IRexp_2 = geneStatus[gene.gid]
                gene.IRGTested = True
                gene.IRQvals = [ list() for _ in gene.introns]
                gene.IRZ     = [ list() for _ in gene.introns]
                gene.IRPvals = [ list() for _ in gene.introns]
            else:
                gene.IRGTested = False

    indicator.finish()
    testIR(geneRecords, aVals, nspace)

                 
def procGeneStatsSE( tasks, test_status):
    """Compute exon skipping statistics for gene

    Compute exon skipping statistics for all SEs in gene
    by calculating the adjusted :math:`\log`-fold change 
    statistic.  This function can be run in parallel or
    serial.  

    See Also
    --------
    computeSEStatistics : Wrapper function
    getReadDepths : Get read depths from bamfiles
    
    Parameters
    ----------
    gene : iDiffIR.IntronModel.IntronModel 
           Gene to compute SE statistics
    nspace : argparse.ArgumentParser
             Command line arguments from **idiffir.py** 
    f1LNorm, f2LNorm : list
                       Contains library size normalization
                       scalars of type float
    Returns
    -------
    status : bool
             True if gene has tested SE events
    gene : Reference to gene object

    """
    for gene, nspace, f1LNorm, f2LNorm in iter(tasks.get, 'STOP'):
        f1Depths, f2Depths, f1Juncs, f2Juncs =  getDepthsFromBamfiles( gene, 
                                                                       nspace.factor1bamfiles, 
                                                                       nspace.factor2bamfiles 
                                                                   )
        SEnodes = gene.flavorDict['SE']
        aVals = range(nspace.krange[0], nspace.krange[1]+1)
        # nothing to test 
        # .. todo:: change to search for undetected splice juntions
        #           if requested
        if not SEnodes: 
            test_status.put((gene.gid, False))
            continue
        # filter genes with no expression 
        if numpy.sum( numpy.array(f1Depths).mean(0) ) == 0 or \
                numpy.sum( numpy.array(f2Depths).mean(0) ) == 0:
            test_status.put((gene.gid, False))
            continue

        # create stats arrays for each 
        # SE event
        nSEs = len(SEnodes)
        SEFC = [ ]
        SEfc = [ ]
        SEexp = [ ]
        SEserrs = []
        SEstat = []
        SETested = [ ]
        f1Depths += 1
        f2Depths += 1
        serr, f1Norm,f2Norm, f1exp, f2exp, f1expr, f2expr, fc = computeSE( gene, f1Depths, 
                                                                           f2Depths, f1LNorm, 
                                                                           f2LNorm)
        gene.serr = serr
        gene.f1Norm = f1Norm
        gene.f2Norm = f2Norm

        # iterate through each event
        for k,event  in enumerate(SEnodes):
            l,se,r = event
            ls, le = l
            s, e = se
            e += 1
            rs, re = r
            if ls+1 >= le or s+1 >= e or rs +1 >= re:
                SETested.append(False)
                SEFC.append(None)
                SEfc.append(None)
                SEexp.append(None)
                SEstat.append(None)
                continue
            if False:#biasAdj:
                f1SE = numpy.array([f1ExpV[i][(s):(e)] * f1LNorm[i] * get_weight(s,e,gene,f1ExonWts[i])\
                                        for i in xrange(len(f1LNorm))]).mean(0)
                f2SE = numpy.array([f2ExpV[i][(s):(e)] * f2LNorm[i] * get_weight(s,e,gene,f2ExonWts[i])\
                                        for i in xrange(len(f2LNorm))]).mean(0)
            else:
                f1SE = numpy.array([f1Depths[i][s:e] * f1LNorm[i]\
                                        for i in xrange(len(f1LNorm))]).mean(0)
                f2SE = numpy.array([f2Depths[i][s:e] * f2LNorm[i]\
                                        for i in xrange(len(f2LNorm))]).mean(0)

                f1SER = numpy.array([(f1Depths[i][s:e]-1) * f1LNorm[i]\
                                        for i in xrange(len(f1LNorm))]).mean(0)
                f2SER = numpy.array([(f2Depths[i][s:e]-1) * f2LNorm[i]\
                                        for i in xrange(len(f2LNorm))]).mean(0)
            if False:#biasAdj:
                f1ExonL = numpy.array([f1ExpV[i][ls:le]*f1LNorm[i]*get_weight( ls, le, f1ExonWts[i]) \
                                           for i in xrange(len(f1LNorm))]).mean(0)
                f1ExonR = numpy.array([f1ExpV[i][rs:re]*f1LNorm[i]*get_weight( rs, re, f1ExonWts[i]) \
                                           for i in xrange(len(f1LNorm))]).mean(0)
                f2ExonL = numpy.array([f2ExpV[i][ls:le]*f2LNorm[i]*get_weight( ls, le, f2ExonWts[i]) \
                                           for i in xrange(len(f2LNorm))]).mean(0)
                f2ExonR = numpy.array([f2ExpV[i][rs:re]*f2LNorm[i]*get_weight( rs, re, f2ExonWts[i]) \
                                           for i in xrange(len(f2LNorm))]).mean(0)
            else:
                f1ExonL = numpy.array([f1Depths[i][ls:le]*f1LNorm[i] \
                                           for i in xrange(len(f1LNorm))]).mean(0)
                f1ExonR = numpy.array([f1Depths[i][rs:re]*f1LNorm[i] \
                                           for i in xrange(len(f1LNorm))]).mean(0)
                f2ExonL = numpy.array([f2Depths[i][ls:le]*f2LNorm[i] \
                                           for i in xrange(len(f2LNorm))]).mean(0)
                f2ExonR = numpy.array([f2Depths[i][rs:re]*f2LNorm[i] \
                                           for i in xrange(len(f2LNorm))]).mean(0)

            tested = True
            if len(f1SE) == 0 :
                tested = False
                f1SEm = EPS
                f2SEm = EPS
            else:
                f1SEm = numpy.mean(f1SE)
                f2SEm = numpy.mean(f2SE)

            if tested:
                f1LSS = numpy.array(list(f1ExonL[-10:]) + list(f1SE[:11]))
                f1RSS = numpy.array(list(f1SE[-10:]) + list(f1ExonR[:11]))

                f2LSS = numpy.array(list(f2ExonL[-10:]) + list(f2SE[:11]))
                f2RSS = numpy.array(list(f2SE[-10:]) + list(f2ExonR[:11]))

                cosL = numpy.dot( f1LSS, f2LSS) / (EPS+numpy.linalg.norm(f1LSS) * numpy.linalg.norm(f2LSS))
                cosR = numpy.dot( f1RSS, f2RSS) / (EPS+numpy.linalg.norm(f1RSS) * numpy.linalg.norm(f2RSS))

                f1exon = 0.5 * ( numpy.mean(f1ExonL) + numpy.mean(f1ExonR) )
                f2exon = 0.5 * ( numpy.mean(f2ExonL) + numpy.mean(f2ExonR) )

                f1Coverage = numpy.sum(f1SER > 0 )/float(len(f1SER))
                f2Coverage = numpy.sum(f2SER > 0 )/float(len(f2SER))

                # check if event has read coverage
                if f1Coverage < nspace.cverage and f2Coverage < nspace.coverage:

                    tested = False

            # check if gene is DE and if sufficient read 
            # depth exists in the lower-expressed condition
            if tested:
                if gene.f1Norm < gene.f2Norm:
                    if  max(f1Norm, f2Norm) > nspace.dexpThresh or \
                        numpy.sum(f2LSS > 1 ) / float(len(f2LSS)) < 1 or \
                        numpy.sum(f2RSS > 1 ) / float(len(f2RSS)) < 1:
                        tested = False
                    if f1SEm > f2SEm * f2Norm and f2SEm == 0 and (f2exon+EPS) * min(1.0, f1SEm / (f1exon+EPS) ) < 1.0/len(f1SE):

                        tested = False

                elif gene.f2Norm < gene.f1Norm:
                    if  max(f1Norm, f2Norm) > nspace.dexpThresh or \
                        numpy.sum(f1LSS > 1 ) / float(len(f1LSS)) < 1 or \
                        numpy.sum(f1RSS > 1 ) / float(len(f1RSS)) < 1:
                        tested = False

                    if f2SEm > f1SEm * f1Norm and f1SEm == 0 and (f1exon+EPS)* min(1.0, f2SEm / (f2exon+EPS) ) < 1.0/len(f2SE):            
                        tested = False

            if not tested:
                SETested.append(False)
                SEFC.append(None)
                SEfc.append(None)
                SEexp.append(None)
                SEstat.append(None)
            else:
                SETested.append(True)
                numer = [numpy.log2( 2**a + f1SEm * f1Norm  ) \
                             for a in aVals ]
                denom = [numpy.log2( 2**a + f2SEm * f2Norm ) \
                             for a in aVals ]


                SEfc.append( (numpy.log2(f1SEm * f1Norm+EPS)) -  \
                                  (numpy.log2(f2SEm * f2Norm+EPS)))
                SEexp.append(0.5 * (numpy.log2(f1SEm * f1Norm+EPS)\
                                         + numpy.log2(f2SEm *\
                                                          f2Norm+EPS)))

                SEFC.append( [(numer[i] - denom[i]) \
                                   for i in xrange(len(numer))])

                SEstat.append( [ ( numer[i] - denom[i] ) / gene.serr \
                                      for i in xrange(len(numer))])

        if SETested.count(True) > 0:
            test_status.put((gene.gid, True, SEexp, SEfc, SETested, SEstat))
        else:
            test_status.put((gene.gid, False))

def testSE(geneRecords, aVals, nspace):
    """Hypothesis testing using Z-scores of 
    modified $\log FC$ statisitc
    """
    #compute Z-score distribution parameters
    for aidx in xrange(len(aVals)):
        X = list(chain.from_iterable([ numpy.array([geneRecords[i].SEstat[t][aidx]\
                                                    for t in xrange(len(geneRecords[i].SEstat)) \
                                                    if geneRecords[i].SETested[t] and geneRecords[i].SEGTested]) \
                                       for i in xrange(len(geneRecords)) if geneRecords[i].SEGTested]) )
        mu = numpy.mean(X)
        sigma = numpy.std(X)
        pvals = [ ]
        N = len(geneRecords)
        # Assign z-scores and pvalues 
        for gene in geneRecords:
            if not gene.SEGTested: continue
            for i in xrange( len( gene.SETested)):
                if gene.SETested[i]:
                    z = zscore(gene.SEstat[i][aidx], mu, sigma)
                    gene.SEZ[i].append( z )
                    gene.SEPvals[i].append(min(1.0,2*sNorm.cdf( -abs(z), loc=0, scale=1.0)))
                    pvals.append(gene.SEPvals[i][-1])
                else:
                    gene.SEZ[i].append(0)
                    gene.SEPvals[i].append(2)
        assert len(pvals) == len(X)
        if nspace.multTest == 'BH':
            qvals = bh(pvals)
        elif nspace.multTest == 'BF':
            qvals = bonferroni(pvals)
        elif nspace.multTest == 'QV':
            qvals = qvalues(pvals)[0]

        i = 0
        for gene in geneRecords:
            if not gene.SEGTested: continue
            ci = 0
            tested = gene.SETested.count(True)
            for idx in xrange(len(gene.flavorDict['SE'])):
                if gene.SETested[idx]:
                    gene.SEQvals[idx].append(qvals[i+ci])
                    ci += 1
                else:
                    gene.SEQvals[idx].append(2)
            i = i + tested

def computeSEStatistics(geneRecords, nspace, validChroms, f1LNorm, f2LNorm):
    """Compute SE statistics

    Parameters
    ----------
    geneRecords : list
                  List of iDiffIR.IntronModel.IntronModel objects to 
                  process.
    nspace : argparse.ArgumentParser
            Command line arguments for **idiffir.py**
    validChroms : set
                  Set of valid chromosomes to test. All depth files
                  contain them.
    f1Norm : list
             List of effective library size norm scaling values for 
             factor 1
    f2Norm : list
             List of effective library size norm scaling values for 
             factor 2

    Returns
    -------
    testedGenes : list
                  List of genes containing tested IR events
    aVals : list
            List of pseudo-count parameters, :math:`a`, used

    """
    indicator = ProgressIndicator(10000, verbose=nspace.verbose)
    exonExp = [ ] # list of exonic expression across all genes
    intronExp = [ ] # list of exonic expression across all genes
    if nspace.verbose:
        sys.stderr.write('Computing differential splicing scores...\n')
    status_queue = Queue()
    aVals = range(nspace.krange[0], nspace.krange[1]+1)
    geneStatus = { }
    # parallel call
    if True:
        #freeze_support()
        task_queue = Queue()
        nTasks = 0
        for gene in geneRecords:
            if gene.chrom in validChroms:
                task_queue.put( (gene, nspace, f1LNorm, f2LNorm) )
                nTasks += 1
            else:
                gene.SEGTested = False
        for _ in xrange(nspace.procs):
            Process(target=procGeneStatsSE, 
                    args=(task_queue, status_queue)).start()

        for _ in xrange(nspace.procs):
            task_queue.put('STOP')

        for _ in xrange(nTasks):
            result = status_queue.get()
            if len(result) > 2:
                geneStatus[result[0]] = result[2:]
            indicator.update()
        for gene in geneRecords:
            if gene.gid in geneStatus:
                gene.SEexp, gene.SEfc, gene.SETested, gene.SEstat = geneStatus[gene.gid]
                gene.SEGTested = True
                gene.SEQvals = [ list() for _ in gene.flavorDict['SE']]
                gene.SEZ     = [ list() for _ in gene.flavorDict['SE']]
                gene.SEPvals = [ list() for _ in gene.flavorDict['SE']]

            else:
                gene.SEGTested = False
                
    else:
        # serial call
        for gene in geneRecords:
            # skip genes not in bamfiles
            if gene.chrom not in validChroms: 
                gene.SEGTested = False
                continue

            res = procGeneStatsSE((gene, nspace, f1LNorm, f2LNorm), status_queue)
            indicator.update()
    indicator.finish()
    testSE(geneRecords, aVals, nspace)


def geoStd( x ):
    """
    Compute geometric standard deviation
    """
    x = numpy.array(x)
    gMean = gmean(x)
    return 2**numpy.sqrt(numpy.sum(numpy.log2( x / gMean )**2 ) / len( x ) )

def avg_two_variances(n1, n2, s1, s2, mu1, mu2):
    """
    Average two standard deviations
    """
    numer = n1**2*s1**2 + n2**2*s2**2 - n2*s1**2 - n2*s2**2 - \
            n1*s2**2 - n1*s1**2 + n1*n2*s2**2 + n1*n2*s1**2  + n1*n2*(mu1-mu2)**2
    denom = (n1+n2-1)*(n1+n2)
    return numer / denom
    
def computeSE_sensitive( gene, f1EV, f2EV, f1LNorm, f2LNorm):
    """
    Compute gene-wise standard error estimate
    """
    f1ExonExp = [ [f1EV[i][s:(e+1)]*f1LNorm[i] + EPS \
                   for s,e in gene.exonsI] for i in xrange(len( f1EV))]
    # collapse over replicates
    F1C = numpy.array([ [(f1EV[i][s:(e+1)]*f1LNorm[i]).mean()+EPS \
                                    for s,e in gene.exonsI] for i in xrange(len( f1EV))]).mean(0)
               
    f2ExonExp = [ [f2EV[i][s:(e+1)]*f2LNorm[i] + EPS \
                   for s,e in gene.exonsI] for i in xrange(len( f2EV))]

    # collapse over replicates
    F2C = numpy.array([ [(f2EV[i][s:(e+1)]*f2LNorm[i]).mean()+EPS \
                                    for s,e in gene.exonsI] for i in xrange(len( f2EV))]).mean(0)

    vars1 = numpy.array( [ geoStd( x ) for x in chain.from_iterable( f1ExonExp) ] )
    vars2 = numpy.array( [ geoStd( x ) for x in chain.from_iterable( f2ExonExp) ] )
 
    mu1 = numpy.array( [ gmean( x ) for x in chain.from_iterable( f1ExonExp) ] ).mean()
    mu2 = numpy.array( [ gmean( x ) for x in chain.from_iterable( f2ExonExp) ] ).mean()

    up = max(F1C.sum(), F2C.sum())
    fc = numpy.log2(F1C.mean() / F2C.mean())
    f1Exp = F1C.mean()
    f2Exp = F2C.mean()
    
    f1Norm = up / F1C.sum()
    f2Norm = up / F2C.sum()

    s1 = numpy.sqrt(numpy.array( [x**2 for x in vars1] ).sum() / len( vars1 ))
    s2 = numpy.sqrt(numpy.array( [x**2 for x in vars2] ).sum() / len( vars2 ))

    #return numpy.sqrt(0.5*(geostd1+geostd2) )
    return numpy.sqrt( avg_two_variances(len(vars1), len(vars2), s1, s2, mu1, mu2) ), f1Norm, f2Norm, f1Exp, f2Exp, fc

def computeSE( gene, f1EV, f2EV, f1LNorm, f2LNorm):
    """
    Compute gene-wise standard error estimate
    """
    f1ExonExp = [ [(f1EV[i][s:(e+1)]*f1LNorm[i]).mean()+EPS \
                       for s,e in gene.exonsI] for i in xrange(len( f1EV))]
    f1RepExps = [ numpy.mean(x) for x in f1ExonExp ]
    mexp = max( f1RepExps )
    f1RepExps = [ mexp/x for x in f1RepExps ]
    # collapse over replicates
    F1C = numpy.array([ [(f1EV[i][s:(e+1)]*f1LNorm[i]).mean() \
                                    for s,e in gene.exonsI] for i in xrange(len( f1EV))]).mean(0)+EPS

    F1Cr = numpy.array([ [(f1EV[i][s:(e+1)] \
                           for s,e in gene.exonsI] for i in xrange(len( f1EV))]).mean(0)+EPS

    f2ExonExp = [ [(f2EV[i][s:(e+1)]*f2LNorm[i]).mean()+EPS \
                       for s,e in gene.exonsI] for i in xrange(len( f2EV))]
    f2RepExps = [ numpy.mean(x) for x in f2ExonExp ]
    mexp = max( f2RepExps )
    f2RepExps = [ mexp/x for x in f2RepExps ]

    # collapse over replicates
    F2C = numpy.array([ [(f2EV[i][s:(e+1)]*f2LNorm[i]).mean() \
                                    for s,e in gene.exonsI] for i in xrange(len( f2EV))]).mean(0)+EPS
    F2Cr = numpy.array([ [(f2EV[i][s:(e+1)] \
                                    for s,e in gene.exonsI] for i in xrange(len( f2EV))]).mean(0)+EPS

    F1 = numpy.array(list(chain.from_iterable( [f1ExonExp[i]*f1RepExps[i] for i in xrange(len(f1ExonExp))]))) + EPS
    F2 = numpy.array(list(chain.from_iterable( [f2ExonExp[i]*f2RepExps[i] for i in xrange(len(f2ExonExp))]))) + EPS
    up = max(F1C.sum(), F2C.sum())
    fc = numpy.log2(F1C.mean() / F2C.mean())
    f1Exp = F1C.mean()
    f2Exp = F2C.mean()
    f1Expr = F1Cr.mean()
    f2Expr = F2Cr.mean()
    
    f1Norm = up / F1C.sum()
    f2Norm = up / F2C.sum()
    geomean1 = gmean(F1)
    geomean2 = gmean(F2)
    gene.geneExp = numpy.mean([F1C.mean(), F2C.mean()])
    geostd1 = 2**numpy.sqrt(numpy.sum(numpy.log2(F1/geomean1 )**2 ) / len( F1 ))
    geostd2 = 2**numpy.sqrt(numpy.sum(numpy.log2(F2/geomean2 )**2 ) / len( F2 ))

    #return numpy.sqrt(0.5*(geostd1**2+geostd2**2) ) * numpy.sqrt(2.0/(len(f1EV)+len(f2EV))), f1Norm, f2Norm
                       return numpy.sqrt(0.5*(geostd1**2+geostd2**2) ), f1Norm, f2Norm, f1Exp, f2Exp, f1Expr, f2Expr, fc




def computeSE_old( gene, f1EV, f2EV, f1LNorm, f2LNorm):
    """
    Compute gene-wise standard error estimate
    """
    f1ExonExp = [ [(f1EV[i][s:(e+1)]*f1LNorm[i]).mean() \
                       for s,e in gene.exonsI] for i in xrange(len( f1EV))]
    # collapse over replicates
    F1C = numpy.array([ [(f1EV[i][s:(e+1)]*f1LNorm[i]).mean() \
                                    for s,e in gene.exonsI] for i in xrange(len( f1EV))]).mean(0)+EPS
               
    f2ExonExp = [ [(f2EV[i][s:(e+1)]*f2LNorm[i]).mean() \
                       for s,e in gene.exonsI] for i in xrange(len( f2EV))]
    # collapse over replicates
    F2C = numpy.array([ [(f2EV[i][s:(e+1)]*f2LNorm[i]).mean() \
                                    for s,e in gene.exonsI] for i in xrange(len( f2EV))]).mean(0)+EPS

    F1 = numpy.array(list(chain.from_iterable( [f1ExonExp[i] for i in xrange(len(f1ExonExp))]))) + EPS
    F2 = numpy.array(list(chain.from_iterable( [f2ExonExp[i] for i in xrange(len(f2ExonExp))]))) + EPS
    up = max(F1C.sum(), F2C.sum())
    fc = numpy.log2(F1C.mean() / F2C.mean())
    f1Exp = F1C.mean()
    f2Exp = F2C.mean()
    
    f1Norm = up / F1C.sum()
    f2Norm = up / F2C.sum()
    geomean1 = gmean(F1)
    geomean2 = gmean(F2)
    gene.geneExp = numpy.mean([F1C.mean(), F2C.mean()])
    geostd1 = 2**numpy.sqrt(numpy.sum(numpy.log2(F1/geomean1 )**2 ) / len( F1 ))
    geostd2 = 2**numpy.sqrt(numpy.sum(numpy.log2(F2/geomean2 )**2 ) / len( F2 ))

    #return numpy.sqrt(0.5*(geostd1**2+geostd2**2) ) * numpy.sqrt(2.0/(len(f1EV)+len(f2EV))), f1Norm, f2Norm
    return numpy.sqrt(0.5*(geostd1**2+geostd2**2) ), f1Norm, f2Norm, f1Exp, f2Exp, fc


def zscore( x, mu=0, sigma=1.0 ):
    """Computes the z-score

    """
    return ( x - mu ) / float(sigma)


if __name__ == '__main__':
    raise 'Supporting module only'
