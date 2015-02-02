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
from SpliceGrapher.SpliceGraph       import *
from scipy.stats import t,gmean, hmean, sem, ttest_1samp
from scipy.stats import norm as sNorm
from scipy.cluster.vq import kmeans,vq
from multi_test import bh,qvalues,bonferroni
from itertools import chain
import numpy, os, sys, multi_test
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
    """
    Removes primer bias by smoothing read depths.
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
    """
    Compute library norm factors for each replicate and
    separately for both introns and exons
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
    """
    Compute library norm factors for each replicate and
    separately for both introns and exons
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

def computeStatistics( geneRecords, f1Dict, f2Dict, nspace ):
    """
    Wrapper to compute statistics for all gene records
    """
    mink=nspace.krange[0]
    maxk=nspace.krange[1]
    c=nspace.coverage
    #biasAdj=nspace.biasAdj

    serrs = [ ]
    if False: #biasAdj:
        f1ExonNorm, f2ExonNorm, f1IntronNorm, f2IntronNorm, k = \
            computeNormFactorsAdj( geneRecords, f1Dict, f2Dict, verbose=nspace.verbose )
        
    else:
        f1ExonNorm, f2ExonNorm, f1IntronNorm, f2IntronNorm, k = \
            computeNormFactors( geneRecords, f1Dict, f2Dict, verbose=nspace.verbose )

    if maxk == None: maxk = max(2, k )
    if mink == None: mink = 2
    aVals = range(mink, maxk+1) 
    if nspace.verbose:
        sys.stderr.write('\tComputing differential splicing scores...\n')
    testedGenes = [ ]
    for gene in geneRecords:
        if gene.gid not in f1Dict or gene.gid not in f2Dict: 
            #print gene.gid
            continue

        f1ExpVr = f1Dict[gene.gid]
        f2ExpVr = f2Dict[gene.gid]
        f1ExpV = [numpy.array(x) for x in f1ExpVr]
        f2ExpV = [numpy.array(x) for x in f2ExpVr]
    
        serr, f1Norm,f2Norm, f1exp, f2exp, fc = computeSE( gene, f1ExpVr, 
                                                           f2ExpVr, f1ExonNorm, 
                                                           f2ExonNorm)

        gene.serr = serr
        gene.f1Norm = f1Norm
        gene.f2Norm = f2Norm
        gene.fc = fc
        gene.f1exp = f1exp
        gene.f2exp = f2exp
        # filter genes with no expression 
        if numpy.sum( numpy.array(f1ExpV).mean(0) ) == 0 or \
                numpy.sum( numpy.array(f2ExpV).mean(0) ) == 0:
            continue
        serrs.append( serr )
        intFC = [ ]
        intTested = []
        intfc = [ ]
        intexp = [ ]
        intserrs = []
        intssdiffs = [ ]
        intstat = [] 
        intQvals = [ list() for _ in gene.introns ]
        intZ     = [ list() for _ in gene.introns ]
        intPvals = [ list() for _ in gene.introns ]

        for k,r  in enumerate(gene.introns):
            s,e = r
            ls,le = gene.exons[k]
            rs,re = gene.exons[k+1]

            if False:#biasAdj:
                f1Int = numpy.array([f1ExpV[i][(s+1):(e)] * f1IntronNorm[i] * ((gene.f1ExonWts[i][k] + gene.f1ExonWts[i][k+1]) / 2.0)\
                                         for i in xrange(len(f1IntronNorm))]).mean(0)
                f1Std = numpy.array( [numpy.sum(f1ExpV[i][(s+1):(e)] * f1IntronNorm[i] * ((gene.f1ExonWts[i][k] + gene.f1ExonWts[i][k+1]) / 2.0))\
                                         for i in xrange(len(f1IntronNorm))]).std()
                f2Int = numpy.array([f2ExpV[i][(s+1):(e)] * f2IntronNorm[i] * ((gene.f2ExonWts[i][k] + gene.f2ExonWts[i][k+1]) / 2.0)\
                                         for i in xrange(len(f2IntronNorm))]).mean(0)
                f2Std = numpy.array( [numpy.sum(f2ExpV[i][(s+1):(e)] * f2IntronNorm[i] * ((gene.f2ExonWts[i][k] + gene.f2ExonWts[i][k+1]) / 2.0))\
                                         for i in xrange(len(f2IntronNorm))]).std()

            else:
                f1Int = numpy.array([f1ExpV[i][(s+1):(e)] * f1IntronNorm[i]\
                                         for i in xrange(len(f1IntronNorm))]).mean(0)
                f1Std = numpy.array( [numpy.sum(f1ExpV[i][(s+1):(e)] * f1IntronNorm[i])\
                                         for i in xrange(len(f1IntronNorm))]).std()
                f2Int = numpy.array([f2ExpV[i][(s+1):(e)] * f2IntronNorm[i]\
                                         for i in xrange(len(f2IntronNorm))]).mean(0)
                f2Std = numpy.array( [numpy.sum(f2ExpV[i][(s+1):(e)] * f2IntronNorm[i])\
                                         for i in xrange(len(f2IntronNorm))]).std()

            if False:#biasAdj:
                f1ExonL = numpy.array([f1ExpV[i][ls:le]*f1ExonNorm[i]*gene.f1ExonWts[i][k] \
                                           for i in xrange(len(f1IntronNorm))]).mean(0)
                f1ExonR = numpy.array([f1ExpV[i][rs:re]*f1ExonNorm[i]*gene.f1ExonWts[i][k+1] \
                                           for i in xrange(len(f1IntronNorm))]).mean(0)
                f2ExonL = numpy.array([f2ExpV[i][ls:le]*f2ExonNorm[i]*gene.f2ExonWts[i][k] \
                                           for i in xrange(len(f2IntronNorm))]).mean(0)
                f2ExonR = numpy.array([f2ExpV[i][rs:re]*f2ExonNorm[i]*gene.f2ExonWts[i][k+1] \
                                           for i in xrange(len(f2IntronNorm))]).mean(0)
            else:
                f1ExonL = numpy.array([f1ExpV[i][ls:le]*f1ExonNorm[i] \
                                           for i in xrange(len(f1IntronNorm))]).mean(0)
                f1ExonR = numpy.array([f1ExpV[i][rs:re]*f1ExonNorm[i] \
                                           for i in xrange(len(f1IntronNorm))]).mean(0)
                f2ExonL = numpy.array([f2ExpV[i][ls:le]*f2ExonNorm[i] \
                                           for i in xrange(len(f2IntronNorm))]).mean(0)
                f2ExonR = numpy.array([f2ExpV[i][rs:re]*f2ExonNorm[i] \
                                           for i in xrange(len(f2IntronNorm))]).mean(0)


            intserrs.append( numpy.sqrt( f1Std**2 / float(len(f1ExpV)) + f2Std**2 / float(len(f2ExpV)) + gene.serr**2) )
            #intserrs.append(serr)
            f1Intm = numpy.mean(f1Int)
            f2Intm = numpy.mean(f2Int)
            
            f1SSdiff = 0.5*(f1ExonL.mean()+EPS +f1ExonR.mean()+EPS)
            f2SSdiff = 0.5*(f2ExonL.mean()+EPS +f2ExonR.mean()+EPS)

            f1LSS = numpy.array(list(f1ExonL[-10:]) + list(f1Int[:11]))
            f1RSS = numpy.array(list(f1Int[-10:]) + list(f1ExonR[:11]))

            f2LSS = numpy.array(list(f2ExonL[-10:]) + list(f2Int[:11]))
            f2RSS = numpy.array(list(f2Int[-10:]) + list(f2ExonR[:11]))

            cosL = numpy.dot( f1LSS, f2LSS) / (EPS+numpy.linalg.norm(f1LSS) * numpy.linalg.norm(f2LSS))
            cosR = numpy.dot( f1RSS, f2RSS) / (EPS+numpy.linalg.norm(f1RSS) * numpy.linalg.norm(f2RSS))

            f1exon = 0.5 * ( numpy.mean(f1ExonL) + numpy.mean(f1ExonR) )
            f2exon = 0.5 * ( numpy.mean(f2ExonL) + numpy.mean(f2ExonR) )


            tested = True
            # if amiguity from antisense gene, don't test!
            if gene.overlap[k]:
                tested = False
                #print gene.gid, gene.intronsR[k]
            f1Coverage = numpy.sum(f1Int > 0 )/float(len(f1Int))
            f2Coverage = numpy.sum(f2Int > 0 )/float(len(f2Int))
            if numpy.sum(f1LSS > 0 ) / float(len(f1LSS)) < 1 or numpy.sum(f1RSS > 0 ) / float(len(f1RSS)) < 1 or \
               numpy.sum(f2LSS > 0 ) / float(len(f2LSS)) < 1 or numpy.sum(f2RSS > 0 ) / float(len(f2RSS)) < 1:
                tested = False
                                                                                                                                                                                             
            if f1Coverage < c and f2Coverage < c: 
                tested = False
            
            # check if gene is DE and if sufficient read 
            # depth exists in the lower-expressed condition
            if tested:
                if gene.f1Norm < gene.f2Norm and f1Intm > f2Intm * f2Norm:
                    if  max(f1Norm, f2Norm) > nspace.dexpThresh:
                        tested = False
                    if f2Intm == 0 and (f2exon+EPS) * min(1.0, f1Intm / (f1exon+EPS) ) < 1.0/len(f1Int):
                        tested = False

                elif gene.f2Norm < gene.f1Norm and f2Intm > f1Intm * f1Norm:
                    if  max(f1Norm, f2Norm) > nspace.dexpThresh:
                        tested = False
                    if f1Intm == 0 and (f1exon+EPS)* min(1.0, f2Intm / (f2exon+EPS) ) < 1.0/len(f2Int):            
                        tested = False

            if not tested:
                intTested.append(False)
                intFC.append(None)
                intfc.append(None)
                intexp.append(None)
                intstat.append(None)

                
            else:
                intTested.append(True)
                numer = [numpy.log2( 2**a + f1Intm * f1Norm  ) \
                             for a in aVals ]
                denom = [numpy.log2( 2**a + f2Intm * f2Norm ) \
                             for a in aVals ]
                numerL = numpy.log2( (2**a + f1Int[:11].mean())/(f1ExonL.mean() + EPS) )
                denomL = numpy.log2( (2**a + f2Int[:11].mean())/(f2ExonL.mean()  + EPS) )
                numerR = numpy.log2( (2**a + f1Int[-10:].mean())/(f1ExonR.mean() + EPS) )
                denomR = numpy.log2( (2**a + f2Int[-10:].mean())/(f2ExonR.mean() + EPS) )

                intfc.append( (numpy.log2(f1Intm * f1Norm+EPS)) -  \
                                  (numpy.log2(f2Intm * f2Norm+EPS)))
                intexp.append(0.5 * (numpy.log2(f1Intm * f1Norm+EPS)\
                                         + numpy.log2(f2Intm *\
                                                          f2Norm+EPS)))
                ssDiff =  0.5*( (numerL - denomL) + (numerR-denomR))

                intFC.append( [(numer[i] - denom[i]) / (cosL + cosR) \
                                   for i in xrange(len(numer))])

                intstat.append( [ ((numer[i] - denom[i] ) / (cosL + cosR))/gene.serr \
                                  for i in xrange(len(numer))])
                #intstat.append( [ ((numer[i] - denom[i] ) / (cosL + cosR)) \
                #                  for i in xrange(len(numer))])

                #intstat.append( [ ((numer[i] - denom[i] )  ) / gene.serr \
                #                      for i in xrange(len(numer))])

        if intTested.count(True) > 0:
            gene.intexp = intexp
            gene.intfc  = intfc
            gene.intFC = intFC
            gene.intTested = intTested
            gene.intserrs = intserrs
            gene.intstat = intstat
            gene.intZ = intZ
            gene.intQvals = intQvals
            gene.intPvals = intPvals
            testedGenes.append(gene)
    #compute Z-score distribution parameters
    for aidx in xrange(len(aVals)):
        if nspace.verbose:
            sys.stderr.write('\tComputing significance for alpha=%d...\n'% aVals[aidx])

        X = list(chain.from_iterable([ numpy.array([testedGenes[i].intstat[t][aidx]\
                                                        for t in xrange(len(testedGenes[i].intFC)) \
                                                        if testedGenes[i].intTested[t] ]) \
                                           for i in xrange(len(testedGenes))])) 
        mu = numpy.mean(X)
        sigma = numpy.std(X)


        # Assign z-scores and pvalues 
        for gene in testedGenes:
            for i in xrange( len( gene.intTested)):
                if gene.intTested[i]:
                    z = zscore(gene.intstat[i][aidx], mu, sigma)
                    gene.intZ[i].append( z )
                    gene.intPvals[i].append(min(1.0,2*sNorm.cdf( -abs(z), loc=0, scale=1.0)))

                else:
                    gene.intZ[i].append(0)
                    gene.intPvals[i].append(2)

        # compute q-values
        N =  len(testedGenes)
        pvals = list(chain.from_iterable([ testedGenes[i].intPvals[t][aidx] \
                                               for t in xrange(
                        len(testedGenes[i].intTested)) \
                                               if testedGenes[i].intTested[t]] \
                                               for i in xrange(N)) )
        if nspace.multTest == 'BH':
            qvals = bh(pvals)
        elif nspace.multTest == 'BF':
            qvals = bonferroni(pvals)
        elif nspace.multTest == 'QV':
            qvals = qvalues(pvals)[0]
        #qvals = bh(pvals)
        i = 0
        for gene in testedGenes:
            ci = 0
            tested = gene.intTested.count(True)
            for idx in xrange(len(gene.introns)):
                if gene.intTested[idx]:
                    gene.intQvals[idx].append(qvals[i+ci])
                    ci += 1
                else:
                    gene.intQvals[idx].append(2)
            i = i + tested
    return testedGenes, aVals

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

def computeSEStatistics( geneRecords, f1Dict, f2Dict, nspace ):
    """
    Wrapper to compute statistics for all gene records
    """

    mink=nspace.krange[0]
    maxk=nspace.krange[1]
    c=nspace.coverage
    #biasAdj=nspace.biasAdj

    serrs = [ ]
    f1ExonNorm, f2ExonNorm, f1IntronNorm, f2IntronNorm, k = computeNormFactors( geneRecords, 
                                                                                f1Dict, 
                                                                                f2Dict, 
                                                                                verbose=nspace.verbose )

    if maxk == None: maxk = max(2, k )
    if mink == None: mink = 2
    aVals = range(mink, maxk+1) 
    if nspace.verbose:
        sys.stderr.write('\tComputing differential splicing scores...\n')
    testedGenes = [ ]
    for gene in geneRecords:
        if gene.gid not in f1Dict or gene.gid not in f2Dict: 
            #print gene.gid
            continue

    for gene in geneRecords:
        SEnodes = gene.flavorDict['SE']

        if not SEnodes: continue
        nSEs = len(SEnodes)
        SEQvals = [ list() for _ in gene.flavorDict['SE']]
        SEZ     = [ list() for _ in gene.flavorDict['SE']]
        SEPvals = [ list() for _ in gene.flavorDict['SE']]
        f1ExpVr = f1Dict[gene.gid]
        f2ExpVr = f2Dict[gene.gid]
        f1ExpV = [numpy.array(x) for x in f1ExpVr]
        f2ExpV = [numpy.array(x) for x in f2ExpVr]
        serr, f1Norm,f2Norm, f1exp, f2exp, fc = computeSE( gene, f1ExpVr, 
                                                           f2ExpVr, f1ExonNorm, 
                                                           f2ExonNorm)

        gene.serr = serr
        gene.f1Norm = f1Norm
        gene.f2Norm = f2Norm
        # filter genes with no expression 
        if numpy.sum( numpy.array(f1ExpV).mean(0) ) == 0 or \
                numpy.sum( numpy.array(f2ExpV).mean(0) ) == 0:
            continue

        serrs.append( serr )
        SEFC = [ ]
        SETested = []
        SEfc = [ ]
        SEexp = [ ]
        SEserrs = []
        SEstat = []
        for k,event  in enumerate(SEnodes):
            l,se,r = event
            ls, le = l
            s, e = se
            e += 1
            rs, re = r
            if False:#biasAdj:
                f1SE = numpy.array([f1ExpV[i][(s):(e)] * f1ExonNorm[i] * get_weight(s,e,gene,f1ExonWts[i])\
                                        for i in xrange(len(f1ExonNorm))]).mean(0)
                f1Std = numpy.array( [numpy.sum(f1ExpV[i][(s):(e)] * f1ExonNorm[i] * get_weight(s,e,gene,f1ExonWts[i]))\
                                          for i in xrange(len(f1ExonNorm))]).std()
                f2SE = numpy.array([f2ExpV[i][(s):(e)] * f2ExonNorm[i] * get_weight(s,e,gene,f2ExonWts[i])\
                                        for i in xrange(len(f2ExonNorm))]).mean(0)
                f2Std = numpy.array( [numpy.sum(f2ExpV[i][(s):(e)] * f2ExonNorm[i] * get_weight(s,e,gene,f2ExonWts[i]))\
                                          for i in xrange(len(f2ExonNorm))]).std()

            else:
                f1SE = numpy.array([f1ExpV[i][s:e] * f1ExonNorm[i]\
                                        for i in xrange(len(f1ExonNorm))]).mean(0)
                f1Std = numpy.array( [numpy.sum(f1ExpV[i][s:e] * f1ExonNorm[i])\
                                          for i in xrange(len(f1ExonNorm))]).std()
                f2SE = numpy.array([f2ExpV[i][s:e] * f2ExonNorm[i]\
                                        for i in xrange(len(f2ExonNorm))]).mean(0)
                f2Std = numpy.array( [numpy.sum(f2ExpV[i][s:e] * f2ExonNorm[i])\
                                          for i in xrange(len(f2ExonNorm))]).std()


            if False:#biasAdj:
                f1ExonL = numpy.array([f1ExpV[i][ls:le]*f1ExonNorm[i]*get_weight( ls, le, f1ExonWts[i]) \
                                           for i in xrange(len(f1IntronNorm))]).mean(0)
                f1ExonR = numpy.array([f1ExpV[i][rs:re]*f1ExonNorm[i]*get_weight( rs, re, f1ExonWts[i]) \
                                           for i in xrange(len(f1IntronNorm))]).mean(0)
                f2ExonL = numpy.array([f2ExpV[i][ls:le]*f2ExonNorm[i]*get_weight( ls, le, f2ExonWts[i]) \
                                           for i in xrange(len(f2IntronNorm))]).mean(0)
                f2ExonR = numpy.array([f2ExpV[i][rs:re]*f2ExonNorm[i]*get_weight( rs, re, f2ExonWts[i]) \
                                           for i in xrange(len(f2IntronNorm))]).mean(0)
            else:
                f1ExonL = numpy.array([f1ExpV[i][ls:le]*f1ExonNorm[i] \
                                           for i in xrange(len(f1IntronNorm))]).mean(0)
                f1ExonR = numpy.array([f1ExpV[i][rs:re]*f1ExonNorm[i] \
                                           for i in xrange(len(f1IntronNorm))]).mean(0)
                f2ExonL = numpy.array([f2ExpV[i][ls:le]*f2ExonNorm[i] \
                                           for i in xrange(len(f2IntronNorm))]).mean(0)
                f2ExonR = numpy.array([f2ExpV[i][rs:re]*f2ExonNorm[i] \
                                           for i in xrange(len(f2IntronNorm))]).mean(0)
                
                
            SEserrs.append( numpy.sqrt( f1Std**2 / float(len(f1ExpV)) + f2Std**2 / float(len(f2ExpV)) + gene.serr**2) )
            tested = True
            #intserrs.append(serr)
            if len(f1SE) == 0 :
                tested = False
                print s,e, gene.gid
                f1SEm = EPS
                f2SEm = EPS
            else:
                f1SEm = numpy.mean(f1SE)
                f2SEm = numpy.mean(f2SE)
            

            
            f1SSdiff = 0.5*(f1ExonL.mean()+EPS +f1ExonR.mean()+EPS)
            f2SSdiff = 0.5*(f2ExonL.mean()+EPS +f2ExonR.mean()+EPS)

            f1LSS = numpy.array(list(f1ExonL[-10:]) + list(f1SE[:11]))
            f1RSS = numpy.array(list(f1SE[-10:]) + list(f1ExonR[:11]))

            f2LSS = numpy.array(list(f2ExonL[-10:]) + list(f2SE[:11]))
            f2RSS = numpy.array(list(f2SE[-10:]) + list(f2ExonR[:11]))

            cosL = numpy.dot( f1LSS, f2LSS) / (EPS+numpy.linalg.norm(f1LSS) * numpy.linalg.norm(f2LSS))
            cosR = numpy.dot( f1RSS, f2RSS) / (EPS+numpy.linalg.norm(f1RSS) * numpy.linalg.norm(f2RSS))

            f1exon = 0.5 * ( numpy.mean(f1ExonL) + numpy.mean(f1ExonR) )
            f2exon = 0.5 * ( numpy.mean(f2ExonL) + numpy.mean(f2ExonR) )

            f1Coverage = numpy.sum(f1SE > 0 )/float(len(f1SE))
            f2Coverage = numpy.sum(f2SE > 0 )/float(len(f2SE))

            #if numpy.sum(f1LSS > 0 ) / float(len(f1LSS)) < 1 or numpy.sum(f1RSS > 0 ) / float(len(f1RSS)) < 1 or \
            #   numpy.sum(f2LSS > 0 ) / float(len(f2LSS)) < 1 or numpy.sum(f2RSS > 0 ) / float(len(f2RSS)) < 1:
            #    tested = False
                                                                                                                                                                                             
            if tested:
                if f1Coverage < c and f2Coverage < c:
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

            if len(f1SE)==0 or len(f2SE)==0:
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
                #intstat.append( [ ((numer[i] - denom[i] )  ) / gene.serr \
                #                      for i in xrange(len(numer))])

        if SETested.count(True) > 0:
            testedGenes.append(gene)
            gene.SEexp = SEexp
            gene.SEfc  = SEfc
            gene.SEFC = SEFC
            gene.SETested = SETested
            gene.SEserrs = SEserrs
            gene.SEstat = SEstat
            gene.SEZ = SEZ
            gene.SEPvals = SEPvals
            gene.SEQvals = SEQvals

    #compute Z-score distribution parameters
    for aidx in xrange(len(aVals)):
        X = list(chain.from_iterable([ numpy.array([testedGenes[i].SEstat[t][aidx]\
                                                        for t in xrange(len(testedGenes[i].SEstat)) \
                                                        if testedGenes[i].SETested[t] ]) \
                                           for i in xrange(len(testedGenes))]) )
        mu = numpy.mean(X)
        sigma = numpy.std(X)

        # Assign z-scores and pvalues 
        for gene in testedGenes:
            for i in xrange( len( gene.SETested)):
                if gene.SETested[i]:
                    z = zscore(gene.SEstat[i][aidx], mu, sigma)
                    gene.SEZ[i].append( z )
                    gene.SEPvals[i].append(min(1.0,2*sNorm.cdf( -abs(z), loc=0, scale=1.0)))
                    
                else:
                    gene.SEZ[i].append(0)
                    gene.SEPvals[i].append(2)

        # compute q-values
        N =  len(testedGenes)
        pvals = list(chain.from_iterable([ testedGenes[i].SEPvals[t][aidx] \
                                               for t in xrange(
                        len(testedGenes[i].SETested)) \
                                               if testedGenes[i].SETested[t]] \
                                               for i in xrange(N)) )

        if nspace.multTest == 'BH':
            qvals = bh(pvals)
        elif nspace.multTest == 'BF':
            qvals = bonferroni(pvals)
        elif nspace.multTest == 'QV':
            qvals = qvalues(pvals)[0]

        i = 0
        for gene in testedGenes:
            ci = 0
            tested = gene.SETested.count(True)
            for idx in xrange(len(gene.flavorDict['SE'])):
                if gene.SETested[idx]:
                    gene.SEQvals[idx].append(qvals[i+ci])
                    ci += 1
                else:
                    gene.SEQvals[idx].append(2)
            i = i + tested
    return testedGenes, aVals

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
    
def computeSE_sensitive( gene, f1EV, f2EV, f1ExonNorm, f2ExonNorm):
    """
    Compute gene-wise standard error estimate
    """
    f1ExonExp = [ [f1EV[i][s:(e+1)]*f1ExonNorm[i] + EPS \
                   for s,e in gene.exonsI] for i in xrange(len( f1EV))]
    # collapse over replicates
    F1C = numpy.array([ [(f1EV[i][s:(e+1)]*f1ExonNorm[i]).mean()+EPS \
                                    for s,e in gene.exonsI] for i in xrange(len( f1EV))]).mean(0)
               
    f2ExonExp = [ [f2EV[i][s:(e+1)]*f2ExonNorm[i] + EPS \
                   for s,e in gene.exonsI] for i in xrange(len( f2EV))]

    # collapse over replicates
    F2C = numpy.array([ [(f2EV[i][s:(e+1)]*f2ExonNorm[i]).mean()+EPS \
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

def computeSE( gene, f1EV, f2EV, f1ExonNorm, f2ExonNorm):
    """
    Compute gene-wise standard error estimate
    """
    f1ExonExp = [ [(f1EV[i][s:(e+1)]*f1ExonNorm[i]).mean()+EPS \
                       for s,e in gene.exonsI] for i in xrange(len( f1EV))]
    f1RepExps = [ numpy.mean(x) for x in f1ExonExp ]
    mexp = max( f1RepExps )
    f1RepExps = [ mexp/x for x in f1RepExps ]
    # collapse over replicates
    F1C = numpy.array([ [(f1EV[i][s:(e+1)]*f1ExonNorm[i]).mean() \
                                    for s,e in gene.exonsI] for i in xrange(len( f1EV))]).mean(0)+EPS
               
    f2ExonExp = [ [(f2EV[i][s:(e+1)]*f2ExonNorm[i]).mean()+EPS \
                       for s,e in gene.exonsI] for i in xrange(len( f2EV))]
    f2RepExps = [ numpy.mean(x) for x in f2ExonExp ]
    mexp = max( f2RepExps )
    f2RepExps = [ mexp/x for x in f2RepExps ]

    # collapse over replicates
    F2C = numpy.array([ [(f2EV[i][s:(e+1)]*f2ExonNorm[i]).mean() \
                                    for s,e in gene.exonsI] for i in xrange(len( f2EV))]).mean(0)+EPS

    F1 = numpy.array(list(chain.from_iterable( [f1ExonExp[i]*f1RepExps[i] for i in xrange(len(f1ExonExp))]))) + EPS
    F2 = numpy.array(list(chain.from_iterable( [f2ExonExp[i]*f2RepExps[i] for i in xrange(len(f2ExonExp))]))) + EPS
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

def computeSE_old( gene, f1EV, f2EV, f1ExonNorm, f2ExonNorm):
    """
    Compute gene-wise standard error estimate
    """
    f1ExonExp = [ [(f1EV[i][s:(e+1)]*f1ExonNorm[i]).mean() \
                       for s,e in gene.exonsI] for i in xrange(len( f1EV))]
    # collapse over replicates
    F1C = numpy.array([ [(f1EV[i][s:(e+1)]*f1ExonNorm[i]).mean() \
                                    for s,e in gene.exonsI] for i in xrange(len( f1EV))]).mean(0)+EPS
               
    f2ExonExp = [ [(f2EV[i][s:(e+1)]*f2ExonNorm[i]).mean() \
                       for s,e in gene.exonsI] for i in xrange(len( f2EV))]
    # collapse over replicates
    F2C = numpy.array([ [(f2EV[i][s:(e+1)]*f2ExonNorm[i]).mean() \
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
    """
    Computes the z-score
    """
    return ( x - mu ) / float(sigma)


def loadTest(dat='sr45'):
    """
    Load data from count files
    """
    from SpliceGrapher.formats.fasta import fasta_itr
    # load factor 1 depths
    class Raw(object): pass
    nspace = Raw()
    if dat == 'sr45':
        nspace.factor1files = ['../../diftext/sr45_counts/sr45_d.cnt.gz','../../diftext/sr45_counts/sr45_e.cnt.gz' ]
        nspace.factor2files = ['../../diftext/sr45_counts/wt_d.cnt.gz', '../../diftext/sr45_counts/wt_e.cnt.gz']
    elif dat == 'ptb':
        nspace.factor1files = ['../../diftext/ptb_mi12_r1.cnt.gz','../../diftext/ptb_mi12_r2.cnt.gz' ]
        nspace.factor2files = ['../../diftext/ptb_wt_r1.cnt.gz', '../../diftext/ptb_wt_r2.cnt.gz']
    elif dat =='hnrnp':
        nspace.factor1files = ['../../diftext/mt_hnrnp.cnt.gz']
        nspace.factor2files = ['../../diftext/wt_hnrnp.cnt.gz']


    else:
        nspace.factor1files = ['../../diftext/mt_met1_MATS.gz']
        nspace.factor2files = ['../../diftext/wt_met1_MATS.gz']

    nspace.verbose = True
    f1Dict = { }
    for i in xrange(len(nspace.factor1files)):
        itr = fasta_itr( nspace.factor1files[i] )
        sys.stderr.write("Reading depths from %s\n" % (nspace.factor1files[i]))

        indicator = ProgressIndicator(10000, verbose=nspace.verbose)
        for rec in itr:
            key = rec.header.split(':')[0].strip()
            l = f1Dict.get(key, [ ] )
            l.append( numpy.array( rec.sequence.strip().split(), int ) )
            f1Dict[key] = l
            indicator.update()
        indicator.finish()
        sys.stderr.write("Read %d gene depths\n" % (indicator.ctr))


    # load factor 2 depths
    f2Dict = { }
    for i in xrange(len(nspace.factor2files)):
        repDict = { }
        itr = fasta_itr( nspace.factor2files[i] )
        sys.stderr.write("Reading depths from %s\n" % ( nspace.factor2files[i] ))
        indicator = ProgressIndicator(10000, verbose=nspace.verbose)
        for rec in itr:
            key = rec.header.split(':')[0].strip()
            l = f2Dict.get(key, [ ] )
            l.append( numpy.array( rec.sequence.strip().split(), int ) )
            f2Dict[key] = l
            indicator.update()
        indicator.finish()
        sys.stderr.write("Read %d gene depths\n" % (indicator.ctr))

    return f1Dict, f2Dict
