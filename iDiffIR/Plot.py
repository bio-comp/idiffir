# Copyright (C) 2013 by Colorado State University
# Contact: Michael Hamilton <hamiltom@cs.colostate.edu>
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307,
# USA.
#
import matplotlib
matplotlib.use('agg')
from SpliceGrapher.shared.utils       import *
from SpliceGrapher.plot.PlotterConfig import *
from SpliceGrapher.view.ViewerUtils   import *
from SpliceGrapher.plot.PlotUtils import *
from sys import maxint as MAXINT
from itertools import chain
import os,sys
sys.setrecursionlimit(10000)
from multiprocessing import Pool, freeze_support, Queue, Process
import warnings  
import numpy, ConfigParser
from iDiffIR.BamfileIO import *

#warnings.filterwarnings('ignore')
from matplotlib import rc
rc('text', usetex=True)
import numpy as np
DEFAULT_FONT   = 12
DEFAULT_HEIGHT = 11.0
DEFAULT_WIDTH  = 8.5

def plotGBCs(f1Codes, f2Codes, f1Cnt, f2Cnt, rmi, odir=os.getcwd() ):
    """
    Plot Global Bias Curves (GBCs) for both
    conditions for each replicate.
    """
    L = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M']
    numClusts = len(f1Codes[0])
    f2map = {}#{0:0,1:1,2:2,3:3,4:4 }
    for i in xrange(numClusts):
        mC = None
        mV = 0
        for j in xrange(numClusts):
            cos = numpy.dot( f1Codes[0][i], f2Codes[0][j]) / ( numpy.linalg.norm( f1Codes[0][i] ) *  numpy.linalg.norm( f2Codes[0][j]))
            if cos > mV:
                mC = j
                mV = cos
                
        f2map[i] = mC
    
    plt.cla()
    plt.figure()
    plt.hold(True)
    for i in xrange(numClusts):
        plt.plot(f1Codes[0][i], '--', marker='.', lw=2, label=L[i])
    plt.legend(ncol=5, loc='upper center')
    plt.xlabel('Bins')
    plt.ylabel('Proportion of reads')
    plt.xticks(range(10), range(1,11))
    plt.savefig(os.path.join(odir, 'figures', 'f1codes.pdf') )
    plt.hold(False)
    plt.cla()
    plt.figure()
    plt.hold(True)
    for i in xrange(numClusts):
        plt.plot(f2Codes[0][f2map[i]], '--', marker='.', lw=2,label=L[i])
    plt.legend(ncol=5, loc='upper center')
    plt.xlabel('Bins')
    plt.ylabel('Proportion of reads')
    plt.xticks(range(10), range(1,11))
    plt.savefig(os.path.join( odir, 'figures', 'f2codes.pdf') )
    plt.hold(False)
    plt.close()
    with open( os.path.join(odir, 'lists', 'curveCounts.txt'), 'w') as fout:
        fout.write( 'f1: %s\n' % ( ','.join( [ '%s:%d' % (L[key], f1Cnt[key]) for key in f1Cnt])))
        fout.write( 'f2: %s\n' % ( ','.join( [ '%s:%d' % (L[f2map[key]], f2Cnt[f2map[key]]) for key in f2Cnt])))
        fout.write( "%f\n" % rmi )
                   
def writeAll( geneRecords, aVals, odir=os.getcwd() ):
    """
    Write a file of all intron values 
    """
    # iterate across all genes
    with open( os.path.join(odir, 'allIntrons.txt'), 'w' ) as fout:
        fout.write('\t'.join(['geneID', 'lowExonCoords', 
                               'intronCoords', 'highExonCoords',
                               'pValue', 'adjPValue', 
                               'logFoldChange','intronExp', 'statistic', 
                               'bestA','known']) + '\n')
        for gene in geneRecords:
            # iterate across all exons
            if not gene.IRGTested: continue
            for i,intron in enumerate(gene.introns):
                if not gene.IRTested[i]: continue
                fout.write('%s\t%s\t%s\t%s\t%f\t%f\t%f\t%f\t%f\t%d\t%s\n' % \
                           ( gene.gid, str(gene.exonsR[i]), str(gene.intronsR[i]), str(gene.exonsR[i+1]), 
                             min(1, min(gene.IRPvals[i])),min( 1, min(gene.IRQvals[i])), gene.IRfc[i], 
                             gene.IRexp[i], gene.IRstat[i][numpy.argmin(gene.IRPvals[i])], 
                             aVals[numpy.argmin(gene.IRPvals[i])], str(gene.retained[i]) )) 
                
def writeAllSE( geneRecords, aVals, odir=os.getcwd() ):
    """
    Write a file of all intron values 
    """
    # iterate across all genes
    with open( os.path.join(odir, 'allSEs.txt'), 'w' ) as fout:

        fout.write('\t'.join(['geneID', 'lowExonCoords', 
                               'SECoords', 'highExonCoords',
                               'pValue', 'adjPValue', 
                               'logFoldChange','SEExp', 'statistic', 
                               'bestA']) + '\n')
        for gene in geneRecords:
            if not gene.SEGTested: continue
            # iterate across all exons
            for i,event in enumerate(gene.flavorDict['SE']):
                if not gene.SETested[i]: continue
                fout.write('%s\t%s\t%s\t%s\t%f\t%f\t%f\t%f\t%f\t%d\n' % \
                           ( gene.gid, str(event[0]), str(event[1]), 
                             str(event[2]), min(1, min(gene.SEPvals[i])),
                             min( 1, min(gene.SEQvals[i])), gene.SEfc[i], 
                             gene.SEexp[i], 
                             gene.SEstat[i][numpy.argmin(gene.SEPvals[i])], 
                             aVals[numpy.argmin(gene.SEPvals[i])] )) 
    
def writeLists( summaryDict, odir=os.getcwd() ):
    """
    Write gene and intron lists to file
    """
    # genes all
    with open( os.path.join(odir, 'allDIRGenes.txt'), 'w' ) as fout:
        allDIRGenes = sorted(set.union( *[summaryDict['known']['upgenes'], summaryDict['known']['downgenes'], 
                                   summaryDict['novel']['upgenes'], summaryDict['novel']['downgenes']]))
        for gid in allDIRGenes:
            fout.write('%s\n' % gid)
    # genes up
    with open( os.path.join(odir, 'upDIRGenes.txt'), 'w' ) as fout:
        upDIRGenes = sorted(set.union( *[summaryDict['known']['upgenes'], summaryDict['novel']['upgenes'] ]))
        for gid in upDIRGenes:
            fout.write('%s\n' % gid)
    # genes down
    with open( os.path.join(odir, 'downDIRGenes.txt'), 'w' ) as fout:
        downDIRGenes = sorted(set.union( *[summaryDict['known']['downgenes'], summaryDict['novel']['downgenes'] ]))
        for gid in downDIRGenes:
            fout.write('%s\n' % gid)
 
    # IRs all
    with open( os.path.join(odir, 'allDIRs.txt'), 'w' ) as fout:
        allDIRs = sorted(set.union( *[summaryDict['known']['upirs'], summaryDict['known']['downirs'], 
                                   summaryDict['novel']['upirs'], summaryDict['novel']['downirs']]))
        for gid in allDIRs:
            fout.write('%s\n' % gid)
    # IRs up
    with open( os.path.join(odir, 'upDIRs.txt'), 'w' ) as fout:
        upDIRs = sorted(set.union( *[summaryDict['known']['upirs'], summaryDict['novel']['upirs'] ]))
        for gid in upDIRs:
            fout.write('%s\n' % gid)
    # IRs down
    with open( os.path.join(odir, 'downDIRs.txt'), 'w' ) as fout:
        downDIRs = sorted(set.union( *[summaryDict['known']['downirs'], summaryDict['novel']['downirs'] ]))
        for gid in downDIRs:
            fout.write('%s\n' % gid)
 
def writeGeneExpression( genes, odir ):
    # genes expression
    with open( os.path.join(odir, 'geneExp.txt'), 'w' ) as fout:
        for gene in genes:
            fout.write( '%s\t%f\t%f\t%f\n' % (gene.gid, gene.f1exp, gene.f2exp, gene.fc))
        
def fullTexTable( summaryTable, odir=os.getcwd() ):
    handle = open( os.path.join(odir, 'dirTable.tex'), 'w' )
    ks = sorted(summaryTable.keys())[:-2]
    # upregulated
    handle.write('& Up')
    for key in ks:
        if type(key) != int: continue
        handle.write( " & %d/%d" % (len(summaryTable[key]['known']['upirs']),
                                    len(summaryTable[key]['novel']['upirs']) ))
    k = [summaryTable[key]['known']['upirs'] for key in ks if type(key) == int]
    n = [summaryTable[key]['novel']['upirs'] for key in ks if type(key) == int]
    handle.write( '& %d/%d' % (len( set.union(*k)), len( set.union(*n))) )
    handle.write('\\\\')
    handle.write('\n')

    # downregulated
    handle.write('& Down')
    for key in ks:
        if type(key) != int: continue
        handle.write( " & %d/%d" % (len(summaryTable[key]['known']['downirs']),
                                    len(summaryTable[key]['novel']['downirs']) ))
    k = [summaryTable[key]['known']['downirs'] for key in ks if type(key) == int ]
    n = [summaryTable[key]['novel']['downirs'] for key in ks if type(key) == int ]
    handle.write( '& %d/%d' % (len( set.union(*k)), len( set.union(*n))) )
    handle.write('\\\\')
    handle.write('\n')

    # both
    handle.write('& Tot')
    for key in ks:
        if type(key) != int: continue
        handle.write( " & %d/%d" % (len(summaryTable[key]['known']['downirs'])+len(summaryTable[key]['known']['upirs']),
                                    len(summaryTable[key]['novel']['downirs'])+len(summaryTable[key]['novel']['upirs']) ))

    k = [summaryTable[key]['known']['downirs'] for key in ks if type(key) == int]+[summaryTable[key]['known']['upirs'] for key in ks if type(key) == int]
    n = [summaryTable[key]['novel']['downirs'] for key in ks if type(key) == int]+[summaryTable[key]['novel']['upirs'] for key in ks if type(key) == int]
    handle.write( '& %d/%d' % (len( set.union(*k)), len( set.union(*n))) )
    handle.write('\\\\')
    handle.write('\n')

    # upregulated
    handle.write('& Up')
    for key in ks:
        handle.write( " & %d/%d" % (len(summaryTable[key]['known']['upgenes']),
                                    len(summaryTable[key]['novel']['upgenes']) ))
    k = [summaryTable[key]['known']['upgenes'] for key in ks]
    n = [summaryTable[key]['novel']['upgenes'] for key in ks]
                           
    handle.write( '& %d/%d' % (len( set.union(*k)), len( set.union(*n))) )
    handle.write('\\\\')
    handle.write('\n')

    # downregulated
    handle.write('& Down')
    for key in ks:
        handle.write( " & %d/%d" % (len(summaryTable[key]['known']['downgenes']),
                                    len(summaryTable[key]['novel']['downgenes']) ))
    k = [summaryTable[key]['known']['downgenes'] for key in ks]
    n = [summaryTable[key]['novel']['downgenes'] for key in ks]
    handle.write( '& %d/%d' % (len( set.union(*k)), len( set.union(*n))) )
    handle.write('\\\\')
    handle.write('\n')

    # both
    handle.write('% Tot')
    for key in ks:
        handle.write( " & %d/%d" % (len(summaryTable[key]['known']['downgenes'])+len(summaryTable[key]['known']['upgenes']),
                                    len(summaryTable[key]['novel']['downgenes'])+len(summaryTable[key]['novel']['upgenes']) ))
    k = [summaryTable[key]['known']['downgenes'] for key in ks]+[summaryTable[key]['known']['upgenes'] for key in ks]
    n = [summaryTable[key]['novel']['downgenes'] for key in ks]+[summaryTable[key]['novel']['upgenes'] for key in ks]
    handle.write( '& %d/%d' % (len( set.union(*k)), len( set.union(*n))) )
    handle.write('\\\\')
    handle.write('\n')


def fullTexTableSE( summaryTable, odir=os.getcwd() ):
    handle = open( os.path.join(odir, 'dirTableSE.tex'), 'w' )
    ks = sorted(summaryTable.keys())[:-1]
    # upregulated
    handle.write('& Up')
    for key in ks:
        if type(key) != int: continue
        handle.write( " & %d" % (len(summaryTable[key]['upirs'])))
                            
    k = [summaryTable[key]['upirs'] for key in ks if type(key) == int]

    handle.write( '& %d' % (len( set.union(*k))))
    handle.write('\\\\')
    handle.write('\n')

    # downregulated
    handle.write('& Down')
    for key in ks:
        if type(key) != int: continue
        handle.write( " & %d" % (len(summaryTable[key]['downirs'])))

    k = [summaryTable[key]['downirs'] for key in ks if type(key) == int ]
    handle.write( '& %d' % (len( set.union(*k))))
    handle.write('\\\\')
    handle.write('\n')

    # both
    handle.write('& Tot')
    for key in ks:
        if type(key) != int: continue
        handle.write( " & %d" % (len(summaryTable[key]['downirs'])+len(summaryTable[key]['upirs'])))


    k = [summaryTable[key]['downirs'] for key in ks if type(key) == int]+[summaryTable[key]['upirs'] for key in ks if type(key) == int]

    handle.write( '& %d' % (len( set.union(*k))))
    handle.write('\\\\')
    handle.write('\n')


    
def summary(geneRecords, aVals, level=0.05):
    summaryTable = {}
    for aidx in xrange(len(aVals)):
        novel = {'upirs':set(), 'downirs':set(),
                 'upgenes':set(), 'downgenes':set()}
        known = {'upirs':set(), 'downirs':set(), 
                 'upgenes':set(), 'downgenes':set()}
        for gene in geneRecords:
            minGene = 1
            if not gene.IRGTested : continue
            for i, qvals in enumerate(gene.IRQvals):
                p = gene.IRPvals[i][aidx]
                if p < minGene:
                    minGene = p
                # significant
                if gene.IRQvals[i][aidx] < level:
                    # upregulated
                    intID = '%s_%s' % ( gene.gid, i)
                    if gene.IRfc[i] > 0:
                        # known
                        if gene.retained[i]: 
                            known['upirs'].add(intID)
                            known['upgenes'].add(gene.gid)
                        # novel
                        else:
                            novel['upirs'].add(intID)
                            novel['upgenes'].add(gene.gid)
                    # downregulated
                    else:
                        # known
                        if gene.retained[i]: 
                            known['downirs'].add(intID)
                            known['downgenes'].add(gene.gid)
                                
                        # novel
                        else:
                            novel['downirs'].add(intID)
                            novel['downgenes'].add(gene.gid)
            gene.minGene = minGene
        summaryTable[aVals[aidx]] = {'novel':novel, 'known':known}
    summaryTable['novel'] = {}
    summaryTable['novel']['upgenes'] = set.union( *[ summaryTable[i]['novel']['upgenes'] for i in aVals])
    summaryTable['novel']['downgenes'] = set.union( *[ summaryTable[i]['novel']['downgenes'] for i in aVals])
    summaryTable['novel']['upirs'] = set.union( *[ summaryTable[i]['novel']['upirs'] for i in aVals])
    summaryTable['novel']['downirs'] = set.union( *[ summaryTable[i]['novel']['downirs'] for i in aVals])
    summaryTable['known'] = {}
    summaryTable['known']['upgenes'] = set.union( *[ summaryTable[i]['known']['upgenes'] for i in aVals])
    summaryTable['known']['downgenes'] = set.union( *[ summaryTable[i]['known']['downgenes'] for i in aVals])
    summaryTable['known']['upirs'] = set.union( *[ summaryTable[i]['known']['upirs'] for i in aVals])
    summaryTable['known']['downirs'] = set.union( *[ summaryTable[i]['known']['downirs'] for i in aVals])


    return summaryTable

def summarySE(geneRecords, aVals, level=0.05):
    summaryTable = {}
    for aidx in xrange(len(aVals)):
        known = {'upirs':set(), 'downirs':set(), 
                 'upgenes':set(), 'downgenes':set()}
        for gene in geneRecords:
            minGene = 1
            if not gene.SEGTested: continue
            for i, qvals in enumerate(gene.SEQvals):
                p = gene.SEPvals[i][aidx]
                if p < minGene:
                    minGene = p
                # significant
                if gene.SEQvals[i][aidx] < level:
                    # upregulated
                    SEID = '%s_%s' % ( gene.gid, i)
                    if gene.SEfc[i] > 0:
                        known['upirs'].add(SEID)
                        known['upgenes'].add(gene.gid)

                        # downregulated
                    else:
                        known['downirs'].add(SEID)
                        known['downgenes'].add(gene.gid)
            gene.minGene = minGene
        summaryTable[aVals[aidx]] = known

    summaryTable['known'] = {}
    summaryTable['known']['upgenes'] = set.union( *[ summaryTable[i]['upgenes'] for i in aVals])
    summaryTable['known']['downgenes'] = set.union( *[ summaryTable[i]['downgenes'] for i in aVals])
    summaryTable['known']['upirs'] = set.union( *[ summaryTable[i]['upirs'] for i in aVals])
    summaryTable['known']['downirs'] = set.union( *[ summaryTable[i]['downirs'] for i in aVals])

    return summaryTable

def plotResults(geneRecords, labels, nspace, geneModel, useLog=True, odir=os.getcwd()):
     for gene in geneRecords:
         if not gene.IRGTested: continue
         plotme = False
         highlights = []
         for i in xrange(len(gene.introns)):
             s,e = gene.introns[i]
             if min(gene.IRQvals[i]) < 0.05:
                 highlights.append( (s,e) )
                 plotme = True
         if plotme:
             f1Depths, f2Depths, f1Juncs, f2Juncs =  getDepthsFromBamfiles( gene, 
                                                                            nspace.factor1bamfiles, 
                                                                            nspace.factor2bamfiles 
                                                                        )
             depths = f1Depths.tolist() + f2Depths.tolist()
             plotDepth( gene, depths, labels,
                        highlights,  
                        os.path.join(odir,gene.gid+'.pdf'), useLog, 
                        geneModel, nspace.shrink_introns )

def plotResultsSE(geneRecords, labels, nspace, geneModel, useLog=True, odir=os.getcwd()):
     for gene in geneRecords:
         if not gene.SEGTested: continue
         plotme = False
         highlights = []
         for i in xrange(len(gene.flavorDict['SE'])):
             s,e = gene.flavorDict['SE'][i][1] 
             if min(gene.SEQvals[i]) < nspace.fdrlevel:
                 highlights.append( (s,e-1) )
                 plotme = True
                     
         if plotme:
             f1Depths, f2Depths, f1Juncs, f2Juncs =  getDepthsFromBamfiles( gene, 
                                                                            nspace.factor1bamfiles, 
                                                                            nspace.factor2bamfiles 
                                                                        )
             depths = f1Depths.tolist() + f2Depths.tolist()
             plotDepth( gene, depths, labels,
                        highlights,  
                        os.path.join(odir,gene.gid+'.pdf'), useLog, 
                        geneModel, nspace.shrink_introns, hcolor="#00EE00" )



def setValue(argv, key, cfgValue, default) :
    """Simple method that performs the logic for assigning a value to a
    parameter.  If there is no configuration value, returns the default.
    If there is a configuration value, the result depends on whether the
    user entered something on the command line."""
    if not cfgValue : return default
    return default if key in argv else cfgValue

def plotDepth_old( geneM, depths, depthIDs, highlights, outname, de, geneModel, shrink_introns=False, hcolor='red'):
    """
    Plot depth for a gene, highlighting differential splicing events
    """
    gene = geneM.gene
    verbose=False
    rc('text', usetex=False)
    minPos = geneM.origGraph.minpos
    maxPos = geneM.origGraph.maxpos
    strand = gene.strand
    maxY = max( [ max(x) for x in depths])

    plt.figure(frameon=False)
    X = range( minPos, maxPos+1)    
    titlePadding = getTitlePadding(16)
    topLine      = 0.99-titlePadding
    c = len(depths)
    if geneM.predicted: c += 1

    height         = topLine * (1.0/(1+c) - titlePadding)
    patchDict = {}
    # Plot gene model 
    curAxes        = axes([AXIS_LEFT, topLine-height, AXIS_WIDTH, height])
    topLine        = topLine - height - titlePadding 

    # plot annotated graph
    graph = makeSpliceGraph(geneM.gene)
    graph.annotate()
    GeneView(geneM.gene, curAxes).plot()
    SpliceGraphView(graph, curAxes, xLimits=(minPos, maxPos)).plot()
    xlim = curAxes.get_xlim()
    xticks = setXticks(int(min(xlim)), int(max(xlim)))
    if shrink_introns:
        ranges = getGraphRanges(geneModelToSpliceGraph(gene), 
                                scaleFactor=plotConfig.shrink_factor)
        introns = [(adjustPosition(x[0]+gene.minpos,ranges)-gene.minpos+1, 
                    adjustPosition(x[1]+gene.minpos,ranges)-gene.minpos+1)\
                        for x in geneM.introns]
        highlights = [ ( adjustPosition(x[0]+gene.minpos,ranges)-gene.minpos+1, 
                    adjustPosition(x[1]+gene.minpos,ranges)-gene.minpos+1)\
                        for x in highlights]
        for r in ranges:
            r.minpos -= gene.minpos
            r.maxpos -= gene.minpos
    else:
        introns = [(x[0], x[1]) for x in geneM.introns]
        exons = geneM.exons

    if geneM.predicted:
        curAxes        = axes([AXIS_LEFT, topLine-height, AXIS_WIDTH, height])
        topLine        = topLine - height - titlePadding         
        GeneView(gene, curAxes).plot()
        SpliceGraphView(graph, curAxes, xLimits=(minPos, maxPos)).plot()

    #plot depths
    for i in xrange( len(depths) ):
        curAxes        = axes([AXIS_LEFT, topLine-height, AXIS_WIDTH, height])
        GeneView(gene, curAxes).plot()
        curAxes.set_ylim((0,maxY))
        if de or shrink_introns :
            curAxes.set_yscale('log',nonposy='clip')
            curAxes.set_ylim((.3,maxY))
        depth = adjustDepths(depths[i], ranges) if shrink_introns else depths[i]
        curAxes.fill_between( X,  0, depth,
                              facecolor='grey', edgecolor='grey',
                              alpha=0.75)
        for intron in introns:
            fcol = '0.65'
            s,e = intron
            curAxes.fill_between( X[s:e],  0, depth[s:e], 
                                  facecolor=fcol, edgecolor=fcol,alpha=1)
        for region in highlights:
            s,e = region
            curAxes.fill_between( X[s:e],  0, depth[s:e], 
                                  facecolor=hcolor, edgecolor=hcolor,alpha=1)
        for s,e in geneM.exonsI:
            curAxes.fill_between( X[s:e+1],  0, depth[s:e+1], 
                                  facecolor='0.45', edgecolor='0.45',alpha=1)
            
        curAxes.set_xticks(xticks)
        curAxes.set_xticklabels([])
        topLine        = topLine - height - titlePadding
        curAxes.set_xlim(xlim)
        
        curAxes.set_title(depthIDs[i])
    #plotLegend(patchDict)                    
    curAxes.set_xticklabels(xticks)
    sys.stderr.write( 'plotting '+ outname)
    savefig(outname, dpi=400)
    plt.close()

def plotDepth( geneM, depths, depthIDs, highlights, outname, de, geneModel, shrink_introns=False, hcolor='red'):
    """
    Plot depth for a gene, highlighting IR events
    """
    verbose=False
    rc('text', usetex=False)
    gene = geneM.gene
    minPos = geneM.minpos
    maxPos = geneM.maxpos
    X = range( minPos, maxPos+1)    

    strand = gene.strand
    maxY = max( [ max(x) for x in depths])

    # Main setup
    cf = ConfigParser.ConfigParser()
    cf.add_section( 'SinglePlotConfig' )
    cf.set( 'SinglePlotConfig', 'legend', 'True' ) 
    cf.set( 'SinglePlotConfig', 'output_file', outname)
    cf.set( 'SinglePlotConfig', 'height', '8.0' )
    cf.set( 'SinglePlotConfig', 'width', '16.0' )
    cf.set( 'SinglePlotConfig', 'fontsize', '16' )
    cf.set( 'SinglePlotConfig', 'shrink_introns', str(shrink_introns) )

    
    cf.add_section('GeneModel')
    cf.set('GeneModel', 'plot_type','gene')
    cf.set('GeneModel', 'hide', 'False')
    cf.set('GeneModel', 'file_format', 'gene_model')
    cf.set('GeneModel', 'gene_name', geneM.gid )
    cf.set('GeneModel', 'source_file', 'none.gtf')
    cf.set('GeneModel', 'title_string', 'Gene Model for %s' % geneM.gid )

    config = SinglePlotConfig()
    config.config = cf
    config.validate()
    config.instantiate()
    displayList = config.getPlotList()
    plotConfig = config.getConfiguration()
    matplotlib.use('agg')

    displayList[0].source_file = 'none.gtf'
    plt.figure(frameon=False)
    
    geneSpecs = loadGeneData(displayList, 
                             plotConfig, models={'none.gtf':geneModel},
                             verbose=verbose)

    width    = setValue(sys.argv, '-W', plotConfig.width, DEFAULT_WIDTH)
    height   = setValue(sys.argv, '-H', plotConfig.height, DEFAULT_HEIGHT)
    fontsize = setValue(sys.argv, '-F', plotConfig.fontsize, DEFAULT_FONT)
    geneSpecs.setFontSize(fontsize)

    initializePlots(geneSpecs, width, height, displayList)
    titlePadding = getTitlePadding(fontsize)
    topLine      = 0.99-titlePadding
    c = 2 if geneM.predicted else 1
    height         = topLine * (1.0/(len(depths)+c) - titlePadding)
    patchDict = {}
    # Plot gene model 
    curAxes        = axes([AXIS_LEFT, topLine-height, AXIS_WIDTH, height])
    topLine        = topLine - height - titlePadding 
    patches,ignore = generatePlot(displayList[0], geneSpecs, curAxes, 
                                  shrink=plotConfig.shrink_introns, 
                                  xTickLabels=(False),
                                  verbose=True)
    patchDict.update(patches)
    curAxes.set_yticks([])
    xlim = curAxes.get_xlim()
    xticks = setXticks(int(min(xlim)), int(max(xlim)))

    curAxes.set_xticklabels([])
    
    # plot prediction
    if geneM.predicted:
        curAxes        = axes([AXIS_LEFT, topLine-height, AXIS_WIDTH, height])
        topLine        = topLine - height - titlePadding         
        GeneView(geneSpecs.gene, curAxes).plot()
        patches,extra = plotSpliceGraph(geneM.origGraph, curAxes,
                                        geneName=geneM.gid,
                                        genes=[],
                                        title='Predicted Graph',
                                        xLimits=[geneSpecs.gene.minpos, geneSpecs.gene.maxpos],
                                        minwidth=displayList[0].edge_weight,
                                        xLabels=displayList[0].x_labels
                                        )
        patchDict.update(patches)

        curAxes.set_yticks([])
        curAxes.set_xlim(xlim)
    #xticks = setXticks(int(min(xlim)), int(max(xlim)))

    curAxes.set_xticklabels([])

    if shrink_introns:
        ranges = getGraphRanges(geneModelToSpliceGraph(gene), 
                                scaleFactor=plotConfig.shrink_factor)
        introns = [(adjustPosition(x[0]+gene.minpos,ranges)-gene.minpos+1, 
                    adjustPosition(x[1]+gene.minpos,ranges)-gene.minpos+1)\
                        for x in geneM.introns]
        highlights = [ ( adjustPosition(x[0]+gene.minpos,ranges)-gene.minpos+1, 
                    adjustPosition(x[1]+gene.minpos,ranges)-gene.minpos+1)\
                        for x in highlights]
        for r in ranges:
            r.minpos -= gene.minpos
            r.maxpos -= gene.minpos
    else:
        introns = [(x[0], x[1]) for x in geneM.introns]
        exons = geneM.exons

    #plot depths
    for i in xrange( len(depths) ):
        curAxes        = axes([AXIS_LEFT, topLine-height, AXIS_WIDTH, height])
        GeneView(geneSpecs.gene, curAxes).plot()
        curAxes.set_ylim((0,maxY))
        if de:
            curAxes.set_yscale('log',nonposy='clip')
            curAxes.set_ylim((.3,maxY))
        depth = adjustDepths(depths[i], ranges) if shrink_introns else depths[i]

#        for s,e in exons:
#        curAxes.fill_between( X[s:e+1],  0, depth[s:e+1], 
#                                  facecolor='grey', alpha=0.5)
        curAxes.fill_between( X,  0, depth,
                              facecolor='grey', edgecolor='grey',
                              alpha=0.75)
        for intron in introns:
            fcol = '0.65'
            s,e = intron
            curAxes.fill_between( X[s:e],  0, depth[s:e], 
                                  facecolor=fcol, edgecolor=fcol,alpha=1)
        for region in highlights:
            s,e = region
            curAxes.fill_between( X[s:e],  0, depth[s:e], 
                                  facecolor=hcolor, edgecolor=hcolor,alpha=1)
        for s,e in geneM.exonsI:
            curAxes.fill_between( X[s:e+1],  0, depth[s:e+1], 
                                  facecolor='0.45', edgecolor='0.45',alpha=1)
            
        curAxes.set_xticks(xticks)
        curAxes.set_xticklabels([])
        topLine        = topLine - height - titlePadding
        curAxes.set_xlim(xlim)
        
        curAxes.set_title(depthIDs[i])
    #plotLegend(patchDict)                    
    curAxes.set_xticklabels(xticks)
    savefig(outname, dpi=400)
    plt.close()

def plotReducedGraph(geneM, odir=os.getcwd(), ext='png', **args):
    """
    Plot original and reduced splice graphs
    """
    reduced_graph = geneM.graph
    original_grapn = geneModelToSpliceGraph(geneM.gene)

    verbose = args.get( 'verbose', False )

    rc('text', usetex=False)
    gene = geneM.gene
    gid = geneM.gid
    minPos = min( gene.start(), gene.end())
    maxPos = max( gene.start(), gene.end())
    strand = gene.strand

    # Main setup
    cf = ConfigParser.ConfigParser()
    cf.add_section( 'SinglePlotConfig' )
    cf.set( 'SinglePlotConfig', 'legend', 'True' ) 
    cf.set( 'SinglePlotConfig', 'output_file', '')
    cf.set( 'SinglePlotConfig', 'height', '8.0' )
    cf.set( 'SinglePlotConfig', 'width', '16.0' )
    cf.set( 'SinglePlotConfig', 'fontsize', '16' )
    cf.set( 'SinglePlotConfig', 'shrink_introns', 'False' )

    # Gene model
    cf.add_section('GeneModel')
    cf.set('GeneModel', 'plot_type','gene')
    cf.set('GeneModel', 'hide', 'False')
    cf.set('GeneModel', 'file_format', 'gene_model')
    cf.set('GeneModel', 'gene_name', gid )
    cf.set('GeneModel', 'source_file', 'none.gtf')
    cf.set('GeneModel', 'title_string', 'Gene Model for %s' % gid )

    # Reduced Graph
#    cf.add_section('ReducedGraph')
#    cf.set('ReducedGraph', 'plot_type', 'splice_graph')
#    cf.set('ReducedGraph', 'hide', 'False')
#    cf.set('ReducedGraph', 'file_format', 'splice_graph')
#    cf.set('ReducedGraph', 'gene_name', gid)
#    cf.set('ReducedGraph', 'source_file', 'none2.gtf')
#    cf.set('ReducedGraph', 'title_string', 'Reduced Graph')
    
    config = SinglePlotConfig()
    config.config = cf
    config.validate()
    config.instantiate()
    displayList = config.getPlotList()
    plotConfig = config.getConfiguration()
    matplotlib.use('agg')

    displayList[0].source_file = 'none.gtf'

    plt.figure(frameon=False)
    

    geneSpecs = loadGeneData(displayList, 
                             plotConfig, models={'none.gtf':geneModel, 'none2.gtf':reduced_graph},
                             verbose=verbose)

    width    = setValue(sys.argv, '-W', plotConfig.width, DEFAULT_WIDTH)
    height   = setValue(sys.argv, '-H', plotConfig.height, DEFAULT_HEIGHT)
    fontsize = setValue(sys.argv, '-F', plotConfig.fontsize, DEFAULT_FONT)
    geneSpecs.setFontSize(fontsize)

    initializePlots(geneSpecs, width, height, displayList)
    titlePadding = getTitlePadding(fontsize)
    topLine      = 0.99-titlePadding
    height         = topLine * (1.0/2 - titlePadding)
    patchDict = {}
    # Plot gene model
    curAxes        = axes([AXIS_LEFT, topLine-height, AXIS_WIDTH, height])
    topLine        = topLine - height - titlePadding 
    geneSpecs.gene=geneM.ogene
    patches,ignore = generatePlot(displayList[0], geneSpecs, curAxes, 
                                  shrink=plotConfig.shrink_introns, 
                                  xTickLabels=(False),
                                  verbose=True)
    patchDict.update(patches)
    curAxes.set_yticks([])
    xlim = curAxes.get_xlim()
    xticks = setXticks(int(min(xlim)), int(max(xlim)))

    curAxes.set_xticklabels([])
    availHeight  = 0.45*getAvailableHeight(displayList, fontsize, legend=True)
    height         = availHeight * displayList[0].relative_size

    # Plot reduced graph
    curAxes        = axes([AXIS_LEFT, topLine-height, AXIS_WIDTH, height])
    GeneView(geneM.ogene, curAxes).plot()
    patches,extra = plotSpliceGraph(reduced_graph, curAxes,
                                    geneName=gid,
                                    genes=[],
                                    title='Reduced Graph',
                                    xLimits=[geneSpecs.gene.minpos, geneSpecs.gene.maxpos],
                                    minwidth=displayList[0].edge_weight,
                                    xLabels=displayList[0].x_labels
                                    )
    patchDict.update(patches)
    plotLegend(patchDict)                    
    curAxes.set_yticks([])
    curAxes.set_xlim(xlim)
    #xticks = setXticks(int(min(xlim)), int(max(xlim)))

    curAxes.set_xticklabels([])
    outname=os.path.join( odir, gid+'.%s'%ext)
    savefig(outname, dpi=400)
    plt.close()



def plotStatistic(geneRecords, nbins=20, odir=os.getcwd(), ext='pdf'):
    X = list(chain.from_iterable([ numpy.array(geneRecords[i].intFC)/ \
                                       geneRecords[i].serr  \
                                       for i in xrange(len(geneRecords))]) )
    fig = plt.figure()
    plot = fig.add_subplot(111)
    plot.tick_params(axis='both', which='major', labelsize=16)
    plot.tick_params(axis='both', which='minor', labelsize=16)
    
    plt.hist( X, bins=nbins )
    plt.savefig(os.path.join(odir,'figures','stat.%s') % (ext) )
    plt.close()


def plotPDist(geneRecords, odir=os.getcwd(), nbins=20, ext='pdf' ):
    """
    Plot p-value distribution
    """
    pvals = []
    for gene in geneRecords:
        if not gene.IRGTested: continue
        for i, IRPvals in enumerate(gene.IRPvals):
            if gene.IRTested[i]:
                pvals.append( numpy.min(IRPvals))

    fig = plt.figure()
    plot = fig.add_subplot(111)
    plot.tick_params(axis='both', which='major', labelsize=16)
    plot.tick_params(axis='both', which='minor', labelsize=16)
    plt.hist( pvals, bins=nbins )
    plt.savefig(os.path.join(odir,'pvalues.%s') % (ext) )
    plt.close()

def plotPDistSE(geneRecords, odir=os.getcwd(), nbins=20, ext='pdf' ):
    """
    Plot p-value distribution
    """
    pvals = []
    for gene in geneRecords:
        if not gene.SEGTested: continue
        for i, SEPvals in enumerate(gene.SEPvals):
            if gene.SETested[i]:
                pvals.append( numpy.min(SEPvals))

    fig = plt.figure()
    plot = fig.add_subplot(111)
    plot.tick_params(axis='both', which='major', labelsize=16)
    plot.tick_params(axis='both', which='minor', labelsize=16)
    plt.hist( pvals, bins=nbins )
    plt.savefig(os.path.join(odir,'pvaluesSE.%s') % (ext) )
    plt.close()


def plotMVAold(geneRecords, odir=os.getcwd(), ext='pdf'):
    """
    Render MA plot
    """
    fig = plt.figure(figsize=(8,8), dpi=200)

    plot = fig.add_subplot(111)
    plot.tick_params(axis='both', which='major', labelsize=12)
    plot.tick_params(axis='both', which='minor', labelsize=12)

    pvals = list(chain.from_iterable([numpy.array([ geneRecords[i].intPvals[t] for t in xrange(
                            len(geneRecords[i].intPvals)) \
                                                        if geneRecords[i].intTested[t]])\
                                          for i in xrange(len(geneRecords))])) 

    N = len(pvals)
    M = list(chain.from_iterable([numpy.array([ geneRecords[i].intFCr[t] for t in xrange(
                            len(geneRecords[i].intPvals)) \
                                                        if geneRecords[i].intTested[t]])\
                                          for i in xrange(len(geneRecords))])) 
    
    A = list(chain.from_iterable([numpy.array([ geneRecords[i].intExp[t] for t in xrange(
                            len(geneRecords[i].intPvals)) \
                                                        if geneRecords[i].intTested[t]])\
                                          for i in xrange(len(geneRecords))])) 

    #plot background points
    CA = [A[i] for i in xrange(0,N) if pvals[i] >= 0.05 and abs(M[i]) < 1]    
    CM = [M[i] for i in xrange(0,N) if pvals[i] >= 0.05 and abs(M[i]) < 1]
    plt.plot( CA, CM, color='0.65', marker='.', linestyle='')

    CA = [A[i] for i in xrange(0,N) if pvals[i] >= 0.05 and abs(M[i]) >= 1]    
    CM = [M[i] for i in xrange(0,N) if pvals[i] >= 0.05 and abs(M[i]) >= 1]
    plt.plot( CA, CM, color='0.65', marker='.', linestyle='')


    #plot pvalue 0.05  points
    CA = [A[i] for i in xrange(N) if pvals[i] < 0.05 and pvals[i] >= 0.001]
    CM = [M[i] for i in xrange(N) if pvals[i] < 0.05 and pvals[i] >= 0.001]
    plt.plot( CA, CM, 'k.', label='0.05')

    #plot pvalue 0.05  points
    CA = [A[i] for i in xrange(N) if pvals[i] < 0.001 ]
    CM = [M[i] for i in xrange(N) if pvals[i] < 0.001 ]

    plt.plot( CA, CM, 'r.', label='0.05')
    
    plt.xlabel( r'$\frac{1}{2}(\log_2 \hat{f}_1 + \log_2 \hat{f}_2)$', size=18)
    plt.ylabel( r'$\log_{2}\hat{f}_1 - \log_2\hat{f}_2$', size=18)
    plt.grid()
    plt.legend()
    plt.savefig(os.path.join(odir,'figures','mva.%s') % (ext))


def plotMVA(geneRecords, aVals, odir=os.getcwd(), ext='pdf'):
    """
    Render MA plot
    """
    plt.clf()
    plt.cla()
    plt.autoscale(True)
    colors = ['blue', 'green',  'brown','cyan', 'orange', 'olive', 'pink', 'yellow', 'black',
              'SpringGreen', 'Coral']
    fig = plt.figure()
    plot = fig.add_subplot(111)
    plot.tick_params(axis='both', which='major', labelsize=12)
    plot.tick_params(axis='both', which='minor', labelsize=12)
 
    #plot background surface
    nsigfcs = []
    nsigexps = []

    for gene in geneRecords:
        if not gene.IRGTested: continue
        for i, qvals in enumerate(gene.IRQvals):
            if  gene.IRTested[i]and numpy.min(qvals) >= 0.0:
                nsigfcs.append( gene.IRfc[i])
                nsigexps.append( gene.IRexp[i])

    H, xedges, yedges = np.histogram2d(nsigfcs, nsigexps, bins=(20, 10))
    ydiff=0#abs(yedges[0]-yedges[1])
    xdiff=0#abs(xedges[0]-xedges[1])
    
    extent = [min(nsigexps)+ydiff, max(nsigexps)+ydiff, -max(nsigfcs), -min(nsigfcs)]
    
    plt.imshow(H, extent=extent, interpolation='gaussian', origin='upper', cmap=get_cmap('binary'), alpha=0.5, 
               rasterized=0, vmax=np.max(H)/25.0, aspect='equal')
    plt.autoscale(False)
    plt.axhline(0, ls='--', lw=1, color='r')
    for aidx in xrange(len(aVals)):
        M = []
        A = []
        for gene in geneRecords:
            if not gene.IRGTested: continue
            for i, qvals in enumerate(gene.IRQvals):
                if gene.IRTested[i] and numpy.argmin(qvals) == aidx and numpy.min(qvals) < 0.05:
                    M.append( gene.IRfc[i])
                    A.append( gene.IRexp[i])
        
        
        plt.plot( A, M, color=colors[aidx], marker='.', linestyle='', label=r'$a = 2^{%s}$' % aVals[aidx] )
    
    
    plt.xlabel( r'$\frac{1}{2}(\log_2 \hat{f}_1 + \log_2 \hat{f}_2)$', size=18)
    plt.ylabel( r'$\log_{2}\hat{f}_1 - \log_2\hat{f}_2$', size=18)

    # Shrink current axis by 20%
    box = plot.get_position()
    plot.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    plot.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=True)
    plot.grid()
    plt.savefig(os.path.join(odir,'mva.%s') % (ext))
    plt.autoscale(True)
    plt.close()

def plotMVASE(geneRecords, aVals, odir=os.getcwd(), ext='pdf'):
    """
    Render MA plot
    """
    fig = plt.figure(figsize=(8,8), dpi=200)

    plot = fig.add_subplot(111)
    plot.tick_params(axis='both', which='major', labelsize=12)
    plot.tick_params(axis='both', which='minor', labelsize=12)
    colors = ['blue', 'green',  'brown','cyan', 'orange', 'olive', 'pink', 'yellow', 'black',
              'SpringGreen', 'Coral']

    pvals = list(chain.from_iterable([numpy.array([ geneRecords[i].SEPvals[t] for t in xrange(
                            len(geneRecords[i].SEPvals)) \
                                                        if geneRecords[i].SETested[t]])\
                                          for i in xrange(len(geneRecords)) if geneRecords[i].SEGTested])) 

    N = len(pvals)
    M = list(chain.from_iterable([numpy.array([ geneRecords[i].SEfc[t] for t in xrange(
                            len(geneRecords[i].SEPvals)) \
                                                        if geneRecords[i].SETested[t]])\
                                          for i in xrange(len(geneRecords)) if geneRecords[i].SEGTested])) 
    
    A = list(chain.from_iterable([numpy.array([ geneRecords[i].SEexp[t] for t in xrange(
                            len(geneRecords[i].SEPvals)) \
                                                        if geneRecords[i].SETested[t]])\
                                          for i in xrange(len(geneRecords)) if geneRecords[i].SEGTested])) 

    #plot background points
    plt.plot( A, M, color='0.65', marker='.', linestyle='')

    # plot significat SEs
    for aidx in xrange(len(aVals)):
        M = []
        A = []
        for gene in geneRecords:
            if not gene.SEGTested: continue
            for i, qvals in enumerate(gene.SEQvals):
                if gene.SETested[i] and numpy.argmin(qvals) == aidx and numpy.min(qvals) < 0.05:
                    M.append( gene.SEfc[i])
                    A.append( gene.SEexp[i])
        
        
        plt.plot( A, M, color=colors[aidx], marker='.', linestyle='', label=r'$a = 2^{%s}$' % aVals[aidx] )
    


    
    plt.xlabel( r'$\frac{1}{2}(\log_2 \hat{f}_1 + \log_2 \hat{f}_2)$', size=18)
    plt.ylabel( r'$\log_{2}\hat{f}_1 - \log_2\hat{f}_2$', size=18)
    plt.grid()
    # Shrink current axis by 20%
    box = plot.get_position()
    plot.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    plot.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=True)

    plt.savefig(os.path.join(odir,'mvaSE.%s') % (ext))


def plotMVASE_smear(geneRecords, aVals, odir=os.getcwd(), ext='pdf'):
    """
    Render MA plot
    """
    plt.clf()
    plt.cla()
    plt.autoscale(True)
    colors = ['blue', 'green',  'brown','cyan', 'orange', 'olive', 'pink', 'yellow', 'black',
              'SpringGreen', 'Coral']
    fig = plt.figure()
    plot = fig.add_subplot(111)
    plot.tick_params(axis='both', which='major', labelsize=12)
    plot.tick_params(axis='both', which='minor', labelsize=12)
 
    #plot background surface
    nsigfcs = []
    nsigexps = []

    for gene in geneRecords:
         for i, qvals in enumerate(gene.SEQvals):
            if  gene.SETested[i]and numpy.min(qvals) >= 0.0:
                nsigfcs.append( gene.SEfc[i])
                nsigexps.append( gene.SEexp[i])

    H, xedges, yedges = np.histogram2d(nsigfcs, nsigexps, bins=(30, 30))
    ydiff=0#abs(yedges[0]-yedges[1])
    xdiff=0#abs(xedges[0]-xedges[1])
    
    extent = [min(nsigexps), max(nsigexps), -max(nsigfcs), -min(nsigfcs)]
    
    plt.imshow(H, extent=extent, interpolation='gaussian', origin='upper', cmap=get_cmap('binary'), alpha=1, rasterized=0, vmax=np.max(H)/25.0, aspect=1)
    plt.autoscale(False)
    plt.axhline(0, ls='--', lw=1, color='r')
    for aidx in xrange(len(aVals)):
        M = []
        A = []
        for gene in geneRecords:
            for i, qvals in enumerate(gene.SEQvals):
                if gene.SETested[i] and numpy.argmin(qvals) == aidx and numpy.min(qvals) < 0.05:
                    M.append( gene.SEfc[i])
                    A.append( gene.SEexp[i])
        
        
        plt.plot( A, M, color=colors[aidx], marker='.', linestyle='', label=r'$a = 2^{%s}$' % aVals[aidx] )
    
    
    plt.xlabel( r'$\frac{1}{2}(\log_2 \hat{f}_1 + \log_2 \hat{f}_2)$', size=18)
    plt.ylabel( r'$\log_{2}\hat{f}_1 - \log_2\hat{f}_2$', size=18)
    plt.legend(ncol=2)
    plt.savefig(os.path.join(odir,'mvaSE.%s') % (ext))
    plt.autoscale(True)
    plt.close()


