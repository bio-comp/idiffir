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
#
import os, sys, numpy
numpy.seterr(invalid='raise')
from SpliceGrapher.formats.GeneModel import *
from SpliceGrapher.SpliceGraph       import *
from SpliceGrapher.shared.GeneModelConverter import *

class IntronModel(object):
    """
    Container for intron-centric gene model and expression
    """
    def __init__(self, gene, graph, newgene, oGraph, exonsI, retained=None,predicted=False ):
        """
        `gene`    - gene model
        `graph`   - (reduced) splice graph
        `reduced` - whether it was reduced from multi-isoform
        """
        self.gid      = gene.id
        self.ogene    = newgene
        self.intronsR  = sorted(list(set([ (x.minpos+1, x.maxpos) for x in edgeSet(graph)])))
        self.exonsR    = sorted([ (e.minpos, e.maxpos) for e in graph.resolvedNodes() ])
        self.graph    = graph
        self.gene     = gene
        self.retained = retained if retained else [False]* len(self.intronsR)
        self.predicted = predicted
        self.strand   = gene.strand
        self.chrom    = gene.chromosome
        self.minpos   = graph.minpos
        self.maxpos   = graph.maxpos
        self.overlap = None
        # convert absolute intron and exon boundaries into relative
        self.origGraph = oGraph
        self.introns = [(s-self.minpos, e-self.minpos) for s,e in self.intronsR]
        self.exons = [(s-self.minpos, e-self.minpos) for s,e in self.exonsR]
        self.exonsI = [(s-self.minpos, e-self.minpos) for s,e in exonsI]
    
        
class ExonModel(IntronModel):
    """
    Container for exon-centric AS analysis

    Currently only supports exon skipping
    """
    def __init__(self, gene, graph, newgene, oGraph, exonsI, flavorDict, predicted=False, nodeEvents=None, retained=None ):
        IntronModel.__init__(self, gene, graph, newgene, oGraph, exonsI, predicted=predicted, retained=retained)
        self.nodeEvents = nodeEvents
        self.flavorDict = flavorDict

def getGraphs( dirList, verbose ):
    """
    Load predictions from SpliceGrapher

    Recursively walks through given directory
    """
    if verbose:
        sys.stderr.write('\tloading splice graph predictions...\n')
    indicator = ProgressIndicator(10000, verbose=verbose)
    spliceGraphs = {}
    # Let's take a walk :)
    for gDir in dirList:
        for root, dirs, files in os.walk( gDir ):
            for fname in files:
                fullPath = os.path.join(root,fname)
                indicator.update()
                try:
                    newgraph = getFirstGraph(fullPath)
                    newgraph.annotate()
                    name = newgraph.name.upper()
                    if name in spliceGraphs:
                        spliceGraphs[name] = spliceGraphs[name].union(newgraph)
                    else:
                        spliceGraphs[name] = newgraph
                except Exception :
                    continue
    indicator.finish()
    return spliceGraphs 


class DecoratedNode(object):
    """
    Containter for a decorated node
    """
    def __init__(self, minpos, maxpos, AS_flavor):
        self.minpos = minpos
        self.maxpos = maxpos
        self.AS_flavor = AS_flavor

def getSELocs( node, offset ):
    """
    Return a 3-tuple of 2-tuples representing (upstream, skipped, and downstream) exons

    """
    if node.strand == '+':
        parents = sorted(list(set([( e.minpos-offset, e.maxpos-offset) for e in node.parents])), key=lambda x: abs(x[0]-x[1]), reverse=True)
        children = sorted(list(set([( e.minpos-offset, e.maxpos-offset) for e in node.children])), key=lambda x: abs(x[0]-x[1]), reverse=True)
    else:
        parents = sorted(list(set([( e.minpos-offset, e.maxpos-offset) for e in node.children])), key=lambda x: abs(x[0]-x[1]), reverse=True)
        children = sorted(list(set([( e.minpos-offset, e.maxpos-offset) for e in node.parents])), key=lambda x: abs(x[0]-x[1]), reverse=True)
        
    return (parents[0],(node.minpos-offset,node.maxpos-offset), children[0]) 

def decorateNodes( reducedGraph, graph ):
    """
    Decorate nodes with as in reduced graph

    """
    flavorDict = {}
    #make Skipped Exons
    SEs = [ getSELocs( node, graph.minpos ) for node in graph.resolvedNodes() \
                if 'SE' in node.altForms() and node.maxpos-node.minpos > 4]
    if SEs:
        SEs.sort( key=lambda x: x[1][0] )
    flavorDict[ 'SE' ] = SEs
    return flavorDict

def makeModels( geneModel, verbose=True, graphDirs=None, exonic=False ):
    """
    Wrapper to make reduced models for 
    all genes in the geneModel.
    """
    # load splice grapher predictions if given
    if graphDirs:
        graphs = getGraphs( graphDirs, verbose )
    models   = [ ]
    indicator = ProgressIndicator(10000, verbose=verbose)
    for gene in geneModel.getAllGenes(geneFilter=gene_type_filter):
        graph = makeSpliceGraph(gene) # splice graph of gene model annotation
        graph.annotate()
        name=graph.name.upper()
        predicted=False # flag if there exists a splice grapher annotation
        # augment gene model graph with predicted graph, if exists
        if graphDirs and name in graphs:
            if graph != graphs[name]:
                predicted = True
                graph=graph.union(graphs[name])
        reducedGraph, irs, newgene = makeReducedGraph(gene, graph)
        exonsI = makeReducedExonModel(gene, graph) # graph where exons exhibit no AS
        # do this if we're looking at exon skipping
        if exonic:
            flavorDict = decorateNodes( reducedGraph, graph)
            model = ExonModel( gene, reducedGraph, newgene, graph, exonsI, 
                               flavorDict, retained=irs, predicted=predicted)
        # looking at IR so do this instead
        else:
            model = IntronModel( gene, reducedGraph, newgene, graph, exonsI,
                                 retained=irs, predicted=predicted)


        # indicate which introns overlap on another strand
        overlappers = [ ]
        if not exonic:
            overlappers = [ x for x in geneModel.getGenesInRange(model.chrom, model.minpos, model.maxpos) if x.id != model.gid and x.strand != model.strand]
        if len(overlappers) == 0: model.overlap = [ False ] * len(model.introns)
        else:
            overlap = [ ]
            ranges = [ ]
            ogenes = [ ]
            for x in overlappers:
                if abs(min(x.start(), x.end()) - model.minpos) < 50 and abs(max(x.start(), x.end()) - model.maxpos) < 50 and x.strand==model.strand: continue
                ranges.append( sorted( (x.start(), x.end()) ) )
                ogenes.append( x.id )
            for intronR in model.intronsR:
                lim = int((abs(intronR[1] - intronR[0]) + 1) * 0.25)
                for i, o_range  in enumerate(ranges):
                    minOverlap, maxOverlap = o_range
                    if minOverlap >= intronR[0] and minOverlap <= intronR[1] and ( minOverlap > intronR[0] + lim or minOverlap < intronR[1] - lim) or \
                       maxOverlap >= intronR[0] and maxOverlap <= intronR[1] and ( maxOverlap > intronR[0] + lim or maxOverlap < intronR[1] - lim) or \
                       minOverlap <= intronR[0] and maxOverlap >= intronR[1]:
                        overlap.append(True)
                        break
                else:
                    overlap.append(False)
            model.overlap = overlap
        models.append(model)
        indicator.update()
    indicator.finish()
    return models
    
def makeReducedExonModel(gene, graph):
    """
    Finds exonic regions that are in common across all isoforms
    """

    minGene = graph.minpos
    maxGene = graph.maxpos
    
    nodesR = set([ (x.minpos, x.maxpos) for x in graph.resolvedNodes()])
    exons = [ ]
    edges = set([ (x.minpos, x.maxpos) for x in edgeSet(graph)])

    roots = set([ (x.minpos, x.maxpos) for x in graph.resolvedNodes() if x.isRoot()])

    leaves = set([ (x.minpos, x.maxpos) for x in graph.resolvedNodes() if x.isLeaf()])
    leafEnd = min( [x[1] for x in leaves] ) if gene.strand == '+' else \
              max( [x[0] for x in leaves])
    if gene.strand == '+':
        for s,e in roots:
            disconnected = False
            if s > leafEnd:
                disconnected = True
                edges.add( (leafEnd, s) )
            if not disconnected and s > minGene:
                edges.add( (minGene, s) )
        for s2,e2 in leaves:
            if e2 <  maxGene:
                edges.add( (e2, maxGene) )                    

    else:
        for s,e in roots:
            disconnected = False
            if e < leafEnd:
                disconnected = True
                edges.add( ( e, leafEnd) )

            if not disconnected and e < maxGene:
                edges.add( (e, maxGene))
        for s2,e2 in leaves:
            if s2 > minGene:
                edges.add( (minGene, s2) )                    
            
    edges = list(edges)
    edges.sort( key=lambda x: x[1] - x[0] + 1, reverse=True )
    exons =  set([ ])
    # remove skipped exons and nodes outside of true min and max
    while len(nodesR) > 0:
        s,e = nodesR.pop()
        for s2, e2 in edges:
        # check if cassette exon
            if s >= s2 and e <= e2:
                break
            # exon subsumes intron (IR)
            elif s2 > s and e2 < e:
                nodesR.add((s2,s))
                nodesR.add((e2,e))
                break
            # overlap left
            elif s2 <= s and e2 > s:
                nodesR.add( (e2,e) )
                break
            # overlap right
            elif s2 < e and e2 >= e:
                nodesR.add( (s2,e))
                break
        else:
            exons.add( (s,e))


    if len(exons) == 0:
        exons = set([ (x.minpos, x.maxpos) for x in graph.resolvedNodes()])
    exons = sorted(list(exons))

    return exons
                
def makeReducedModel( gene, graph ):
    """
    Extracts the exons of the reduced model and annotates
    the resultant introns if they are retained.
    """
    minGene = graph.minpos
    maxGene = graph.maxpos

    intronsR = set([ (x.minpos, x.maxpos) for x in edgeSet(graph)])
    nodes   = set([ (x.minpos, x.maxpos) for x in graph.resolvedNodes()])
    introns = set([])
    while len( intronsR ) > 0:
        s,e = intronsR.pop()
        for s2,e2 in nodes:
            # overlap left
            if s2 <= s and e2 > s and e2 < e:
                intronsR.add( (e2,e) )
                break
            # overlap right
            elif s2 > s and s2 < e and e2 >= e:
                intronsR.add( (s, s2) )
                break
            # node subsumed by intron
            elif s2 > s and e2 < e:
                intronsR.add( (s,s2) )
                intronsR.add( (e2,e) )
                break
            # lacking splice junction support
            #|||||||||||||||-------------|||||||||||||||||
            #              ||||||||||||||||||||----------- 
            #                     or
            #|||||||||||||||--------------||||||||||||||||
            #------------||||||||||||||||||
            elif abs( s- s2) < 3 or abs( e - e2 ) < 3:
                break
        else:
            introns.add( (s,e) )

    introns  = sorted(list(introns))
    irs = [ ]
    exons    = [ ]
    es = minGene
    for s,e in introns:
        exons.append( (es, s) )
        es = e
        for s2, e2 in nodes:
            if s > s2-2 and e < e2 - 2:
                irs.append( True )
                break
        else:
            irs.append( False )
    if es > minGene:
        exons.append( (es,maxGene) )
    if len(introns) == 0:
        return [ (minGene, maxGene)], []
    return exons, irs


def makeReducedGraph(gene, graph):
    exonRanges, irs = makeReducedModel(gene, graph)
    start = graph.minpos if gene.strand == '+' else graph.maxpos
    end = graph.maxpos if gene.strand == '+' else graph.minpos
    newgene = Gene( gene.id+'_RM', gene.note, start, end, gene.chromosome,
                    gene.strand)
    isoName = newgene.id+'I'
    iso  = Isoform( isoName, exonRanges[-1][1], 
                    exonRanges[0][0], gene.chromosome, gene.strand)
    isoAttr   = {PARENT_FIELD:newgene.id, NAME_FIELD:isoName, ID_FIELD:isoName}
    for minpos, maxpos in exonRanges:
        newgene.addExon( iso, Exon( minpos, maxpos, gene.chromosome, gene.strand))

    newgraph = geneModelToSpliceGraph(newgene)
    return newgraph, irs, newgene

