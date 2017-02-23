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
from iDiffIR.SpliceGrapher.formats.GeneModel import *
from iDiffIR.SpliceGrapher.SpliceGraph       import *
from iDiffIR.SpliceGrapher.shared.GeneModelConverter import *
from multiprocessing import Process, Queue, current_process, freeze_support                                                                           

class IntronModel(object):
    """
    Container for an intron-centric gene model
    """
    def __init__(self, gene, graph, newgene, oGraph, exonsI, retained=None, predicted=False ):
        """
        `gene`    - gene model
        `graph`   - (reduced) splice graph
        `reduced` - whether it was reduced from multi-isoform
        """
        self.gid      = gene.id
        self.ogene    = newgene
        self.intronsR  = sorted(list(set([ (x.minpos+1, x.maxpos) for x in edgeSet(graph)])))
        self.exonsR    = sorted([ (e.minpos, e.maxpos) for e in graph.resolvedNodes()+graph.unresolvedNodes() ])
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
    
    def __cmp__(self, other):
        print other.gid
        return cmp(self.gid, other.gid)

    def __eq__(self, other):
        return self.gid == other.gid

    def __lt__(self, other):
        return self.minpos < other.minpos

    def __le__(self, other):
        return self.minpos <= other.minpos


class ExonModel(IntronModel):
    """Container for exon-centric AS analysis

    Currently only supports exon skipping
    """
    def __init__(self, gene, graph, newgene, oGraph, exonsI, flavorDict, predicted=False, nodeEvents=None, retained=None ):
        IntronModel.__init__(self, gene, graph, newgene, oGraph, exonsI, predicted=predicted, retained=retained)
        self.nodeEvents = nodeEvents
        self.flavorDict = flavorDict

def getGraphs( dirList, verbose ):
    """Load predictions from SpliceGrapher

    Recursively walks through given directory and 
    fetches **SpliceGrapher** predictions

    Parameters
    ----------
    dirList : list
              List of directories to process
    verbose : bool
              True if progress should be printed

    Returns
    -------
    spliceGraphs : dict
                   Mapping of gene ID -> splice graphs for all
                   graphs found

    """
    if verbose:
        sys.stderr.write('loading splice graph predictions...\n')
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
        parents = sorted(list(set([( e.minpos-offset, e.maxpos-offset) for e in node.parents])), 
                         key=lambda x: abs(x[0]-x[1]), reverse=True)
        children = sorted(list(set([( e.minpos-offset, e.maxpos-offset) for e in node.children])), 
                          key=lambda x: abs(x[0]-x[1]), reverse=True)
    else:
        parents = sorted(list(set([( e.minpos-offset, e.maxpos-offset) for e in node.children])), 
                         key=lambda x: abs(x[0]-x[1]), reverse=True)
        children = sorted(list(set([( e.minpos-offset, e.maxpos-offset) for e in node.parents])), 
                          key=lambda x: abs(x[0]-x[1]), reverse=True)
        
    return (parents[0],(node.minpos-offset,node.maxpos-offset), children[0]) 

def decorateNodes( reducedGraph, graph ):
    """Decorate nodes with AS in reduced graph

    .. todo:: add Alt 5, Alt 3 splicing

    """
    flavorDict = {}
    #make Skipped Exons
    SEs = [ getSELocs( node, graph.minpos ) for node in graph.resolvedNodes()+graph.unresolvedNodes() \
                if 'SE' in node.altForms() and node.maxpos-node.minpos > 4]
    if SEs:
        SEs.sort( key=lambda x: x[1][0] )
    flavorDict[ 'SE' ] = SEs
    return flavorDict

def procCluster( geneCluster, graphs, exonic, onlyGraphs, clusterFileHandle ):
    """Process cluster of overlapping genes.

    Process genes whose genomic intervals overlap.  This
    function can be called in parallel.

    See Also
    --------
    geneClusters : makes gene clusters

    Parameters
    ----------
    geneCluster : list
                  List of overlapping genes
    graphs : dict
             Mapping of gene.id -> splice graph prediction.  Can be
             empty if no graphs exist for any of the genes in the cluster.
    exonic : bool
             True if analyzing exon skipping
    onlyGraphs : bool
                 True if only using graphs for generating reduced models.  I.e. not using gene annotations.
    Returns
    -------
    models : list 
             List of reduced gene models

    """
    models = [ ]
    for gene in geneCluster:
        graph = makeSpliceGraph(gene) # splice graph of gene model annotation
        graph.annotate()
        name=graph.name.upper()
        predicted=False # flag if there exists a splice grapher annotation
        # augment gene model graph with predicted graph, if exists
        if name in graphs:
            if graph != graphs[name]:
                predicted = True
                if onlyGraphs:
                    graph = graphs[name]
                else:
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
        overlapIntervals = set()
        clusterFileHandle.write('%s' % ( gene.id ) )
        for otherGene in geneCluster:
            # eliminates same gene collision and sense genes
            # that are on the same strand, start, and end near
            # eachother.  
            if len(geneCluster) == 1 or abs(otherGene.minpos - model.minpos) < 50 and \
               abs(otherGene.maxpos - model.maxpos) < 50 and \
               otherGene.strand==model.strand: continue
            if gene.id == otherGene.id: continue
            overlap = (max(model.minpos, otherGene.minpos), 
                       min( model.maxpos, otherGene.maxpos))
            if overlap[1] > overlap[0]:
                overlapIntervals.add(overlap)
                clusterFileHandle.write('\t%s:(%d,%d)' % (otherGene.id, 
                                                          overlap[0], 
                                                          overlap[1] ) )
        clusterFileHandle.write('\n')
        model.overlap = overlapIntervals
        models.append(model)
    return models


def procCluster_parallel( tasks, output_queue):
    """Process cluster of overlapping genes.

    Process genes whose genomic intervals overlap.  This
    function can be called in parallel.

    See Also
    --------
    geneClusters : makes gene clusters

    Parameters
    ----------
    geneCluster : list
                  List of overlapping genes
    graphs : dict
             Mapping of gene.id -> splice graph prediction.  Can be
             empty if no graphs exist for any of the genes in the cluster.
    exonic : bool
             True if analyzing exon skipping

    Returns
    -------
    models : list 
             List of reduced gene models

    """
    for geneCluster, graphs, exonic in iter(tasks.get, 'STOP'):
        models = [ ]
        for gene in geneCluster:
            graph = makeSpliceGraph(gene) # splice graph of gene model annotation
            graph.annotate()
            name=graph.name.upper()
            predicted=False # flag if there exists a splice grapher annotation
            # augment gene model graph with predicted graph, if exists
            if name in graphs:
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

            # indicate which regions overlap on another gene
            overlapIntervals = set()
            for otherGene in geneCluster:
                # eliminates same gene collision and sense genes
                # that are on the same strand, start, and end near
                # eachother.  
                if len(geneCluster) == 1 or abs(otherGene.minpos - model.minpos) < 50 and \
                   abs(otherGene.maxpos - model.maxpos) < 50 and \
                   otherGene.strand==model.strand: continue
                assert gene.id != otherGene.id # previous filtering should prevent this
                overlap = (max(model.minpos, otherGene.minpos), 
                           min( model.maxpos, otherGene.maxpos))
                if overlap[1] > overlap[0]:
                    overlapIntervals.add(overlap)

            model.overlap = overlapIntervals
            models.append(model)

        output_queue.put( models )

def geneClusters( geneModel, graphs, exonic, onlyGraphs):
    """Cluster genes that overlap each other

    Iterator to cluster and package genes that overlap each other in a 
    chromosome.  :math:`O( n \log n )` time: sorting genes plus linear scan of
    genes.

    Parameters
    ----------
    geneModel : SpliceGrapher.formats.GeneModel.GeneModel
                Gene model for given species
    graphs : dict
             Mapping of geneID -> **SpliceGrapher** prediction 
    exonic : bool
             True if performing SE analysis
    clusterFileHandle : file
                        Open file to write overlapping gene clusters
    Returns
    -------
    geneCluster : list
                  List of genes that overlap eachother
    graphCluster : dict
                 **SpliceGrapher** predictions for genes in cluster.  May be
                 empty if no genes the in cluster have a prediction.
    exonic : bool
             True if performing SE analysis
    """

    # process each chromosome separately
    for chrom in geneModel.getChromosomes():
        genes = sorted(geneModel.getGeneRecords(chrom), 
                       key=lambda gene: gene.minpos)
        # initialize first cluster
        geneCluster = [ genes[0] ]
        minClust = genes[0].minpos
        maxClust = genes[0].maxpos
        for gene in genes[1:]:
            # gene is beyond cluster, yield current cluster
            if gene.minpos >= maxClust:
                clusterGraphs = { }
                for geneIC in geneCluster:
                    if geneIC.id in graphs:
                        clusterGraphs[geneIC.id] = graphs[geneIC.id]
                yield geneCluster, clusterGraphs, exonic, onlyGraphs
                # intitalize the next cluster
                del clusterGraphs
                geneCluster = [gene]
                minClust = gene.minpos
                maxClust = gene.maxpos
            # otherwise it must overlap
            else:
                assert gene.minpos >= minClust 
                geneCluster.append(gene)
                maxClust = max(gene.maxpos, maxClust)

        # yield the last cluster
        clusterGraphs = { }
        for geneIC in geneCluster:
            if geneIC.id in graphs:
                clusterGraphs[geneIC.id] = graphs[geneIC.id]
        yield geneCluster, clusterGraphs, exonic, onlyGraphs

def makeModels( geneModel, outdir, verbose=False, graphDirs=None, graphDirsOnly=None, exonic=False, procs=1 ):
    """Make reduced models for all genes in the geneModel.  

    Genes are processed in parallel using a thread pool 
    with the given number of processors or in a loop if number of processors 
    is 1.

    See Also
    --------
    IntronModel : Reduced gene model class for IR analysis
    ExonModel : Reduced gene model class for SE analysis

    Parameters
    ----------
    geneModel : SpliceGrapher.formats.GeneModel.GeneModel
                Gene model for given species
    verbose : bool
              Flag to write verbose output to sys.stdout
    graphDirs : list
                List of **SpliceGrapher** predictions to augment
                gene models
    exonic : bool
             True if predicting differential exon skipping events
    procs : int
            Number of processors to use

    Returns
    -------
    models : list
             List of reduced models

    """
    # load splice grapher predictions if given
    graphs = [ ]
    onlyGraphs = False
    if graphDirs:
        graphs = getGraphs( graphDirs, verbose )
    if graphDirsOnly:
        graphs = getGraphs( graphDirsOnly, verbose )
        onlyGraphs = True
    models   = [ ]
    if verbose:
        sys.stderr.write('Building splicing models\n')
    indicator = ProgressIndicator(10000, verbose=verbose)
    if outdir:
        clusterFileHandle = open(os.path.join(outdir, 'lists', 'geneClusters.txt'), 'w')
    else:
        clusterFileHandle = open(os.devnull,"w")

    def updateIndicator(pModels):
        """Update indicator

        Parameters
        ----------
        pModels : list
                  List of reduced models

        """
        for i in xrange(len(pModels)):             
            indicator.update()

    # run parallel
    if False:
        # create processor pool for parallel calls
        task_queue = Queue()
        status_queue = Queue()
        nTasks = 0
        for clusterTuple in geneClusters(geneModel, graphs, exonic, onlyGraphs):
            task_queue.put( clusterTuple )
            nTasks += 1

        for _ in xrange(procs):
            Process(target=procCluster,
                    args=(task_queue, status_queue)).start()

        for _ in xrange(procs):
            task_queue.put('STOP')

        for _ in xrange(nTasks):
            processedModels = status_queue.get()
            models.extend(processedModels)
            updateIndicator(processedModels)

    # run serial
    else:
        for clusterTuple in geneClusters(geneModel, graphs, 
                                         exonic, onlyGraphs):
            processedModels = procCluster(*clusterTuple, 
                                          clusterFileHandle=clusterFileHandle)
            models.extend(processedModels)
            updateIndicator(processedModels)
    if verbose:
        sys.stderr.write('%d genes' % indicator.count())
    indicator.finish()
    clusterFileHandle.close()
    return models
    
def makeReducedExonModel(gene, graph):
    """
    Finds exonic regions that are in common across all isoforms
    """

    minGene = graph.minpos
    maxGene = graph.maxpos
    #set min gene area
    nodesR = set([ (x.minpos, x.maxpos) for x in graph.resolvedNodes()+graph.unresolvedNodes()])
    exons = [ ]
    edges = set([ (x.minpos, x.maxpos) for x in edgeSet(graph)])

    roots = set([ (x.minpos, x.maxpos) for x in graph.resolvedNodes()+graph.unresolvedNodes() if x.isRoot()])

    leaves = set([ (x.minpos, x.maxpos) for x in graph.resolvedNodes()+graph.unresolvedNodes() if x.isLeaf()])
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
        exons = set([ (x.minpos, x.maxpos) for x in graph.resolvedNodes()+graph.unresolvedNodes()])
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
    nodes   = set([ (x.minpos, x.maxpos) for x in graph.resolvedNodes()+graph.unresolvedNodes()])
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

