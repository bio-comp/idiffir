#!/bin/python
"""Program for generating miso IR examples."""
from iDiffIR.SpliceGrapher.shared.utils        import *
from iDiffIR.SpliceGrapher.shared.GeneModelConverter import *
from iDiffIR.SpliceGrapher.predict.SpliceSite  import* 
from iDiffIR.SpliceGrapher.formats.GeneModel   import *
from iDiffIR.SpliceGrapher.SpliceGraph         import *
from iDiffIR.SpliceGrapher.formats.loader import loadGeneModels

from glob     import glob
from optparse import OptionParser
import os
import sys


def writeEvent(event, outstream, eventID, a5_map_stream, a3_map_stream):
    name = '%s:%d:%d:%s:%s:%s:%d:%s' % (event.geneID, event.minpos,
                                           event.maxpos, event.pathType,
                                           event.chrom, 
                                           event.strand, eventID, event.label)
    if event.AStype == 'A5':
        node1 = event.upPath[0]
        node2 = event.downPath[0]
        positions = sorted((node1.donorEnd(), node2.donorEnd()))
        a5_map_stream.write('%s\t%d,%d\n' % (name, positions[0], positions[1]))
    else:
        node1 = event.upPath[-1]
        node2 = event.downPath[-1]
        positions = sorted((node1.acceptorEnd(), node2.acceptorEnd()))
        a3_map_stream.write('%s\t%d,%d\n' % (name, positions[0], positions[1]))


    # write gene record
    outstream.write( '%s\t%s\tgene\t%d\t%d\t.\t%s\t.\tID=%s;Name=%s\n' % (
        event.chrom, event.AStype, event.minpos, event.maxpos,
        event.strand, name, name ) )

    # write upstream mRNA record
    minpos = min( [ node.minpos for node in event.upPath ] )
    maxpos = max( [ node.maxpos for node in event.upPath ] )
    outstream.write('%s\t%s\tmRNA\t%d\t%d\t.\t%s\t.\tID=%s;Parent=%s\n' % (
        event.chrom, event.AStype, minpos, maxpos,
        event.strand, '%s:upstream' %name, name))

    # write downstream mRNA record
    minpos = min( [ node.minpos for node in event.downPath ] )
    maxpos = max( [ node.maxpos for node in event.downPath ] )
    outstream.write('%s\t%s\tmRNA\t%d\t%d\t.\t%s\t.\tID=%s;Parent=%s\n' % (
        event.chrom, event.AStype, minpos, maxpos,
        event.strand, '%s:downstream' %name, name))

    # write upstream exon records
    for i,node in enumerate(event.upPath):
        outstream.write('%s\t%s\texon\t%d\t%d\t.\t%s\t.\tID=%s;Parent=%s\n' % (
            event.chrom, event.AStype, node.minpos, node.maxpos,
            event.strand, '%s:upstream:exon_%d' % (name, i+1), '%s:upstream' %name))

    # write upstream exon records
    for i,node in enumerate(event.downPath):
        outstream.write('%s\t%s\texon\t%d\t%d\t.\t%s\t.\tID=%s;Parent=%s\n' % (
            event.chrom, event.AStype, node.minpos, node.maxpos,
            event.strand, '%s:downstream:exon_%d' % (name, i+1), '%s:downstream' %name))


def upstream_a5( node1, node2 ):
    """
    Determines if node1's 5' splice site is upstream of node2's
    """
    if node1.strand == '+':
        return node1.donorEnd() < node2.donorEnd()
    else:
        return node1.donorEnd() > node2.donorEnd()


def upstream_a3( node1, node2 ):
    """
    Determines if node1's 3' splice site is upstream of node2's
    """
    if node1.strand == '+':
        return node1.acceptorEnd() < node2.acceptorEnd()
    else:
        return node1.acceptorEnd() > node2.acceptorEnd()

def overlap( node1, node2 ):
    """
    Determines if nodes overlap
    """
    return min( node1.maxpos, node2.maxpos) - max( node1.minpos, node2.minpos) > 0

def resolve_a5( node1, node2 ):
    """
    Return alternative 5p paths representing
    minimal isoform distinction of the event.
    """
    # ignore roots and leaves of splice graphs
    if node1.isRoot() or node1.isLeaf() or node2.isRoot() or node2.isLeaf(): 
        return (None, None, None)
    
    # alternative donors share a child ( cake work )
    for child1 in node1.children:
        for child2 in node2.children:
            if child1 == child2:
                return ( (node1, child1), (node2, child2), 'S' )

    # alternative donors share a child 
    for child1 in node1.children:
        for child2 in node2.children:
            if overlap( child1, child2 ):
                return ( (node1, child1), (node2, child2), 'O' )
    
    # check grandchildren of node1
    for child1 in node1.children:
        if child1.isLeaf(): continue
        for grandchild1 in child1.children:
            for child2 in node2.children:
                if grandchild1 == child2:
                    return ( (node1, child1, grandchild1), (node2, child2), 'S3' )

    # check grandchildren of node1
    # alternative donors overlap grand child - child
    for child1 in node1.children:
        if child1.isLeaf(): continue
        for grandchild1 in child1.children:
            for child2 in node2.children:
                if overlap( grandchild1, child2 ):
                    return ( (node1, child1, grandchild1), (node2, child2), 'O3' )

    # check grandchildren of node2
    # alternative donors same grand child - child
    for child1 in node1.children:
        for child2 in node2.children:
            if child2.isLeaf(): continue
            for grandchild2 in child2.children:
                if child1 == grandchild2:
                    return ( (node1, child1), (node2, child2, grandchild2), 'S3' )

    # check grandchildren of node2
    # alternative donors overlap grand child - child
    for child1 in node1.children:
        for child2 in node2.children:
            if child2.isLeaf(): continue
            for grandchild2 in child2.children:
                if overlap( child1, grandchild2 ):
                    return ( (node1, child1), (node2, child2, grandchild2), 'O3' )

    return (None, None, None)

def resolve_a3( node1, node2 ):
    """
    Return alternative 5p paths representing
    minimal isoform distinction of the event.
    """
    # ignore roots and leaves of splice graphs
    if node1.isRoot() or node1.isLeaf() or node2.isRoot() or node2.isLeaf(): 
        return (None, None, None)
    
    # alternative donors share a parent ( cake work )
    for parent1 in node1.parents:
        for parent2 in node2.parents:
            if parent1 == parent2:
                return ( (parent1, node1), (parent2, node2), 'S' )

    # alternative donors share a parent
    for parent1 in node1.parents:
        for parent2 in node2.parents:
            if overlap(parent1, parent2):
                return ( (parent1, node1), (parent2, node2), 'O' )
    
    # check grandparents of node1
    for parent1 in node1.parents:
        if parent1.isRoot(): continue
        for grandparent1 in parent1.parents:
            for parent2 in node2.parents:
                if grandparent1 == parent2:
                    return ( (grandparent1, parent1, node1), (parent2, node2), 'S3' )

    # check grandchildren of node1
    # alternative donors overlap grand child - child
    for parent1 in node1.parents:
        if parent1.isRoot(): continue
        for grandparent1 in parent1.parents:
            for parent2 in node2.parents:
                if overlap(grandparent1, parent2):
                    return ( (grandparent1, parent1, node1), (parent2, node2), 'O3' )



    # check grandchildren of node2
    # alternative donors same grand child - child
    for parent1 in node1.parents:
        for parent2 in node2.parents:
            if parent2.isRoot(): continue
            for grandparent2 in parent2.parents:
                if parent1 == grandparent2:
                    return ( (parent1, node1), (grandparent2, parent2, node2), 'S3' )

    # check grandchildren of node2
    # alternative donors overlap grand child - child
    for parent1 in node1.parents:
        for parent2 in node2.parents:
            if parent2.isRoot(): continue
            for grandparent2 in parent2.parents:
                if overlap(parent1,grandparent2):
                    return ( (parent1, node1), (grandparent2, parent2, node2), 'O3' )


    return (None, None, None)

    

class AltEvent(object):
    """
    Container for alternative splicing event
    """
    def __init__( self, path1, path2, pathType, AStype, chrom, geneID, label ):
        self.upPath = path1
        self.downPath = path2
        self.pathType = pathType
        self.AStype = AStype
        self.chrom = chrom
        self.strand = path1[0].strand
        self.geneID = geneID
        minpos = min( [node.minpos for node in path1] )
        self.minpos = min( minpos, min( [node.minpos for node in path2] ) )
        maxpos = max( [node.maxpos for node in path1] )
        self.maxpos = max( maxpos, max( [node.maxpos for node in path2] ) )
        self.label = label

def build_parser():
    parser = OptionParser()
    parser.add_option('-m', dest='model',    default=SG_GENE_MODEL, help='Gene model GFF file [default: %default]')
    parser.add_option('-o', dest='outbase',  default='AS',          help='Output file base')
    parser.add_option('-v', dest='verbose',  default=False,         help='Verbose mode [default: %default]', action='store_true')
    parser.add_option('-s', dest='graphPaths', default=None, help='File containing paths to splice graphs to augment gene models')
    return parser


def parse_args(argv=None):
    parser = build_parser()
    opts, args = parser.parse_args(argv)
    errStrings = []
    if not opts.model:
        errStrings.append('** No GFF gene model specified.  Set SPLICEGRAPHER_GENE_MODEL or use the -m option.')
    if errStrings:
        parser.print_help()
        sys.stderr.write('\n%s\n' % '\n'.join(errStrings))
        raise SystemExit(1)
    return opts, args


def load_splice_graphs(graph_paths, verbose):
    spliceGraphs = {}
    if graph_paths is not None:
        if verbose:
            sys.stderr.write('Loading aux. splice graphs\n')
        indicator = ProgressIndicator(10000, description=' files', verbose=verbose)
        with open(graph_paths, 'r') as fin:
            for line in fin:
                fname = line.strip()
                if not os.path.exists(fname):
                    if verbose:
                        sys.stderr.write('Missing splice graph file: %s\n' % fname)
                    continue
                indicator.update()
                graph = getFirstGraph(fname)
                geneName = graph.name.upper()
                if geneName in spliceGraphs:
                    spliceGraphs[geneName] = spliceGraphs[geneName].union(graph)
                else:
                    spliceGraphs[geneName] = graph
        indicator.finish()
        if verbose:
            sys.stderr.write('Loaded %d aux. splice graphs\n' % indicator.ctr)
    return spliceGraphs


def main(argv=None):
    opts, _ = parse_args(argv)
    geneModel = loadGeneModels(opts.model, verbose=opts.verbose, alltypes=True)
    spliceGraphs = load_splice_graphs(opts.graphPaths, opts.verbose)

    with open(opts.outbase + '_a5.gff', 'w') as outStream_a5, \
         open(opts.outbase + '_a3.gff', 'w') as outStream_a3, \
         open('a5Maps.txt', 'w') as a5f, \
         open('a3Maps.txt', 'w') as a3f:
        for chrm in geneModel.getChromosomes():
            if opts.verbose:
                sys.stderr.write('Processing genes from chromosome: %s\n' % chrm)
            indicatorG = ProgressIndicator(10000, description=' genes', verbose=opts.verbose)
            genes = geneModel.getGeneRecords(chrm, geneFilter=gene_type_filter)
            genes.sort()

            for g in genes:
                if opts.verbose:
                    indicatorG.update()
                geneGraph = makeSpliceGraph(g)
                if geneGraph.name.upper() in spliceGraphs:
                    geneGraph = geneGraph.union(spliceGraphs[geneGraph.name.upper()])
                geneGraph.annotate()

                a5 = [node for node in geneGraph.resolvedNodes() if node.isAltDonor()]
                a3 = [node for node in geneGraph.resolvedNodes() if node.isAltAcceptor()]
                a5Events = []
                a3Events = []
                for node1 in a5:
                    for node2 in a5:
                        if upstream_a5(node1, node2) and overlap(node1, node2):
                            path1, path2, pathType = resolve_a5(node1, node2)
                            if pathType is not None:
                                label = 'P' if node1.isPredicted() or node2.isPredicted() else 'K'
                                a5Events.append(AltEvent(path1, path2, pathType, 'A5', chrm, geneGraph.name, label))

                for i, event in enumerate(a5Events):
                    writeEvent(event, outStream_a5, i + 1, a5f, a3f)

                for node1 in a3:
                    for node2 in a3:
                        if upstream_a3(node1, node2) and overlap(node1, node2):
                            path1, path2, pathType = resolve_a3(node1, node2)
                            if pathType is not None:
                                label = 'P' if node1.isPredicted() or node2.isPredicted() else 'K'
                                a3Events.append(AltEvent(path1, path2, pathType, 'A3', chrm, geneGraph.name, label))

                for i, event in enumerate(a3Events):
                    writeEvent(event, outStream_a3, i + 1, a5f, a3f)

            if opts.verbose:
                indicatorG.finish()
                sys.stderr.write('%d genes\n' % indicatorG.ctr)

    return 0


if __name__ == '__main__':
    raise SystemExit(main())
