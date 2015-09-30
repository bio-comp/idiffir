#! /usr/bin/env python
"""Program for generating known splice junction sequences."""
from SpliceGrapher.shared.utils        import *
from SpliceGrapher.shared.GeneModelConverter import *
from SpliceGrapher.predict.SpliceSite  import *
from SpliceGrapher.formats.GeneModel   import *
from SpliceGrapher.SpliceGraph         import *
from SpliceGrapher.formats.loader import loadGeneModels

from glob     import glob
from optparse import OptionParser
import os, sys, warnings

parser = OptionParser(usage='Generate GFF file for MISO of exon skipping events')
parser.add_option('-m', dest='model',    default=SG_GENE_MODEL, help='Gene model GFF file [default: %default]')
parser.add_option('-o', dest='outfile',  default=None,          help='Output file [default: %default]')
parser.add_option('-v', dest='verbose',  default=False,         help='Verbose mode [default: %default]', action='store_true')
parser.add_option( '-s', dest='graphPaths', default=None, help='File containing paths to splice graphs to augment gene models')
opts, args = parser.parse_args(sys.argv[1:])

#-------------------------------------------------------
# Main program
#-------------------------------------------------------
opts, args = parser.parse_args(sys.argv[1:])

errStrings = []
if not opts.model : errStrings.append('** No GFF gene model specified.  Set SPLICEGRAPHER_GENE_MODEL or use the -m option.')
if errStrings :
    parser.print_help()
    sys.stderr.write('\n%s\n' % '\n'.join(errStrings))
    sys.exit(1)

geneModel = loadGeneModels(opts.model, verbose=opts.verbose, alltypes=True)
genes     = geneModel.getAllGenes(geneFilter=gene_type_filter)
genes.sort()

outStream = open(opts.outfile, 'w') if opts.outfile else sys.stdout

totalRecs  = 0
uniqueEvents = 0
exceptions = set([])
chrCounter = {}
indicator = ProgressIndicator(10000, description=' genes', verbose=opts.verbose)

if opts.graphPaths != None:
    if opts.verbose: sys.stderr.write('Loading aux. splice graphs\n')
    indicator = ProgressIndicator(10000, description=' files', verbose=opts.verbose)
    with open( opts.graphPaths, 'r' ) as fin:
        for line in fin:
            fname = line.strip()
            if not os.path.exists(fname):
                if opts.verbose: 
                    sys.stderr.write('Missing splice graph file: %s\n' % fname)
                    continue
            indicator.update()
            graph = getFirstGraph( fname )
            geneName = graph.name.upper()
            if geneName in spliceGraphs:
                spliceGraphs[geneName] = spliceGraphs[geneName].union( graph )
            else:
                spliceGraphs[geneName] = graph
            #for node in graph.unresolvedNodes():
            #    spliceGraphs[geneName].addNode( node.id, node.start, node.end)
        
    indicator.finish()
    if opts.verbose : sys.stderr.write('Loaded %d aux. splice graphs\n' % (indicator.ctr))

def getEventLocs( node ):
    if node.strand == '+':
        parents = list(set([( e.minpos, e.maxpos) for e in node.parents]))
        children = list(set([( e.minpos, e.maxpos) for e in node.children]))
    else:
        parents = list(set([( e.minpos, e.maxpos) for e in node.children]))
        children = list(set([( e.minpos, e.maxpos) for e in node.parents]))
        
    return [ (parents[i],(node.minpos,node.maxpos), children[j])\
                 for i in xrange(len(parents)) for j in xrange(len(children))]
for gene in genes :
    if gene.chromosome not in chrCounter :
        chrCounter[gene.chromosome] = 0
    #if opts.verbose : indicator.update()


    # Create splice graph for gene:
    try :
        geneGraph = makeSpliceGraph(gene)#geneModelToSpliceGraph(g)
    except ValueError, ve:
        sys.stderr.write("Unable to create graph for %s\n" % gene.name)
        continue 

    if geneGraph.name.upper() in spliceGraphs:
        geneGraph = geneGraph.union(spliceGraphs[ geneGraph.name.upper() ])
    geneGraph.annotate()


    if not geneGraph.hasAS(): continue
    if 'SE' not in geneGraph.altForms(): continue

    seNodes = [node for node in geneGraph.resolvedNodes() if node.isSkippedExon()]

            
    for nNum,node in enumerate(seNodes):
        uniqueEvents += 1
        nInt = 0
        for exon5, exonS, exon3 in getEventLocs(node):
            exon5, exonS, exon3
            totalRecs += 1
            nInt+=1
            eventID = '%s:%d-%d:%d:%d' % (gene.id, exonS[0], exonS[1], nNum+1, nInt)            
            outStream.write('%s\tIR\tgene\t%d\t%d\t.\t%s\t.\tID=%s;Name=%s\n' % (gene.chromosome, 
                                                                                 exon5[0], 
                                                                                 exon3[1], 
                                                                                 node.strand, 
                                                                                 eventID, eventID))

            irID = '%s:%d-%d:%d:%d' % (gene.id, exonS[0], exonS[1], nNum+1, nInt)
            #make mRNA record for SE
            outStream.write('%s\tSE\tmRNA\t%d\t%d\t.\t%s\t.\tID=%s:SE;Parent=%s\n' % \
                                (gene.chromosome, exon5[0], exon3[1], node.strand, 
                                 irID, eventID))
            #make mRNA record for CS
            outStream.write('%s\tSE\tmRNA\t%d\t%d\t.\t%s\t.\tID=%s:CS;Parent=%s\n' % \
                                (gene.chromosome,  exon5[0], exon3[1], node.strand, 
                                 irID, eventID))                
            #make 5' exon record for SE
            outStream.write('%s\tSE\texon\t%d\t%d\t.\t%s\t.\tID=%s:se5Exon;Parent=%s:SE\n' % \
                                (gene.chromosome, exon5[0], exon5[1], node.strand, 
                                 irID, irID))
            #make 3' exon record for SE
            outStream.write('%s\tSE\texon\t%d\t%d\t.\t%s\t.\tID=%s:se3Exon;Parent=%s:SE\n' % \
                                (gene.chromosome, exon3[0], exon3[1], node.strand, 
                                 irID, irID))

            #make 5' exon record for CS
            outStream.write('%s\tSE\texon\t%d\t%d\t.\t%s\t.\tID=%s:se5Exon;Parent=%s:CS\n' % \
                                (gene.chromosome, exon5[0], exon5[1], node.strand, 
                                 irID, irID))
            #make 3' exon record for CS
            outStream.write('%s\tSE\texon\t%d\t%d\t.\t%s\t.\tID=%s:se3Exon;Parent=%s:CS\n' % \
                                (gene.chromosome, exon3[0], exon3[1], node.strand, 
                                 irID, irID))
            #make cassette exon record for CS
            outStream.write('%s\tSE\texon\t%d\t%d\t.\t%s\t.\tID=%s:cassExon;Parent=%s:CS\n' % \
                                (gene.chromosome, exonS[0], exonS[1], node.strand, 
                                 irID, irID))





        
outStream.flush()
outStream.close()
if opts.verbose : indicator.finish()

sys.stderr.write("Wrote %d total SE events\n" % totalRecs)
sys.stderr.write("Wrote %d unique SE events\n" % uniqueEvents)



