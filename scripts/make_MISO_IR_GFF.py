#!/usr/bin/env python
"""Program for generating miso IR examples."""
from iDiffIR.SpliceGrapher.shared.utils        import *
from iDiffIR.SpliceGrapher.shared.GeneModelConverter import *
from iDiffIR.SpliceGrapher.predict.SpliceSite  import* 
from iDiffIR.SpliceGrapher.formats.GeneModel   import *
from iDiffIR.SpliceGrapher.SpliceGraph         import *
from iDiffIR.SpliceGrapher.formats.loader import loadGeneModels

from glob     import glob
import argparse
import os, sys, warnings, numpy
from pathlib import Path

def writeIntrons( gene, introns, label, outstream ):
    i = 0
    for intron in introns:
        i += 1
        if intron.strand == '+':
            start = sorted(intron.p5)[0]
            stop = sorted(intron.p3)[1]
        else:
            start = sorted(intron.p3)[0]
            stop = sorted(intron.p5)[1]
            
        name = '%s:%d:%d:%d:%d:%s:%s:%d:%s' % (gene.id, start,
                                               stop, intron.minpos, intron.maxpos,
                                               intron.chrom,
                                               intron.strand, i, label)

        # write gene record
        outstream.write( '%s\tIR\tgene\t%d\t%d\t.\t%s\t.\tID=%s;Name=%s\n' % (
            intron.chrom, start, stop,
            intron.strand, name, name ) )
        
        # write IR mRNA record
        outstream.write('%s\tIR\tmRNA\t%d\t%d\t.\t%s\t.\tID=%s;Parent=%s\n' % (
            intron.chrom, start, stop,
            intron.strand, '%s:IR' %name, name))

        # write IE mRNA record
        outstream.write('%s\tIR\tmRNA\t%d\t%d\t.\t%s\t.\tID=%s;Parent=%s\n' % (
            intron.chrom, start, stop,
            intron.strand, '%s:IE' %name, name))

        # write IR exon record
        outstream.write('%s\tIR\texon\t%d\t%d\t.\t%s\t.\tID=%s;Parent=%s\n' % (
            intron.chrom, start, stop,
            intron.strand, '%s:IR:exon' %name, '%s:IR' %name))

        if intron.strand == '+':
            # write IE 5p exon record
            outstream.write('%s\tIR\texon\t%d\t%d\t.\t%s\t.\tID=%s;Parent=%s\n' % (
                intron.chrom, sorted(intron.p5)[0], sorted(intron.p5)[1],
                intron.strand, '%s:IE:5pexon' %name, '%s:IE' %name))

            # write IE 3p exon record
            outstream.write('%s\tIR\texon\t%d\t%d\t.\t%s\t.\tID=%s;Parent=%s\n' % (
                intron.chrom, sorted(intron.p3)[0], sorted(intron.p3)[1],
                intron.strand, '%s:IE:3pexon' %name, '%s:IE' %name))
        else:
             # write IE 3p exon record
            outstream.write('%s\tIR\texon\t%d\t%d\t.\t%s\t.\tID=%s;Parent=%s\n' % (
                intron.chrom, sorted(intron.p3)[0], sorted(intron.p3)[1],
                intron.strand, '%s:IE:3pexon' %name, '%s:IE' %name))

             # write IE 5p exon record
            outstream.write('%s\tIR\texon\t%d\t%d\t.\t%s\t.\tID=%s;Parent=%s\n' % (
                intron.chrom, sorted(intron.p5)[0], sorted(intron.p5)[1],
                intron.strand, '%s:IE:5pexon' %name, '%s:IE' %name))


class InferredIntron(object) :
    """Simple class that stores information about an inferred intron."""
    def __init__(self, chrom, strand, parent, child) :
        assert(type(parent) != int)
        self.parent = parent
        self.child  = child
        self.chrom  = chrom
        self.strand = strand
        self.p5 = (self.parent.acceptorEnd(), self.parent.donorEnd())
        self.p3 = (self.child.acceptorEnd(), self.child.donorEnd())
        self.irange = (self.parent.donorEnd(), self.child.acceptorEnd())
        self.minpos = min(self.irange)
        self.maxpos = max(self.irange)
        self.maxPexp = 0.0
        self.minPexp = 1.0
        self.stages = [ ]
        self.unresolved = False

    def adjust(self, other):
        # adjust 5p exon
        if abs(other.p5[0] - other.p5[1]) < abs(self.p5[0] - self.p5[1]):
            self.parent = other.parent
            self.p5 = other.p5
            self.irange = (self.parent.donorEnd(), self.child.acceptorEnd())
            self.minpos = min(self.irange)
            self.maxpos = max(self.irange)

        # adjust 3p exon
        if abs(other.p3[0] - other.p3[1]) < abs(self.p3[0] - self.p3[1]):
            self.child = other.child
            self.p3 = other.p3
            self.irange = (self.parent.donorEnd(), self.child.acceptorEnd())
            self.minpos = min(self.irange)
            self.maxpos = max(self.irange)
            

    def acceptor(self) :
        return self.child.acceptorEnd()

    def donor(self) :
        return self.parent.donorEnd()

    def __eq__(self, other) :
        return self.chrom == other.chrom and self.strand == other.strand and self.minpos == other.minpos and self.maxpos == other.maxpos

    def __hash__(self) :
        return self.__str__().__hash__()

    def __str__(self) :
        return "%s,%d,%d,%s" % (self.chrom, self.minpos, self.maxpos, self.strand)

def getGraphIntrons(graph, chrom) :
    """Uses a splice graph to infer all the introns in a gene."""
    strand   = graph.strand
    resIR    = set([ ])
    resIE    = set([ ])
    for node in graph.resolvedNodes():
        for c in node.children :
            intron = InferredIntron(chrom, strand, node, c)
            minpos = intron.minpos
            maxpos = intron.maxpos
            containers = [o for o in graph.resolvedNodes() if o.minpos < minpos - 10 and maxpos + 10 < o.maxpos]
            if containers :
                resIR.add(intron)
            else:
                resIE.add(intron)

    return resIR, resIE

def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', dest='model', default=SG_GENE_MODEL, help='Gene model GFF file [default: %(default)s]')
    parser.add_argument('-o', dest='outfile', default='out.gff', help='Output file base')
    parser.add_argument('-v', dest='verbose', default=False, help='Verbose mode [default: %(default)s]', action='store_true')
    return parser


def parse_args(argv=None):
    parser = build_parser()
    opts = parser.parse_args(argv)
    errStrings = []
    if not opts.model:
        errStrings.append('** No GFF gene model specified.  Set SPLICEGRAPHER_GENE_MODEL or use the -m option.')
    if errStrings:
        parser.print_help()
        sys.stderr.write('\n%s\n' % '\n'.join(errStrings))
        raise SystemExit(1)
    return opts, []


def main(argv=None):
    """Generate MISO intron-retention GFF records from a gene model."""
    opts, _ = parse_args(argv)
    gene_model = loadGeneModels(
        opts.model,
        verbose=opts.verbose,
        alltypes=True,
        outdir=Path(opts.outfile).parent,
    )

    irs = 0
    ies = 0
    with open(opts.outfile, 'w') as out_stream:
        for chrm in gene_model.getChromosomes():
            if opts.verbose:
                sys.stderr.write('Processing genes from chromosome: %s\n' % chrm)
            indicator_g = ProgressIndicator(10000, description=' genes', verbose=opts.verbose)
            genes = gene_model.getGeneRecords(chrm, geneFilter=gene_type_filter)
            genes.sort()

            for g in genes:
                chrom = g.gffString().split('\t')[0]
                if opts.verbose:
                    indicator_g.update()

                gene_graph = makeSpliceGraph(g)
                gene_graph.annotate()

                # get introns from gene models
                ir_gm, ie_gm = getGraphIntrons(gene_graph, chrom)
                writeIntrons(g, ir_gm, 'K', out_stream)
                writeIntrons(g, ie_gm, 'P', out_stream)

                irs += len(ir_gm)
                ies += len(ie_gm)

            if opts.verbose:
                indicator_g.finish()
                sys.stderr.write('%d genes\n' % indicator_g.ctr)

    if opts.verbose:
        sys.stderr.write('Found %d retained introns and %d excised introns\n' % (irs, ies))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
