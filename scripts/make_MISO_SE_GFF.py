#! /usr/bin/env python
"""Program for generating known splice junction sequences."""
from iDiffIR.SpliceGrapher.shared.utils        import *
from iDiffIR.SpliceGrapher.shared.GeneModelConverter import *
from iDiffIR.SpliceGrapher.predict.SpliceSite  import *
from iDiffIR.SpliceGrapher.formats.GeneModel   import *
from iDiffIR.SpliceGrapher.SpliceGraph         import *
from iDiffIR.SpliceGrapher.formats.loader import loadGeneModels

from glob     import glob
import argparse
import os, sys, warnings
from pathlib import Path

def build_parser():
    parser = argparse.ArgumentParser(description='Generate GFF file for MISO of exon skipping events')
    parser.add_argument('-m', dest='model', default=SG_GENE_MODEL, help='Gene model GFF file [default: %(default)s]')
    parser.add_argument('-o', dest='outfile', default=None, help='Output file [default: %(default)s]')
    parser.add_argument('-v', dest='verbose', default=False, help='Verbose mode [default: %(default)s]', action='store_true')
    parser.add_argument('-s', dest='graphPaths', default=None, help='File containing paths to splice graphs to augment gene models')
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

def getEventLocs( node ):
    """Return parent/cassette/child exon tuples for an SE node."""
    if node.strand == '+':
        parents = list(set([( e.minpos, e.maxpos) for e in node.parents]))
        children = list(set([( e.minpos, e.maxpos) for e in node.children]))
    else:
        parents = list(set([( e.minpos, e.maxpos) for e in node.children]))
        children = list(set([( e.minpos, e.maxpos) for e in node.parents]))

    return [ (parents[i],(node.minpos,node.maxpos), children[j])\
                 for i in range(len(parents)) for j in range(len(children))]


def load_splice_graphs(graph_paths, verbose):
    """Load optional auxiliary splice graphs keyed by upper-case gene name."""
    splice_graphs = {}
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
                gene_name = graph.name.upper()
                if gene_name in splice_graphs:
                    splice_graphs[gene_name] = splice_graphs[gene_name].union(graph)
                else:
                    splice_graphs[gene_name] = graph
        indicator.finish()
        if verbose:
            sys.stderr.write('Loaded %d aux. splice graphs\n' % indicator.ctr)
    return splice_graphs


def main(argv=None):
    """Generate MISO SE GFF records from gene models and optional graph overlays."""
    opts, _ = parse_args(argv)
    gene_model = loadGeneModels(
        opts.model,
        verbose=opts.verbose,
        alltypes=True,
        outdir=Path(opts.outfile).parent if opts.outfile else Path.cwd(),
    )
    genes = gene_model.getAllGenes(geneFilter=gene_type_filter)
    genes.sort()
    splice_graphs = load_splice_graphs(opts.graphPaths, opts.verbose)

    total_recs = 0
    unique_events = 0
    chr_counter = {}
    indicator = ProgressIndicator(10000, description=' genes', verbose=opts.verbose)

    out_stream = open(opts.outfile, 'w') if opts.outfile else sys.stdout
    try:
        for gene in genes:
            if gene.chromosome not in chr_counter:
                chr_counter[gene.chromosome] = 0

            # Create splice graph for gene:
            try:
                gene_graph = makeSpliceGraph(gene)
            except ValueError:
                sys.stderr.write("Unable to create graph for %s\n" % gene.name)
                continue

            if gene_graph.name.upper() in splice_graphs:
                gene_graph = gene_graph.union(splice_graphs[gene_graph.name.upper()])
            gene_graph.annotate()

            if not gene_graph.hasAS():
                continue
            if 'SE' not in gene_graph.altForms():
                continue

            se_nodes = [node for node in gene_graph.resolvedNodes() if node.isSkippedExon()]

            for node_num, node in enumerate(se_nodes):
                unique_events += 1
                node_intron_count = 0
                for exon5, exonS, exon3 in getEventLocs(node):
                    total_recs += 1
                    node_intron_count += 1
                    event_id = '%s:%d-%d:%d:%d' % (gene.id, exonS[0], exonS[1], node_num + 1, node_intron_count)
                    out_stream.write('%s\tIR\tgene\t%d\t%d\t.\t%s\t.\tID=%s;Name=%s\n' % (
                        gene.chromosome, exon5[0], exon3[1], node.strand, event_id, event_id))

                    ir_id = '%s:%d-%d:%d:%d' % (gene.id, exonS[0], exonS[1], node_num + 1, node_intron_count)
                    # make mRNA record for SE
                    out_stream.write('%s\tSE\tmRNA\t%d\t%d\t.\t%s\t.\tID=%s:SE;Parent=%s\n' % (
                        gene.chromosome, exon5[0], exon3[1], node.strand, ir_id, event_id))
                    # make mRNA record for CS
                    out_stream.write('%s\tSE\tmRNA\t%d\t%d\t.\t%s\t.\tID=%s:CS;Parent=%s\n' % (
                        gene.chromosome, exon5[0], exon3[1], node.strand, ir_id, event_id))
                    # make 5' exon record for SE
                    out_stream.write('%s\tSE\texon\t%d\t%d\t.\t%s\t.\tID=%s:se5Exon;Parent=%s:SE\n' % (
                        gene.chromosome, exon5[0], exon5[1], node.strand, ir_id, ir_id))
                    # make 3' exon record for SE
                    out_stream.write('%s\tSE\texon\t%d\t%d\t.\t%s\t.\tID=%s:se3Exon;Parent=%s:SE\n' % (
                        gene.chromosome, exon3[0], exon3[1], node.strand, ir_id, ir_id))

                    # make 5' exon record for CS
                    out_stream.write('%s\tSE\texon\t%d\t%d\t.\t%s\t.\tID=%s:se5Exon;Parent=%s:CS\n' % (
                        gene.chromosome, exon5[0], exon5[1], node.strand, ir_id, ir_id))
                    # make 3' exon record for CS
                    out_stream.write('%s\tSE\texon\t%d\t%d\t.\t%s\t.\tID=%s:se3Exon;Parent=%s:CS\n' % (
                        gene.chromosome, exon3[0], exon3[1], node.strand, ir_id, ir_id))
                    # make cassette exon record for CS
                    out_stream.write('%s\tSE\texon\t%d\t%d\t.\t%s\t.\tID=%s:cassExon;Parent=%s:CS\n' % (
                        gene.chromosome, exonS[0], exonS[1], node.strand, ir_id, ir_id))
    finally:
        if out_stream is not sys.stdout:
            out_stream.close()

    if opts.verbose:
        indicator.finish()

    sys.stderr.write("Wrote %d total SE events\n" % total_recs)
    sys.stderr.write("Wrote %d unique SE events\n" % unique_events)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
