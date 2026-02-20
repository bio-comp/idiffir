#!/usr/bin/env python

#
#    iDiffIR- identifying differential intron retention (IR) from RNA-seq
#    Copyright (C) 2015  Michael Hamilton
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

"""
Generate differential isoform reads

.. moduleauthor:: Mike Hamilton <hamiltom@cs.colostate.edu>
"""

from iDiffIR.SpliceGrapher.shared.utils        import *
from iDiffIR.SpliceGrapher.shared.GeneModelConverter import *
from iDiffIR.SpliceGrapher.formats.FastaLoader import FastaLoader
from iDiffIR.SpliceGrapher.formats.GeneModel   import *
from iDiffIR.SpliceGrapher.formats.fasta       import FastaRecord

from iDiffIR.SpliceGrapher.formats.loader import loadGeneModels

import argparse
import os, sys, warnings, numpy, itertools, random, math
from pathlib import Path

USAGE="""%(prog)s [options]

"""

opts = None
gene_model = None
depths = {}
genes = []
loader = None
diff_isos = []
diff_genes = []


def build_parser():
    parser = argparse.ArgumentParser(usage=USAGE)
    parser.add_argument('-f', dest='fasta', default=SG_FASTA_REF, help='FASTA reference file [default: %(default)s]')
    parser.add_argument('-m', dest='model', default=SG_GENE_MODEL, help='Gene model GFF file [default: %(default)s]')
    parser.add_argument('-n', dest='qname', default='DISimulator', help='QNAME for read alignment records [default: %(default)s]')
    parser.add_argument('-o', dest='outfile', default=None, help='Output file [default: %(default)s]')
    parser.add_argument('-g', dest='diffGExp', default=0.25, type=float, help='Frequency of differential gene expression [default: %(default)s]')
    parser.add_argument('-i', dest='diffIExp', default=0.10, type=float, help='Frequency of differential isoform expression [default: %(default)s]')
    parser.add_argument('-D', dest='diffIso', default=2.0, type=float, help='Relative differential isoform expression')
    parser.add_argument('-l', dest='read_length', default=80, type=int, help='Length of simulated reads [default: %(default)s]')
    parser.add_argument('-a', dest='minAnchor', default=8, type=int, help='Minimum anchor length for a spliced alignment [default: %(default)s]')
    parser.add_argument('-t', dest='test', default=False, action='store_true', help='Run test [default: %(default)s]')
    parser.add_argument('-v', dest='verbose', default=False, help='Verbose mode [default: %(default)s]', action='store_true')
    parser.add_argument('-d', dest='depth', default=50, type=int, help='Depth of up-regulated and equal expression')
    parser.add_argument('-e', dest='expfile', default=None, type=str, help='path of expression file')
    parser.add_argument('-s', dest='summary_file', default='summary.txt', type=str, help='name of summary file')
    parser.add_argument('-r', dest='reps', default=1, type=int, help='number of replicates for each condition')
    parser.add_argument('-p', dest='procs', default=1, type=int, help='number of processors to use for sorting sam files')
    return parser


def parse_args(argv=None):
    parser = build_parser()
    parsed_opts = parser.parse_args(argv)
    errStrings = []
    if not parsed_opts.model:
        errStrings.append('** No GFF gene model specified.  Set SPLICEGRAPHER_GENE_MODEL or use the -m option.')
    if not parsed_opts.fasta:
        errStrings.append('** No FASTA reference specified.  Set SPLICEGRAPHER_FASTA_REF or use the -f option.')
    if errStrings:
        parser.print_help()
        sys.stderr.write('\n%s\n' % '\n'.join(errStrings))
        raise SystemExit(1)
    return parsed_opts, []


def initialize_state(argv=None):
    """Load models and simulation state before generating reads."""
    global opts, gene_model, depths, genes, loader, diff_isos, diff_genes
    opts, _ = parse_args(argv)
    gene_model = loadGeneModels(
        opts.model,
        verbose=opts.verbose,
        alltypes=True,
        outdir=Path(opts.outfile).parent if opts.outfile else Path.cwd(),
    )
    depths = {}
    genes = gene_model.getAllGenes(geneFilter=gene_type_filter)
    if depths:
        genes = [g for g in genes if not g.isSingleExon() and g.id in depths]
    else:
        genes = [g for g in genes if not g.isSingleExon()]
    genes.sort()
    loader = FastaLoader(opts.fasta, verbose=opts.verbose)
    multi_isos = [gene.id for gene in genes if len(gene.isoforms) > 1]
    gene_ids = [gene.id for gene in genes]
    # select differential isoform genes
    diff_isos = random.sample(multi_isos, int(len(multi_isos) * opts.diffIExp))
    # select differential genes
    diff_genes = random.sample(gene_ids, int(len(gene_ids) * opts.diffGExp))


def loadExp( ):
    with open( opts.expfile, 'r') as fin:
        _ = fin.readline()
        for line in fin:
            line = line.strip()
            gid = line.split('\t')[0]
            mtg, wtg, mti, wti = line.split('\t')[1:]
            mtg = int(math.ceil(float(mtg)))
            wtg = int(math.ceil(float(wtg)))
            if mtg == 0 or wtg == 0: continue
            if wti == '-':
                depths[gid] = (max(1,mtg), max(1,wtg), 0)
            else:
                wti = float(wti)
                depths[gid] = ( max(1,mtg), max(1,wtg),  min(1.0,max(1.0,wti)/max(1,wtg) ))

def makeCIGAR( positions):
    cigar = ""
    match_run = 1
    prev = positions[ 0 ]
    for idx in positions[1:]:
        if idx > prev + 1:
            cigar += '%dM' % match_run
            cigar += '%dN' % (idx - prev - 1)
            match_run = 1
        else:
            match_run += 1
        prev = idx
    cigar += '%dM' % match_run
    return cigar

def makeSAMRecord( rid, chrom, strand, start, positions, sequence ):
    quality = 'F'
    cigar = makeCIGAR( positions )
    minAnchor = numpy.inf
    if 'N' in cigar:
        for x in cigar.split('N'):
            anchor= int(x.split('M')[0])
            if anchor < minAnchor:
                minAnchor = anchor

    record = '%s\t%d\t%s\t%d\t50\t%s\t*\t0\t0\t%s\t%s\tAS:i:0\tXN:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tNM:i:0\tMD:Z:%d\tNH:i:1' % (rid,
                                                                                                                        0 if strand =='+' else 16,
                                                                                                                        chrom, start, cigar,
                                                                                                                        sequence,
                                                                                                                        quality*len(sequence),
                                                                                                                        len(sequence))
    return record, minAnchor

def makeSAMHeader( loader, outstream ):
    """
    Write header for sam file.

    :param loader: FastaLoader instance
    :param outstream: sam file stream
    """
    outstream.write('@HD\tVN:1.0\tSO:coordinate\n')
    for sid in sorted(loader.keys()):
        outstream.write('@SQ\tSN:%s\tLN:%d\n' % ( sid, loader.sequenceLength(sid)))
    outstream.write('@PG\tID:DISimulator\tVN:1.0\tCL:%s\n' % ( ' '.join( sys.argv ) ) )

def simulateReads( positions, sequence, chrom, strand, depth, length, qname, outstream, minAnchor=8):
    """
    Simulate isoform reads.

    :param positions: list of positions of the CDS w.r.t genome reference
    :param sequence: coding sequence
    :param chrom: chromosome name
    :param strand: '+' or '-'
    :param depth: depth of the read coverage
    :param length: length of the reads
    :param qname: name of the read sample (prepended to each read number)
    """
    assert( len(positions) == len(sequence) )
    L = len(sequence)
    nReads = depth * int(L/float(length))
    #nReads = depth * L
    for _ in range( nReads ):
        simulateReads.n += 1
        rid = '%s.%d' % ( qname, simulateReads.n)
        idx = random.randint( 0, L-length )
        start = positions[idx]
        rec, minA = makeSAMRecord( rid, chrom, strand, start,
                             positions[idx:(idx+length)], sequence[idx:(idx+length)] )
        if 1:
            outstream.write('%s\n' % rec )

def run_test( ):
    alt_events = { }
    as_map = {}
    gene_dict = {}
    for gene in genes:
        gene_dict[gene.id] = gene
    gene = gene_dict[diff_isos[0]]
    print(gene.id, gene.strand)

    g = makeSpliceGraph(gene)
    g.annotate()
    for ev in g.distinctAltEvents():
        key = tuple(sorted( [ev[0], ev[1]]))
        l = as_map.get(key, set([]))
        l.add( ev[2] )
        as_map[key] = l
    if opts.verbose:
        sys.stderr.write('Creating differential isoform reads for %s\n' %  gene.id )
    diffIso_n  = random.choice( range( len( gene.isoforms)))
    mt_exp     = numpy.random.poisson( 50, 1 )
    mt_diffExp = int(mt_exp * 0.75)
    mt_eqExp   = max(1, int((mt_exp*0.25)/ (len(gene.isoforms)-1) ))
    mt_IsoDepths = [ mt_diffExp if i == diffIso_n else mt_eqExp for i in range( len(gene.isoforms ) ) ]

    wt_exp     = numpy.random.poisson( 50, 1 )
    wt_eqExp   = max(1, int(wt_exp/ float(len(gene.isoforms) )))


    wt_IsoDepths = [ wt_eqExp ] * len( gene.isoforms )

    for i,iso in enumerate(gene.isoforms):
        exons = sorted( [sorted( [x.start(), x.end()] ) for x in gene.isoforms.values()[i].sortedExons() ] )
        sequence = ''.join([ loader.subsequence( gene.chromosome,
                                                 exon[0]-1, exon[1]-1,
                                                 reverse=gene.strand=='-') for exon in exons])
        positions = list(itertools.chain( *[range( x[0], x[1]+1) for x in exons ] ))
        simulateReads( positions, sequence, gene.chromosome, gene.strand,
                       mt_IsoDepths[i], opts.read_length, 'mutant_1', mt_OutStream)
        simulateReads( positions, sequence, gene.chromosome, gene.strand,
                       wt_IsoDepths[i], opts.read_length, 'wildtype_1', wt_OutStream)
    mt_OutStream.close()
    wt_OutStream.close()

    os.system("""awk '$1 ~ /^@/{print $0;next}{print $0 | "sort -k 3,3 -k 4,4n"}' wildtype.sam > tmpFile && mv tmpFile wildtype.sam""")
    os.system("""awk '$1 ~ /^@/{print $0;next}{print $0 | "sort -k 3,3 -k 4,4n"}' mutant.sam > tmpFile && mv tmpFile mutant.sam""")
    evs = []
    for start,end in sorted( [sorted( [x.start(), x.end()] ) for x in gene.isoforms.values()[i].sortedExons() ] ):
        key2 = (start,end)
        for key in as_map:
            if key == key2:
                print(key, as_map[key])
            elif key2[0] >= key[0] and key2[1] >= key[1]:
                print(key, as_map[key])
    print(as_map)
    template = \
               """[SinglePlotConfig]
legend        = True
output_file   = %s.pdf
height        = 8.0
width         = 16.0
fontsize      = 16
shrink_introns= False

[GeneModel]
plot_type     = gene
hide          = False
source_file   = %s
file_format   = gene_model
gene_name     = %s
title_string  = Gene Model for %s

[CoverageMT]
plot_type     = read_depth
source_file   = mutant.sam
title_string  = Mutant
log_threshold = 100
labels        = False

[JunctionsMT]
plot_type     = junctions
source_file   = mutant.sam
title_string  = Mutant
log_threshold = 100
labels        = True


[CoverageWT]
plot_type     = read_depth
source_file   = wildtype.sam
title_string  = Wildtype
log_threshold = 100
labels        = False

[JunctionsWT]
plot_type     = junctions
source_file   = wildtype.sam
title_string  = Wildtype
log_threshold = 100
labels        = True
               """ % (gene.id, opts.model, gene.id, gene.id )
    with open('plot.config', 'w') as fout:
        fout.write( '%s' % template )
    os.system( 'plotter.py plot.config' )

def makeSingleIsoform( gene, depths1, label1, outstream1, depths2, label2, outstream2, minAnchor ):
    exons = sorted( [sorted( [x.start(), x.end()] ) for x in gene.isoforms.values()[0].sortedExons() ] )
    sequence = ''.join([ loader.subsequence( gene.chromosome,
                                             exon[0]-1, exon[1]-1,
                                             reverse=gene.strand=='-') for exon in exons])
    positions = list(itertools.chain( *[range( x[0], x[1]+1) for x in exons ] ))
    simulateReads( positions, sequence, gene.chromosome, gene.strand,
                   depths1, opts.read_length, label1, outstream1, minAnchor)
    simulateReads( positions, sequence, gene.chromosome, gene.strand,
                   depths2, opts.read_length, label2, outstream2, minAnchor)

def makeMultiIsoform( gene,  depths1, label1, outstream1, depths2, label2, outstream2, minAnchor ):
    for i, iso in enumerate(gene.isoforms):
        exons = sorted( [sorted( [x.start(), x.end()] ) for x in gene.isoforms.values()[i].sortedExons() ] )
        sequence = ''.join([ loader.subsequence( gene.chromosome,
                                                 exon[0]-1, exon[1]-1,
                                                 reverse=gene.strand=='-') for exon in exons])
        positions = list(itertools.chain( *[range( x[0], x[1]+1) for x in exons ] ))
        simulateReads( positions, sequence, gene.chromosome, gene.strand,
                       depths1[i], opts.read_length, label1, outstream1, minAnchor)
        simulateReads( positions, sequence, gene.chromosome, gene.strand,
                       depths2[i], opts.read_length, label2, outstream2, minAnchor)

def annotateIRs(gene):
    isoIRs = [ ]
    exons = set([( min(e.start(), e.end()), max(e.start(), e.end())) for e in gene.sortedExons()])
    for iso in gene.isoforms.values():
        hasIR = False
        introns = [ sorted(x) for x in iso.sortedIntrons() ]
        for s,e in introns:
            for s2, e2 in exons:
                if s >= s2 and e <= e2:
                    hasIR = True
                    break
        isoIRs.append(hasIR)
    return isoIRs



def main():
    summary_stream = open(opts.summary_file, 'w')
    summary_stream.write( 'geneID\tmt_depth\twt_depth\tDE\tDAS\tEvents\n')
    IR_events = { }
    mt_OutStreams = [ ]
    wt_OutStreams = [ ]
    for r in range(opts.reps):
        mt_OutStreams.append( open('mutant_%d.sam' % (r+1), 'w'))
        wt_OutStreams.append( open('wildtype_%d.sam' % (r+1), 'w'))

        makeSAMHeader(loader, mt_OutStreams[-1])
        makeSAMHeader(loader, wt_OutStreams[-1])

    if opts.expfile:
        loadExp()
    indicator = ProgressIndicator(10000, description=' genes', verbose=opts.verbose)
    for gene in genes:

        indicator.update()
        # use given expressions
        if depths:
            if gene.id not in depths:continue
            mte, wte, ire = depths[gene.id]
            ce = 1.0 - ire
            if mte > wte:
                up_outstreams = mt_OutStreams
                up_label = 'mutant'
                up_depth = int(mte)
                dn_outstreams = wt_OutStreams
                dn_label = 'wildtype'
                dn_depth = int(wte)

            # pick wildtype as up-regulated
            else:
                dn_outstreams = mt_OutStreams
                dn_label = 'mutant'
                dn_depth = int(mte)
                up_outstreams = wt_OutStreams
                up_label = 'wildtype'
                up_depth = int(wte)

            summary_stream.write( '%s\t%d\t%d\t%s\t' % (gene.id, up_depth, dn_depth, 'True' ))
        #otherwise pick random expression
        else:
        # differential gene expression
            if gene.id in diff_genes:
                # pick mutant as up-regulated
                if random.randint(0,1):
                    up_outstreams = mt_OutStreams
                    up_label = 'mutant'
                    up_depth = numpy.random.poisson( opts.depth, 1 )
                    dn_outstreams = wt_OutStreams
                    dn_label = 'wildtype'
                    dn_depth = numpy.random.poisson( int(opts.depth/2.0), 1 )
                    summary_stream.write( '%s\t%d\t%d\t%s\t' % (gene.id, up_depth, dn_depth, 'True' ))
                # pick wildtype as up-regulated
                else:
                    dn_outstreams = mt_OutStreams
                    dn_label = 'mutant'
                    dn_depth = numpy.random.poisson( int(opts.depth/2.0), 1 )
                    up_outstreams = wt_OutStreams
                    up_label = 'wildtype'
                    up_depth = numpy.random.poisson( opts.depth, 1 )
                    summary_stream.write( '%s\t%d\t%d\t%s\t' % (gene.id, dn_depth, up_depth, 'True' ))
            # equal gene expression
            else:
                up_outstreams = mt_OutStreams
                up_label = 'mutant'
                up_depth = numpy.random.poisson( opts.depth, 1 )
                dn_outstreams = wt_OutStreams
                dn_label = 'wildtype'
                dn_depth = numpy.random.poisson( opts.depth, 1 )
                summary_stream.write( '%s\t%d\t%d\t%s\t' % (gene.id, up_depth, dn_depth, 'False' ))
        # only one isoform
        if len( gene.isoforms ) < 2:
            for r in range(opts.reps):
                makeSingleIsoform( gene, up_depth, up_label, up_outstreams[r], dn_depth, dn_label, dn_outstreams[r], opts.minAnchor )
            summary_stream.write('False\t\n')
            continue
        # multi-isoforms
        as_map = {}
        g = makeSpliceGraph(gene)
        g.annotate()
        irIsos = annotateIRs(gene)

        for ev in g.distinctAltEvents():
            key = tuple(sorted( [ev[0], ev[1]]))
            l = as_map.get(key, set([]))
            l.add( ev[2] )
            as_map[key] = l
        diffIso_n  = random.choice( range( len( gene.isoforms)))
        diso = gene.isoforms.values()[diffIso_n]
        # get differntial IR events, if there are any in the selected isoform
        evs = set([])
        diso_start, diso_end = sorted( (diso.start(), diso.end()) )
        did = list(gene.isoforms.keys())[diffIso_n]
        for start,end in sorted( [sorted( [x.start(), x.end()] ) for x in diso.sortedExons() ] ):
            key2 = (start,end)
            for key in as_map:
                if key == key2:
                    if 'IR' in as_map[key]:
                        evs.add( (start, end, 'P') )

        for j in range( len ( gene.isoforms ) ):
            if j == diffIso_n: continue
            for start,end in sorted( [sorted( [x.start(), x.end()] ) for x in gene.isoforms.values()[j].sortedExons() ] ):
                key2 = (start,end)
                for key in as_map:
                    if key == key2:
                        if 'IR' in as_map[key] and (start >= diso_start and start <= diso_end or end >= diso_start and end <= diso_end):
                            evs.add( (start, end, 'S') )

        lambdasEq = numpy.array([ max(0.15,ire) if irIsos[i] else max(0.15,1.0-ire) for i in range( len( gene.isoforms ) ) ])
        lambdasEq = lambdasEq/sum(lambdasEq)
        isoScalar = random.uniform(1.5,5)
        if random.randint(0,1):
            isoScalar = 1.0 / isoScalar
        lambdasDiff = numpy.array( [ isoScalar*lambdasEq[i] if i == diffIso_n else lambdasEq[i] for i in range( len(gene.isoforms))])
        lambdasDiff = lambdasDiff/sum(lambdasDiff)

        if gene.id in diff_isos:
            # pick up-regulated condition as differential isoform-containing condition
            if random.randint(0,1):
                diso_label = up_label
                diso_outstreams = up_outstreams
                diso_depths = [ int(math.ceil(x)) for x in lambdasDiff*up_depth]
                eiso_label = dn_label
                eiso_outstreams = dn_outstreams
                eiso_depths = [ int(math.ceil(x)) for x in lambdasEq*dn_depth]
                summary_stream.write('True\t%s\t' % (did))
                summary_stream.write( '%s\n' % ( ';'.join( ['%d,%d,%s' % ( start, stop, t) for start, stop, t in evs ] ) ) )


            else:
                diso_label = dn_label
                diso_outstreams = dn_outstreams
                diso_depths = [ int(math.ceil(x)) for x in lambdasDiff*dn_depth]
                eiso_label = up_label
                eiso_outstreams = up_outstreams
                eiso_depths = [ int(math.ceil(x)) for x in lambdasEq*up_depth]
                summary_stream.write('True\t%s\t' % (did))
                summary_stream.write( '%s\n' % ( ';'.join( ['%d,%d,%s' % ( start, stop, t) for start, stop, t in evs ] ) ) )
            summary_stream.flush()
        else:
            diso_label = up_label
            diso_outstreams = up_outstreams
            diso_depths =  [ int(x) for x in lambdasEq*up_depth]
            eiso_label = dn_label
            eiso_outstreams = dn_outstreams
            eiso_depths =  [ int(x) for x in lambdasEq*dn_depth]
            summary_stream.write('False\t\n')
        for r in range(opts.reps):
            makeMultiIsoform( gene,  diso_depths, diso_label, diso_outstreams[r], eiso_depths, eiso_label, eiso_outstreams[r], opts.minAnchor )
    indicator.finish()
    for r in range(opts.reps):
        mt_OutStreams[r].close()
        wt_OutStreams[r].close()
    summary_stream.close()
    # sort sam files
#    for r in range(opts.reps):
#        os.system("""awk '$1 ~ /^@/{print $0;next}{print $0 | "sort --parallel=%d -k 3,3 -k 4,4n"}' wildtype_%d.sam > tmpFile && mv tmpFile wildtype_%d.sam""" % (opts.procs, r+1, r+1))
#        os.system("""awk '$1 ~ /^@/{print $0;next}{print $0 | "sort --parallel=%d -k 3,3 -k 4,4n"}' mutant_%d.sam > tmpFile && mv tmpFile mutant_%d.#sam""" % (opts.procs, r+1, r+1))

if __name__ == '__main__':
    initialize_state()
    simulateReads.n = 0
    if opts.test:
        run_test()
    else:
        main()
