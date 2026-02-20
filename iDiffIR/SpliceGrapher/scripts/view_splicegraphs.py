#! /usr/bin/env python
# Copyright (C) 2010 by Colorado State University
# Contact: Mark Rogers <rogersma@cs.colostate.edu>
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
"""
Script for viewing 1 or more splice graphs.
"""
from iDiffIR.SpliceGrapher.shared.config             import *
from iDiffIR.SpliceGrapher.shared.adjust             import *
from iDiffIR.SpliceGrapher.SpliceGraph               import *
from iDiffIR.SpliceGrapher.shared.utils              import *
from iDiffIR.SpliceGrapher.shared.GeneModelConverter import *

from pylab      import *
import argparse
from sys import maxsize as MAXINT
import sys, os

FILE_TAG    = 'file'
GENE_TAG    = 'gene'
DEFAULT_PDF = 'view_splicegraphs.pdf'

#==========================================================================
# Incompatibilities in some matplotlib versions may yield runtime warnings
import warnings
warnings.filterwarnings('ignore')
#==========================================================================

#==========================================================================
from iDiffIR.SpliceGrapher.view.ViewerUtils  import *
from iDiffIR.SpliceGrapher.view.GeneView     import GeneView
from iDiffIR.SpliceGrapher.formats.loader    import *
#==========================================================================

EXT_TO_DISPLAY = {'.gff':PREDICTED_GRAPH}
GZIP_EXTS      = ['.gz', '.gzip']

# Arbitrary Y limit for splice graphs
SG_Y_LIMIT     = 100.0
X_PAD_FRACTION = 0.01

#==========================================================================
# Main program
#==========================================================================
USAGE = """%(prog)s [options] files-or-genes

Displays one or more splice graphs given a list of files or genes.  A name
that is not a file is assumed to be a gene.  When gene names are given,
gene models must also be provided (using SG_GENE_MODEL or the -m option).

Examples:

    %(prog)s gene_1.gff                    (plot a single graph file)
    %(prog)s gene*.gff                     (display all graphs in the local directory that begin with 'gene')
    %(prog)s AT1G01060 AT1G01260 AT1G01910 (plot several genes)
    %(prog)s AT1G01060 gene_[123].gff      (mixture of gene model graphs and splice graph files)

For more than four graphs you may need to adjust the output height (-H)."""

#==========================================================================
# Initialize command-line options:
parser = argparse.ArgumentParser(usage=USAGE)
parser.add_argument('-a', dest='annotate', default=False,         help='(Re)annotate graphs [default: %(default)s]', action='store_true')
parser.add_argument('-m', dest='model',    default=SG_GENE_MODEL, help='GFF gene model reference [default: %(default)s]')
parser.add_argument('-o', dest='output',   default=None,          help='Output file [default: screen]')
parser.add_argument('-x', dest='xLabel',   default=False,         help='Show genomic positions [default: %(default)s]', action='store_true')
parser.add_argument('-v', dest='verbose',  default=False,         help='Use verbose output [default: %(default)s]', action='store_true')
parser.add_argument('-E', dest='edge',     default=1,             help='Intron edge weight [default: %(default)s]', type=int)
parser.add_argument('-L', dest='legend',   default=False,         help='Add legend to splice graph plot [default: %(default)s]', action='store_true')
parser.add_argument('-F', dest='fontsize', default=12,            help='Font size for plot titles [default: %(default)s]', type=int)
parser.add_argument('-H', dest='height',   default=8,             help='Display height, in inches [default: %(default)s]', type=float)
parser.add_argument('-W', dest='width',    default=15,            help='Display width, in inches [default: %(default)s]', type=float)
parser.add_argument('-S', dest='shrink',   default=False,         help='Shrink introns [default: %(default)s]', action='store_true')
parser.add_argument('-X', dest='xrange',   default=False,         help='Use the same genomic position range for all plots [default: %(default)s]', action='store_true')

# Deprecated to simplify interface:
#parser.add_argument('-t', dest='titles',   default=None,          help='List of graph titles [default: %(default)s]')
#parser.add_argument('-A', dest='adjust',   default=False,         help='Adjust graphs to match gene boundaries [default: %(default)s]', action='store_true')
#parser.add_argument('-G', dest='grid',     default=False,         help='Add axis grids to plots [default: %(default)s]', action='store_true')
## parser.add_argument('-U', dest='urmargin', default=0,             help='Margin for subsuming unresolved nodes into resolved exons [default: %(default)s]', type=int)
#parser.add_argument('-l', dest='labels',   default=False,         help='Show exon labels on plots [default: %(default)s]', action='store_true')
#parser.add_argument('--shrink-factor', dest='shrinkfactor', default=MIN_INTRON_SIZE, help='Factor for shrinking introns (-S option) [default: %(default)s]', type=int)

#==========================================================================
# Process command-line options:
def _parse_opts_and_args(parser, argv):
    parser.add_argument('args', nargs='*')
    opts = parser.parse_args(argv)
    args = opts.args
    delattr(opts, 'args')
    return opts, args

opts, args = _parse_opts_and_args(parser, sys.argv[1:])
#----------------------------------
# Deprecated to simplify interface:
opts.adjust       = False
opts.grid         = False
opts.labels       = False
opts.titles       = None
opts.urmargin     = sys.maxsize
opts.shrinkfactor = MIN_INTRON_SIZE
#----------------------------------

if len(args) < 1 :
    parser.print_help()
    sys.stderr.write('No files given.\n')
    sys.exit(1)

if opts.model : validateFile(opts.model)

# Infer the type for each plot based on whether a file exists with that name.
graphType = {}
for f in args :
    graphType[f] = FILE_TAG if os.path.exists(f) else GENE_TAG

writeStartupMessage()

titles = None
if opts.titles :
    titles = opts.titles.split(',')
    if len(titles) != len(args) :
        raise ValueError('You must enter as many titles as there are files.')

geneModel = None
if opts.model :
    geneModel = loadGeneModels(opts.model, verbose=opts.verbose)

#==========================================================================
# Draw the graphs:
# Standard widths for all subplots
axLeft      = 0.05
axWidth     = 0.90
# Reference heights for subplots
pageTop     = 0.98
totalHeight = 0.85 if opts.legend else 0.90

# Initialize matplotlib values:
rcParams['font.size']      = opts.fontsize
rcParams['font.weight']    = 'bold'
rcParams['figure.figsize'] = opts.width, opts.height

numPlots    = len(args)
graphHeight = totalHeight/numPlots
dispHeight  = 0.8*graphHeight

#--------------------------------------------------
# Plot all splice graphs in the order given
patchDict  = {}
numPlotted = 0
allAxes    = {}
allMinpos  = MAXINT
allMaxpos  = 0
origMinpos = MAXINT
origMaxpos = 0
for i in range(len(args)) :
    f = args[i]
    if graphType[f] == FILE_TAG :
        g = getFirstGraph(f)
    else :
        gene = geneModel.getGeneByName(f)
        if gene :
            g = geneModelToSpliceGraph(gene, minexon=1, minintron=4, verbose=opts.verbose)
        else :
            sys.stderr.write('** No gene model found for %s; skipping\n' % f)
            continue

    if opts.annotate : g.annotate()

    if opts.shrink :
        origGraph = g
        ranges    = getGraphRanges(g, scaleFactor=opts.shrinkfactor)
        g         = adjustGraph(g, ranges=ranges)
        if opts.verbose : sys.stderr.write('Shrinking introns in graph %s; old length %d/new length %d\n' % (g.name, len(origGraph), len(g)))

    title   = titles[i] if titles else '%s (%s)' % (g.name, g.strand)
    offset  = (i+1)*graphHeight
    curAxes = axes([axLeft, pageTop-offset, axWidth, dispHeight])

    # Establish graph X boundaries (plus small padding) based on gene
    padding = int(X_PAD_FRACTION*(g.maxpos-g.minpos))

    # Try to place gene view in background of all graphs
    adjustment  = 0
    nearbyGenes = []
    if geneModel :
        gene = geneModel.getGeneByName(g.name)
        if gene :
            if opts.adjust :
                adjustment = gene.minpos-g.minpos
                g.adjust(adjustment)

            if opts.shrink :
                gene = adjustGene(gene, ranges=ranges)
                nearbyGenes = [gene]
            else :
                nearbyGenes = geneModel.getGenesInRange(gene.chromosome, g.minpos, g.maxpos, strand=gene.strand)

            gv = GeneView(nearbyGenes, curAxes)
            gv.plot()

    if adjustment != 0 and opts.verbose :
        sys.stderr.write('Adjusting graph in %s by %d\n' % (f,adjustment))

    # Plot view
    tmpPatches,extraPatches = plotSpliceGraph(g, curAxes,
                                              labels=opts.labels,
                                              xLabels=opts.xLabel,
                                              minwidth=opts.edge,
                                              title=title,
                                              unresolved=False,
                                              urmargin=opts.urmargin,
                                              adjustment=adjustment,
                                              genes=nearbyGenes)

    if tmpPatches :
        try :
            patchDict.update(tmpPatches)
        except Exception :
            sys.stderr.write('** Warning: invalid legend labels for graph %d\n' % i)

    curAxes.grid(opts.grid)
    if g.strand == '+' :
        curAxes.set_xlim(g.minpos-padding, g.maxpos+padding)
    else :
        curAxes.set_xlim(g.maxpos+padding, g.minpos-padding)

    allMinpos = min(allMinpos, g.minpos)
    allMaxpos = max(allMaxpos, g.maxpos)

    # Display x positions only on last (bottom) plot
    if opts.shrink :
        origMinpos = min(origMinpos, origGraph.minpos)
        origMaxpos = max(origMaxpos, origGraph.maxpos)
        curAxes.set_xticks([g.minpos,g.maxpos])
        curAxes.set_xticklabels(['%d'%x for x in [origGraph.minpos,origGraph.maxpos]])
    else :
        xvalues = setXticks(g.minpos-padding, g.maxpos+padding)
        curAxes.set_xticks(xvalues)
        curAxes.set_xticklabels(['%d'%x for x in xvalues])

    # Y labels are meaningless for splice graph plots
    curAxes.set_yticklabels([])
    allAxes[curAxes] = g
    numPlotted += 1

if opts.xrange :
    # Reset X range to be the same for all graphs
    padding = int(X_PAD_FRACTION*(allMaxpos-allMinpos))
    minPos  = allMinpos-padding
    maxPos  = allMaxpos+padding
    xvalues = setXticks(minPos, maxPos)
    for a in allAxes :
        if allAxes[a].strand == '+' :
            a.set_xlim(minPos, maxPos)
        else :
            a.set_xlim(maxPos, minPos)
        a.set_xticks(xvalues)
        if opts.shrink :
            a.set_xticklabels(['%d'%x for x in [origMinpos,origMaxpos]])
        else :
            a.set_xticklabels(['%d'%x for x in xvalues])

if numPlotted == 0 :
    sys.stderr.write('** No graphs were plotted; exiting\n')
    sys.exit(1)

if opts.legend and patchDict :
    # Values found via trial and error:
    rcParams['legend.axespad']       = 0.01
    rcParams['legend.handlelen']     = 0.02
    rcParams['legend.handletextsep'] = 0.01
    rcParams['legend.labelsep']      = 0.008
    rcParams['legend.fontsize']      = max(8,0.75*opts.fontsize)
    numCols = max(3,len(patchDict))
    figlegend(patchDict.values(), patchDict.keys(), 'lower center', ncol=numCols)

if opts.output :
    savefig(opts.output, dpi=100)
else :
    show()
close()
