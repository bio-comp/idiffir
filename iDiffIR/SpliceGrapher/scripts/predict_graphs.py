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
from iDiffIR.SpliceGrapher.shared.config     import *
from iDiffIR.SpliceGrapher.shared.utils      import *
from iDiffIR.SpliceGrapher.shared.ShortRead  import *
from iDiffIR.SpliceGrapher.formats.loader    import *
from iDiffIR.SpliceGrapher.formats.alignment_io       import *
from iDiffIR.SpliceGrapher.predict.SpliceGraphPredictor import *

import argparse
import os,sys

def depthPredictions(chromosomes, dDict, jDict, model, opts) :
    """Run predictions based on alignment information from a SpliceGrapher depths file."""
    nGenes = nPred  = 0
    sys.stderr.write('Generating predictions:\n')
    chromList = sorted(chromosomes)
    for c in chromList :
        sys.stderr.write('  starting chromosome %s' % c)
        geneList  = model.getGeneRecords(c)
        nGenes   += len(geneList)
        geneList.sort()
        depths = dDict[c] if c in dDict else []
        jcts   = jDict[c] if c in jDict else []
        outDir = os.path.join(opts.outdir, c)
        if not os.path.isdir(outDir) : os.makedirs(outDir)
        nPred += predictGraphs(geneList, depths, jcts, outDir, depthThreshold=opts.mindepth, novelOnly=opts.novel)
    return nPred, nGenes

def getCommonChromosomes(geneSet, depthSet, depthsFile, opts) :
    result = geneSet & depthSet
    if result :
        sys.stderr.write('Found %d chromosomes in common between %s and %s:\n' % \
                (len(result), opts.model, depthsFile))
        sys.stderr.write('  %s\n' % ', '.join(sorted(result)))
    else :
        sys.stderr.write('Found no chromosomes in common between %s and %s:\n' % (opts.model, depthsFile))
        sys.stderr.write('  Gene models: %s\n' % ', '.join(sorted(geneSet)))
        sys.stderr.write('  Depths file: %s\n' % ', '.join(sorted(depthSet)))
        sys.exit(1)
    return result

USAGE = """%(prog)s SAM-file [options]

Predicts splice graphs for all genes based on evidence in the given SAM file.
Predicted graphs will be put into separate files under the given output directory,
creating a subdirectory for each chromosome."""

# Establish command-line options:
parser = argparse.ArgumentParser(usage=USAGE)
parser.add_argument('-a', dest='anchor',   default=8,   help='Minimum anchor required for junctions [default: %(default)s]', type=int)
parser.add_argument('-d', dest='outdir',   default='.', help='Output file [default: current directory]')
parser.add_argument('-j', dest='minjct',   default=2,   help='Minimum coverage required for splice junctions [default: %(default)s]', type=int)
parser.add_argument('-D', dest='mindepth', default=1,   help='Minimum average read depth required for clusters [default: %(default)s]', type=float)
parser.add_argument('-m', dest='model',    default=SG_GENE_MODEL, help='Gene models file [default: %(default)s]')
parser.add_argument('-N', dest='novel',    default=False, help='Only write splice graphs that were modified [default: %(default)s]', action='store_true')
parser.add_argument('-v', dest='verbose',  default=False, help='Verbose mode [default: %(default)s]', action='store_true')
def _parse_opts_and_args(parser, argv):
    parser.add_argument('args', nargs='*')
    opts = parser.parse_args(argv)
    args = opts.args
    delattr(opts, 'args')
    return opts, args

opts, args = _parse_opts_and_args(parser, sys.argv[1:])
MIN_ARGS = 1
if len(args) != MIN_ARGS :
    parser.print_help()
    if args : sys.stderr.write('\nExpected %d parameters; received %d:\n  %s\n' % (MIN_ARGS, len(args), '\n  '.join(args)))
    sys.exit(1)

if not opts.model :
    parser.print_help()
    if args : sys.stderr.write('\nYou must supply gene models (-m option) to make predictions\n')
    sys.exit(1)

dataFile = args[0]
validateFile(dataFile)
validateFile(opts.model)

# Load reference information from gene models
model        = loadGeneModels(opts.model, verbose=opts.verbose)
geneChromSet = set(model.getChromosomes())

if isDepthsFile(dataFile) :
    if opts.verbose : sys.stderr.write('Loading depth information from %s\n' % dataFile)
    depthDict,jctDict = readDepths(dataFile, minjct=opts.minjct, minanchor=opts.anchor, verbose=opts.verbose)
    depthChromSet     = set(depthDict.keys())
    commonChromSet    = getCommonChromosomes(geneChromSet, depthChromSet, dataFile, opts)
    predCount, totalGenes = depthPredictions(commonChromSet, depthDict, jctDict, model, opts)
else :
    # SAM-based predictions:
    if opts.verbose :
        sys.stderr.write('Initializing SAM input from %s\n' % dataFile)
    depth_chrom_set = set([c.lower() for c in getSamSequences(dataFile).keys()])

    # Find chromosomes common to both gene models and alignment records:
    commonChromSet = getCommonChromosomes(geneChromSet, depth_chrom_set, dataFile, opts)

    totalGenes = predCount = 0
    sys.stderr.write('Generating predictions:\n')
    for c in sorted(commonChromSet):
        sys.stderr.write('  starting chromosome %s\n' % c)
        sys.stderr.write('    loading alignment records from %s\n' % dataFile)
        depthDict, jctDict = getSamReadData(
            dataFile,
            minjct=opts.minjct,
            minanchor=opts.anchor,
            chromosomes=[c],
        )

        sys.stderr.write('    predicting splice graphs ')
        geneList = model.getGeneRecords(c)
        totalGenes += len(geneList)
        geneList.sort()

        depths = depthDict[c] if c in depthDict else []
        jcts = jctDict[c] if c in jctDict else []
        outDir = os.path.join(opts.outdir, c)
        if not os.path.isdir(outDir) :
            os.makedirs(outDir)
        predCount += predictGraphs(geneList, depths, jcts, outDir, depthThreshold=opts.mindepth, novelOnly=opts.novel)

sys.stderr.write('Finished: splice graphs were modified for %s / %s genes.\n' % \
        (commaFormat(predCount), commaFormat(totalGenes)))
