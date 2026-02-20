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
from iDiffIR.SpliceGrapher.formats.loader import *
import argparse
import sys

USAGE = """%(prog)s gene-model-file [options]

Converts a gene model file in GTF or GFF3 format into an output
file in GTF or GFF3 format."""

# Establish command-line options:
parser = argparse.ArgumentParser(usage=USAGE)
parser.add_argument('--gff3', dest='gff3',    default=None,  help='Output file for GFF3 format [default: stdout]')
parser.add_argument('--gtf',  dest='gtf',     default=None,  help='Output file for GTF format [default: stdout]')
parser.add_argument('-v',     dest='verbose', default=False, help='Verbose mode [default: %(default)s]', action='store_true')
def _parse_opts_and_args(parser, argv):
    parser.add_argument('args', nargs='*')
    opts = parser.parse_args(argv)
    args = opts.args
    delattr(opts, 'args')
    return opts, args

opts, args = _parse_opts_and_args(parser, sys.argv[1:])
if len(args) != 1 :
    parser.print_help()
    sys.exit(1)

if not (opts.gff3 or opts.gtf) :
    parser.print_help()
    sys.stderr.write('\nNo output file specified (--gff3 or --gtf); exiting.\n')
    sys.exit(1)

inputFile = args[0]
validateFile(inputFile)

if not opts.verbose : sys.stderr.write('Loading gene models from %s\n' % inputFile)
model = loadGeneModels(inputFile, verbose=opts.verbose, alltypes=True)

if opts.gff3 :
    sys.stderr.write('Writing GFF3 output to %s\n' % opts.gff3)
    model.writeGFF(opts.gff3, verbose=opts.verbose)

if opts.gtf :
    sys.stderr.write('Writing GTF output to %s\n' % opts.gtf)
    model.writeGTF(opts.gtf, verbose=opts.verbose)
