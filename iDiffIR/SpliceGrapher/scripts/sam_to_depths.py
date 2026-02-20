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
from iDiffIR.SpliceGrapher.shared.config import *
from iDiffIR.SpliceGrapher.shared.utils  import *
from iDiffIR.SpliceGrapher.formats.alignment_io   import *
from iDiffIR.SpliceGrapher.shared.ShortRead import *
import argparse
import os,sys

USAGE = """%(prog)s SAM-file [options]

Converts a SAM/BAM file into a SpliceGrapher depth file."""

# Establish command-line options:
parser = argparse.ArgumentParser(usage=USAGE)
parser.add_argument('-o', dest='output',  default=None,  help='Output file [default: same file with .depths extension]')
parser.add_argument('-v', dest='verbose', default=False, help='Verbose mode [default: %(default)s]', action='store_true')
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

samFile = args[0]
validateFile(samFile)

if not opts.output :
    pfx,ignore  = os.path.splitext(samFile)
    opts.output = '%s.depths' % pfx

if os.path.isfile(opts.output) :
    sys.stderr.write('%s will be over-written after SAM data has been loaded\n' % opts.output)

outStream = open(opts.output, 'w')
if opts.verbose :
    sys.stderr.write('converting alignment records to depths\n')
depthDict,jctDict = getSamReadData(samFile, verbose=opts.verbose)
if opts.verbose :
    sys.stderr.write('writing to %s\n' % opts.output)
writeDepths(outStream, depthDict=depthDict, jctDict=jctDict)

if opts.verbose : sys.stderr.write('Finished.\n')
