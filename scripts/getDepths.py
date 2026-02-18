#! /usr/bin/env python

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

import sys, os, numpy, pysam, gzip
from subprocess import Popen, PIPE
from argparse import ArgumentParser
from iDiffIR.SpliceGrapher.shared.utils      import *


MATCH    = 0   #M
INSERT   = 1   #I
DELETE   = 2   #D
GAP      = 3   #N
SOFT_CLIP= 4   #4
HARD_CLIP= 5   #H
PAD      = 6   #P
EQUAL    = 7   #=
DIFF     = 8   #X

def parseArgs():
    parser = ArgumentParser(description='Get chromosomal read depths and junctions')

    parser.add_argument('-o', '--output-dir',dest='outdir', action='store',
                        default='iDiffIR_counts', help="output file directory name")

    parser.add_argument('-v', '--verbose',dest='verbose', action='store_true',
                        default=False, help="verbose output")

    parser.add_argument( 'bamfile_in', action='store', type=str, default=None,
                         help='Name of sorted, indexed BAM file')


    args = parser.parse_args()

    # make output directory
    if not os.path.exists( args.outdir ):
        os.makedirs( args.outdir )
    return args

def processRead( read, depths, junctions):
    rpos = 0
    pos = read.pos

    for i,rec in enumerate(read.cigar):
        t,l = rec
        if t == MATCH:
            depths[pos:(pos+l)] += 1
            pos = pos+ l
            rpos = rpos + l

        #insertion
        elif t == INSERT:
            rpos += l
        #deletion
        elif t == DELETE:
            depths[pos] += 1
            pos += l
        # junction
        elif t == GAP:
            jct = (pos, pos+l)
            junctions[jct] = junctions.get(jct, 0) + 1
            pos += l

def writeDepths( depths, chromosome, out_dir, chrom_len, verbose ):
    """
    write depths to file
    """
    if verbose:
        sys.stderr.write('Compressing and writing depths (%s positions)\n' % (commaFormat(chrom_len)))
    indicator = ProgressIndicator(10000000)
    with open( os.path.join(out_dir, '%s.cnt.gz' % (chromosome.lower())), 'wb') as fout:
        pipe = Popen('gzip', stdin=PIPE, stdout=fout)
        for d in depths:
            pipe.stdin.write('%d ' % d )
            indicator.update()
        indicator.finish()
        pipe.communicate()


def main( ):
    args = parseArgs()
    bamfile = pysam.Samfile(args.bamfile_in, 'rb')
    lmap = dict(zip(bamfile.references, bamfile.lengths))
    tot_nonzero   = 0
    tot_positions = 0
    tot_junctions = 0
    read = next(bamfile)
    chromosome = bamfile.getrname(read.tid)
    if args.verbose:
        sys.stderr.write('-'*70+'\n')
        sys.stderr.write('Processing chromosome %s alignments: \n' % ( chromosome ))
    indicator = ProgressIndicator(1000000, verbose=args.verbose)
    depths = numpy.zeros(lmap[chromosome], int)
    junctions = { }
    nonzero = 0
    processRead( read, depths, junctions )
    indicator.update()
    for read in bamfile:
        rchrom = bamfile.getrname(read.tid)
        # found new chromosome
        if rchrom != chromosome:
            if args.verbose:
                sys.stderr.write( '%s reads' % commaFormat(indicator.count()))
            indicator.finish()

            writeDepths( depths, chromosome, args.outdir, lmap[chromosome], args.verbose)
            nonzero = sum(depths > 0)

            if args.verbose:
                sys.stderr.write('Writing junctions\n')
            with gzip.open( os.path.join(args.outdir, '%s.jct.gz' % (chromosome.lower())), 'wb') as fout:
                for key in junctions:
                    fout.write('%d\t%d\t%d\n' % (key[0], key[1], junctions[key] ) )

            tot_nonzero   += nonzero
            tot_junctions += len(junctions)
            if args.verbose:
                sys.stderr.write('Coverage: %s / %s (%0.2f%%) positions have non-zero depth\n' % \
                                 (commaFormat(nonzero), commaFormat(lmap[chromosome]), float(nonzero)/lmap[chromosome] * 100))
                sys.stderr.write('Junctions: %s\n' % commaFormat(len(junctions)) )
                sys.stderr.write('+'*70+'\n')

            chromosome = bamfile.getrname(read.tid)
            if args.verbose:
                sys.stderr.write('-'*70+'\n')
                sys.stderr.write('Processing chromosome %s alignments: \n' % ( chromosome ))
            indicator = ProgressIndicator(1000000, verbose=args.verbose)
            depths = numpy.zeros( lmap[chromosome], int)
            junctions = { }
        processRead(read, depths, junctions)
        indicator.update()
    if args.verbose:
        sys.stderr.write( '%s reads' % commaFormat(indicator.count()))
    indicator.finish()

    writeDepths( depths, chromosome, args.outdir, lmap[chromosome], args.verbose)
    nonzero = sum(depths > 0)

    if args.verbose:
        sys.stderr.write('Writing junctions\n')
    with gzip.open( os.path.join(args.outdir, '%s.jct.gz' % (chromosome.lower())), 'wb') as fout:
        for key in junctions:
            fout.write('%d\t%d\t%d\n' % (key[0], key[1], junctions[key] ) )

    if args.verbose:
        sys.stderr.write('Coverage: %s / %s (%0.2f%%) positions have non-zero depth\n' % \
                         (commaFormat(nonzero),
                          commaFormat(lmap[chromosome]),
                          float(nonzero)/lmap[chromosome] * 100))
        sys.stderr.write('Junctions: %s\n' % commaFormat(len(junctions)) )
        sys.stderr.write('+'*70+'\n')

    tot_nonzero   += nonzero
    tot_positions = sum(lmap.values())
    tot_junctions += len(junctions)

    if args.verbose:
        sys.stderr.write('*'*70 + '\n')
        sys.stderr.write('Total Coverage: %d / %d (%0.2f%%) positions have non-zero depth\n' % \
                         (tot_nonzero, tot_positions, float(tot_nonzero)/tot_positions* 100))
        sys.stderr.write('Total Junctions: %d\n' % tot_junctions)
        sys.stderr.write('*'*70 + '\n')


if __name__ == '__main__':
    main()
