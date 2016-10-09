#! /usr/bin/env python
import sys, subprocess, os
from argparse import ArgumentParser, ArgumentTypeError

GFF_CMD   = 'make_MISO_IR_GFF.py -m %s -o %s'
INDEX_CMD = 'index_gff --index %s %s'
RUN_CMD   = 'miso --run %s %s --output-dir %s --read-len %d -p %d'
COMP_CMD  = 'compare_miso --compare-samples %s %s %s'
def parseArgs():
    """Parse command line arguments

    Returns
    -------
    a : argparse.ArgumentParser

    """
    parser = ArgumentParser(description='Run MISO for identifying differential intron retention.')
    parser.add_argument('-v', '--verbose',dest='verbose', action='store_true', 
                        default=False, help="verbose output [default is quiet running]")
    parser.add_argument('-n', '--noplot', dest='noplot', action='store_true', 
                        default=False, help="Do not plot figures [default is to make figures]")
    parser.add_argument('-l', '--factorlabel', dest="factorlabels",action='store', nargs=2,
                        type=str, default=['factor1','factor2'], 
                        help="factor labels, example:  -f Mutant Wildtype", metavar='FACTORLABEL')
    parser.add_argument('-o', '--output-dir',dest='outdir', action='store', 
                        default='miso_ir_output', help="output file directory name")
    parser.add_argument('-s', '--shrink_introns', dest='shrink_introns', action='store_true', 
                        default=False, help='shrink introns for depth plots [default is no shrinking]')
    parser.add_argument('-p', '--procs', dest='procs', action='store', default=1, 
                        type=int, help='Number of processing cores to use, [default = 1]')
    parser.add_argument('-r', '--read-length', dest='readlength', action='store', default=100, 
                        type=int, help='Number of processing cores to use, [default = 100]')

    #parser.add_argument('-g', '--graph-dirs', dest='graphDirs', type=fileList, 
    #                    help='colon-separated list of directories to recursively search for SpliceGrapher predictions')
    parser.add_argument('genemodel', type=str,
                        help="gene model file: NAME.gtf[.gz] | NAME.gff[.gz]")
    parser.add_argument('factor1bamfile', type=str,
                        help="bamfile for 1st factor")
    parser.add_argument('factor2bamfile', type=str,
                        help="bamfile for 2nd factor")
    
    return parser.parse_args()

def main():
    args = parseArgs()
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    gffloc = os.path.join(args.outdir, 'miso_ir.gff')
    # Make MISO GFF
    if args.verbose:
        sys.stderr.write('Generating MISO GFF of IR events\n')
    cmd = GFF_CMD % (args.genemodel, gffloc)
    if args.verbose:
        cmd += ' -v'
    status = subprocess.call(cmd, shell=1)
    # Index MISO GFF
    indexloc = os.path.join(args.outdir, 'ir_index')
    cmd = INDEX_CMD % (gffloc, indexloc)
    status = subprocess.call(cmd, shell=1)
    # Run miso on sample 1
    factor1loc = os.path.join(args.outdir, 'factor1')
    cmd = RUN_CMD % ( indexloc, args.factor1bamfile, factor1loc, args.readlength, args.procs)
    status = subprocess.call(cmd, shell=1)
    # Run miso on sample 2
    factor2loc = os.path.join(args.outdir, 'factor2')
    cmd = RUN_CMD % ( indexloc, args.factor2bamfile, factor2loc, args.readlength, args.procs)
    status = subprocess.call(cmd, shell=1)
    # Run differential comparisons
    comploc = os.path.join(args.outdir, 'comparisons')
    cmd = COMP_CMD % (factor1loc, factor2loc, comploc)
if __name__ == '__main__':
    main()
