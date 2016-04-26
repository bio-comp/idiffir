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
iDiffIR plotting
"""
import matplotlib
matplotlib.use('agg')
from iDiffIR.Plot import *

def parseArgs():
    parser = ArgumentParser(description='Identify differentially expressed introns.')
    parser.add_argument('-v', '--verbose',dest='verbose', action='store_true', 
                        default=False, help="verbose output [default is quiet running]")

    parser.add_argument('-l', '--factorlabel', dest="factorlabels",action='store', nargs=2,
                        type=str, default=['factor1','factor2'], 
                        help="factor labels, example:  -f Mutant Wildtype", metavar='FACTORLABEL')
    parser.add_argument('-o', '--output-dir',dest='outdir', action='store', 
                        default='iDiffIR_output', help="output file directory name")
    parser.add_argument('-s', '--shrink_introns', dest='shrink_introns', action='store_true', 
                        default=False, help='shrink introns for depth plots [default is no shrinking]')
    parser.add_argument('-g', '--graph-dirs', dest='graphDirs', type=fileList, 
                        help='colon-separated list of directories to recursively search for SpliceGrapher predictions')
    parser.add_argument('genemodel', type=str,
                        help="gene model file: NAME.gtf[.gz] | NAME.gff[.gz]")
    parser.add_argument('factor1Dirs', type=fileList,
                        help="colon-separated list of directories: PATH-TO-REPLICATE_1[:PATH-TO-REPLICATE_2,...]")
    parser.add_argument('factor2Dirs', type=fileList,
                        help="colon-separated list of directories: PATH-TO-REPLICATE_1[:PATH-TO-REPLICATE_2,...]")

    args = parser.parse_args()
    if not validateArgs( args ):
        raise Exception("Argument Errors: check arguments and usage!")
    return args
