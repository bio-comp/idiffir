#! /usr/bin/python
from iDiffIR.SpliceGrapher.shared.config             import *
from iDiffIR.SpliceGrapher.shared.adjust             import *
from iDiffIR.SpliceGrapher.shared.GeneModelConverter import *
from iDiffIR.SpliceGrapher                           import SpliceGraph
from iDiffIR.SpliceGrapher.shared.utils              import *

g = SpliceGraph.getFirstGraph('AT2G04700.gff')
