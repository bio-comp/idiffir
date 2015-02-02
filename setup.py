#!/usr/bin/env python
"""
Distribution script for DinrSR, a program to identify
differntial intron retention from short read data.

\author Mike Hamilton
"""
from distutils.core import setup
import sys
SCRIPTS  = ['scripts/analyze_introns.py',
            'scripts/getDepths.py',
            'scripts/convertSam.py',
            'scripts/simulate_IR.py',
            'scripts/make_MISO_IR_GFF.py',
            'scripts/make_MISO_AS_GFF.py'
            ]

PACKAGES = [ 'iDiffIR'
             ]


REQUIRES = ['SpliceGrapher',
            'scipy',
            'numpy',
            'matplotlib',
            'pysam'
            ]

# check for dependencies if we're building or installing
installable = True
if 'build' in sys.argv or 'install' in sys.argv:
    sys.stderr.write('Checking for dependencies...\n')
    for package in REQUIRES:
        try:
            p = __import__( package )
            sys.stderr.write('package "%s" found (version %s)\n' % (package, p.__version__))
        except ImportError:
            installable = False
            sys.stderr.write('\tRequired package `%s` not found\n' % package )
    if not installable:
        sys.exit('Missing dependencies--terminating installation')

setup(name='iDiffIR',
      version='0.0.1',
      description='Identifying differential intron retention from RNA-seq',
      author='Michael Hamilton',
      author_email='hamiltom@cs.colostate.edu',
      url='combi.cs.colostate.edu/idiffir',
      packages=PACKAGES,
      requires=REQUIRES,
      scripts=SCRIPTS,
      license='GNU General Public License'
     )
