#! /usr/bin/env python

from argparse import ArgumentParser, ArgumentTypeError
import os, sys, subprocess, shutil
"""
Create a release package for distribution
"""

def main():
    parser = ArgumentParser(description='Build iDiffIR release')
    parser.add_argument('version', type=str, help='version number')
    parser.add_argument('-v', '--verbose', dest='verbose', 
                        action='store_true', default='False',
                      help='Print verbose output')

    parser.add_argument('-f', '--force-overwrite', dest='force', 
                        action='store_true', default='False',
                      help='Force overwrite the directory if it exists')
    
    nspace = parser.parse_args()

    nspace.outdir = 'iDiffIR_v%s/' % nspace.version

    if os.path.exists(nspace.outdir):
        if not nspace.force:
            sys.exit('%s exists, use -f to overwrite\n' % nspace.outdir)
        else:
            shutil.rmtree(nspace.outdir)
            os.makedirs(nspace.outdir)
            os.makedirs(os.path.join(nspace.outdir, 'doc'))
            os.makedirs(os.path.join(nspace.outdir, 'doc', 'pdf'))
    else:
        os.makedirs(nspace.outdir)
        os.makedirs(os.path.join(nspace.outdir, 'doc', 'pdf'))
    cmd = 'cp -R %s %s' 

    # copy license
    status_copy = subprocess.call(cmd % ('LICENSE', nspace.outdir),
                                  shell=True)
    # copy README
    status_copy = subprocess.call(cmd % ('README.md', nspace.outdir),
                                  shell=True)
    # copy contributing guide
    status_copy = subprocess.call(cmd % ('CONTRIBUTING.md', nspace.outdir),
                                  shell=True)
    # copy code of conduct
    status_copy = subprocess.call(cmd % ('CODE_OF_CONDUCT.md', nspace.outdir),
                                  shell=True)
    # copy CONTRIBUTORS
    status_copy = subprocess.call(cmd % ('CONTRIBUTORS.md', nspace.outdir),
                                  shell=True)

    # copy setup.py
    status_copy = subprocess.call(cmd % ('setup.py', nspace.outdir),
                                  shell=True)

    # copy scripts dir
    status_copy = subprocess.call(cmd % ('scripts', nspace.outdir),
                                  shell=True)

    # copy iDiffIR dir
    status_copy = subprocess.call(cmd % ('iDiffIR', nspace.outdir),
                                  shell=True)
    

    # copy documentation build
    if os.path.exists('doc/_build/html') and \
       os.path.exists('doc/_build/latex/iDiffIR.pdf'):
        status_copy = subprocess.call(cmd % ('doc/_build/html', 
                                             os.path.join(nspace.outdir, 
                                                          'doc')),
                                      shell=True)
        status_copy = subprocess.call(cmd % ('doc/_build/latex/iDiffIR.pdf', 
                                             os.path.join(nspace.outdir, 
                                                          'doc', 'pdf')),
                                      shell=True)
    else:
        sys.stderr.write('Missing documentation directories\n')


if __name__ == "__main__":
    main()
