#! /usr/bin/env python

"""Create a release package for distribution."""

import os
import shutil
import subprocess
import sys
from argparse import ArgumentParser


def copy_path(command: str, source: str, target: str) -> None:
    """Copy one file or directory to the given release target directory."""
    subprocess.call(command % (source, target), shell=True)


def main() -> None:
    """Build the release directory layout and copy expected project artifacts."""
    parser = ArgumentParser(description='Build iDiffIR release')
    parser.add_argument('version', type=str, help='version number')
    parser.add_argument('-v', '--verbose', dest='verbose', 
                        action='store_true', default='False',
                      help='Print verbose output')

    parser.add_argument('-f', '--force-overwrite', dest='force', 
                        action='store_true', default='False',
                      help='Force overwrite the directory if it exists')
    
    nspace = parser.parse_args()

    nspace.outdir = f'iDiffIR_v{nspace.version}/'

    if os.path.exists(nspace.outdir):
        if not nspace.force:
            sys.exit(f'{nspace.outdir} exists, use -f to overwrite\n')
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
    copy_path(cmd, 'LICENSE', nspace.outdir)
    # copy README
    copy_path(cmd, 'README.md', nspace.outdir)
    # copy contributing guide
    copy_path(cmd, 'CONTRIBUTING.md', nspace.outdir)
    # copy code of conduct
    copy_path(cmd, 'CODE_OF_CONDUCT.md', nspace.outdir)
    # copy CONTRIBUTORS
    copy_path(cmd, 'CONTRIBUTORS.md', nspace.outdir)

    # copy modern build metadata
    copy_path(cmd, 'pyproject.toml', nspace.outdir)
    # copy lockfile for reproducible environments
    copy_path(cmd, 'uv.lock', nspace.outdir)

    # copy scripts dir
    copy_path(cmd, 'scripts', nspace.outdir)

    # copy iDiffIR dir
    copy_path(cmd, 'iDiffIR', nspace.outdir)
    

    # copy documentation build
    if os.path.exists('doc/_build/html') and \
       os.path.exists('doc/_build/latex/iDiffIR.pdf'):
        copy_path(cmd, 'doc/_build/html', os.path.join(nspace.outdir, 'doc'))
        copy_path(
            cmd,
            'doc/_build/latex/iDiffIR.pdf',
            os.path.join(nspace.outdir, 'doc', 'pdf'),
        )
    else:
        sys.stderr.write('Missing documentation directories\n')


if __name__ == "__main__":
    main()
