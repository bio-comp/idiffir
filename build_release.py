from argparse import ArgumentParser, ArgumentTypeError
import os, sys, subprocess

def main():
    parser = ArgumentParser(description='Build iDiffIR release')
    parser.add_argument('version', type=str, help='version number')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', default='False',
                      help='Print verbose output')

    parser.add_argument('-f', '--force-overwrite', dest='force', action='store_true', default='False',
                      help='Force overwrite the directory if it exists')
    
    nspace = parser.parse_args()

    nspace.outdir = '../v%s' % nspace.version

    if os.path.exists(nspace.outdir):
        if not nspace.force:
            sys.exit('%s exists, use -f to overwrite\n' % nspace.outdir)
    else:
        os.makedirs(nspace.outdir)

    cmd = 'cp -R %s %s' 

    # copy license
    status_sort = subprocess.call(cmd % ('LICENSE', nspace.outdir),shell=True)
    # copy README
    status_sort = subprocess.call(cmd % ('README', nspace.outdir),shell=True)
    


if __name__ == "__main__":
    main()
