#!/usr/bin/env python3


import argparse, sys, os, subprocess


def parse_args():
    parser = argparse.ArgumentParser(prog='prepare-directory.py')
    parser.add_argument('-o', '--outdir', default='samples')
    parser.add_argument('infiles', nargs='+')
    args = parser.parse_args()
    return args


def prep_directory(args):

    if (len(args.infiles) % 2) != 0:
        print('Error: odd number of infiles (%d)' % len(args.infiles))
        sys.exit(1)

    infiles = sorted(args.infiles)
    infiles_1, infiles_2 = infiles[::2], infiles[1::2]

    counter = 0

    for (infile_1, infile_2) in zip(infiles_1, infiles_2):
        prefix_1 = os.path.basename(infile_1).split('_', 1)[0]
        prefix_2 = os.path.basename(infile_2).split('_', 1)[0]
        if prefix_1 != prefix_2:
            print('Error: unmatch pair [%s] - [%s]' % (prefix_1, prefix_2))
            sys.exit(1)

        cmd1 = 'ln -s ../../../%s %s/%s/reads/raw.R1.fastq.gz' % (infile_1, args.outdir, prefix_1)
        cmd2 = 'ln -s ../../../%s %s/%s/reads/raw.R2.fastq.gz' % (infile_2, args.outdir, prefix_2)

        os.makedirs('%s/%s/reads' % (args.outdir, prefix_1))
        os.makedirs('%s/%s/sanger' % (args.outdir,prefix_1))
        print('Preparing for [%s]' % prefix_1)
        subprocess.run(cmd1, shell=True)
        subprocess.run(cmd2, shell=True)
        counter += 1

    print("Preparing directories for %d samples." % counter)
        

def main():
    args = parse_args()
    prep_directory(args)


if __name__ == '__main__':
    main()
