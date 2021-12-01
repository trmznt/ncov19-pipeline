#!/usr/bin/env python3

# (c) Hidayat Trimarsanto <anto@eijkman.go.id>

import argparse, sys


def parse_args():

    p = argparse.ArgumentParser('polish_consensus.py - polish consensus sequence')

    p.add_argument('-o', '--outfile', default='consensus.fas')
    p.add_argument('-l', '--label', default='consensus')
    p.add_argument('infile')

    return p.parse_args()


def main():

    args = parse_args()
    polish_consensus( args )


def polish_consensus( args ):

    # open
    lines = open(args.infile).read().split('\n')
    sequence = ''
    label = None
    for line in lines:
        if line.startswith('>'):
            if label is None:
                label = line[1:].strip()
            else:
                return -1, -1
            continue
        sequence += line.strip()

    # remove leading and trailing N
    sequence = sequence.strip('N')

    with open(args.outfile, 'w') as fout:
        fout.write(f'>{args.label}\n')
        fout.write(sequence)
        fout.write('\n')

    sys.stderr.write(f'[Writing polished consensus to {args.outfile}]\n')


if __name__ == '__main__':
    main()

