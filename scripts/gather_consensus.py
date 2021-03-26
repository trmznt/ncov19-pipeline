#!/usr/bin/env spcli

# (c) Hidayat Trimarsanto <anto@eijkman.go.id>

from seqpy import cout, cerr, cexit
from seqpy.cmds import arg_parser
from seqpy.core.bioio import load, save, multisequence, biosequence
from seqpy.core.funcs import profiles


def init_argparser(p=None):

    if p is None:
        p = arg_parser('gather_consensus.py - gather consensus')

    p.add_argument('--add', default=None)
    p.add_argument('-o', '--outfile', default='consensus.fas')
    p.add_argument('-i', '--infile', default='/cons/trimmed.fa')
    p.add_argument('indirs', nargs='+')

    return p

def main( args ):

    gather_consensus( args )


def gather_consensus( args ):

    # open input file
    cons = multisequence()

    if args.add:
        seqs = load(args.add)
        cons.extend( seqs )

    for indir in args.indirs:
        seqs = load(indir + args.infile)

        cons.append( biosequence(indir, seqs[0].seq))

    save( cons, args.outfile )
