#!/usr/bin/env spcli

# (c) Hidayat Trimarsanto <anto@eijkman.go.id>

from seqpy import cout, cerr, cexit
from seqpy.cmds import arg_parser
from seqpy.core.bioio import load, save, multisequence, biosequence
from seqpy.core.funcs import profiles
import os


def init_argparser(p=None):

    if p is None:
        p = arg_parser('gather_consensus.py - gather consensus')

    p.add_argument('--add', default=None)
    p.add_argument('-o', '--outdir', default='')
    p.add_argument('-c', '--consfile', default='cons/consensus.fa')
    p.add_argument('-s', '--statfile', default='logs/stats.tsv')
    p.add_argument('indir')

    return p

def main( args ):

    gather_consensus( args )


def gather_consensus( args ):

    # set output directory
    args.outdir = args.indir + '-results' if not args.outdir else args.outdir

    # open input file
    cons = multisequence()
    header = None
    stat_lines = []

    if args.add:
        seqs = load(args.add)
        cons.extend( seqs )

    for indir in sorted(os.listdir(args.indir)):

        seqpath = os.path.join(args.indir, indir, args.consfile)
        print(args.indir, indir, args.consfile, seqpath)
        try:
            seqs = load(seqpath)
        except FileNotFoundError:
            cerr('[WARN: no such file: %s]' % (seqpath) )
            continue

        cons.append( biosequence(indir, seqs[0].seq))

        statpath = os.path.join(args.indir, indir, args.statfile)
        with open(statpath) as fin:
            lines = fin.read().split('\n')
            if not header:
                header = lines[0].strip()
            stat_lines.append( lines[1].strip() )

    try:
        os.mkdir(args.outdir)
    except:
        pass

    save( cons, os.path.join(args.outdir, 'consensus.fas' ) )
    with open( os.path.join(args.outdir, 'stats.tsv'), 'w') as fout:
        fout.write(header)
        fout.write('\n')
        fout.write('\n'.join(stat_lines))

    cerr(f'[Writing results to directory {args.outdir}]')



