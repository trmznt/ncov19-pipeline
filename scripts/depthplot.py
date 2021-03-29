#!/usr/bin/env python3

import argparse, sys, gzip

try:
    from matplotlib import pyplot as plt
except:
    print('ERR: require properly installed matplotlib')
    sys.exit(1)

import numpy as np

def init_argparser():
    p = argparse.ArgumentParser("Create depth plot")
    p.add_argument('-t', '--title', default='depth')
    p.add_argument('-d', '--mindepth', type=int, default=3)
    p.add_argument('--outgap', default='')
    p.add_argument('--statfile', default='')
    p.add_argument('-o', '--outfile', default='outplot.png')

    p.add_argument('infiles', nargs='+')

    return p

def main( args ):

    depthplot( args )


import itertools
import yaml

colors = itertools.cycle(['forestgreen', 'royalblue', 'lightcoral'])


def ranges(i):
    for a, b in itertools.groupby(enumerate(i), lambda pair: pair[1] - pair[0]):
        b = list(b)
        yield b[0][1], b[-1][1]


def depthplot( args ):


    if args.statfile:
        fig, axs = plt.subplots(1, 2, figsize=(30,5),
                gridspec_kw={'width_ratios': [5, 1]})
        ax = axs[0]
    else:
        fig, ax = plt.subplots(figsize=(25,5))


    for depthfile in args.infiles:

        color = next(colors)

        # read data
        depth_list = []
        with gzip.open(depthfile, 'rt') as fin:
            for line in fin:
                tokens = line.split()
                depth_list.append( (int(tokens[1]), int(tokens[2])) )

        # plot data
        length = depth_list[-1][0] + 1
        x = np.arange(length)
        y = np.zeros(length)

        for (pos, depth) in depth_list:
            y[pos] = depth

        # check for holes (depth < 1):
        holes = []
        for i in range(length):
            if y[i] < args.mindepth:
                holes.append(i)

        print('Missing region(s):')
        regions = []
        for region in list(ranges(holes)):
            regions.append(list(region))
            print(region)
        if args.outgap:
            yaml.dump({args.title: regions}, open(args.outgap, 'w'))

        ax.fill_between(x, y, facecolor=color, color=color, alpha=0.5)

    ax.set_yscale("log", base=10)
    ax.set_ylim(ymin=0.7)
    ax.set_title(args.title)

    if args.statfile:

        # read data
        insert_size = []
        insert_count = []
        with gzip.open(args.statfile, 'rt') as fin:
            for line in fin:
                if not line.startswith('IS'):
                    continue
                tokens = line.split()
                insert_size.append( int(tokens[1]) )
                insert_count.append( int(tokens[3]) )

        # plot data
        axs[1].hist(insert_size, len(insert_size), weights=insert_count, color='plum')
        axs[1].set_yscale("log", base=10)

    fig.tight_layout()
    fig.savefig(args.outfile)

def main():
    depthplot( init_argparser().parse_args() )

if __name__ == '__main__':
    main()

