#!/usr/bin/env python3

import argparse
import sys
import itertools

try:
    from matplotlib import pyplot as plt
except:
    print('ERR: require properly installed matplotlib')
    sys.exit(1)

import numpy as np
import pysam
import yaml

MAX_INSERT = 5000

def init_argparser():
    p = argparse.ArgumentParser("Create depth plot")
    p.add_argument('-t', '--title', default='depth')
    p.add_argument('-d', '--mindepth', type=int, default=3)
    p.add_argument('--maxinsert', type=int, default=MAX_INSERT)
    p.add_argument('--stat_insert', default=False, action='store_true')
    p.add_argument('--log_insert', default=False, action='store_true')
    p.add_argument('--pool', type=int, default=-1)
    p.add_argument('--outgap', default='')
    p.add_argument('--outdepth', default='')
    p.add_argument('-o', '--outfile', default='outplot.png')

    p.add_argument('infile')

    return p


def ranges(i):
    for a, b in itertools.groupby(enumerate(i), lambda pair: pair[1] - pair[0]):
        b = list(b)
        yield b[0][1], b[-1][1]


def depthplot( args ):

    if args.stat_insert:
        fig, axs = plt.subplots(1, 2, figsize=(30,5),
                gridspec_kw={'width_ratios': [5, 1]})
        ax = axs[0]
    else:
        fig, ax = plt.subplots(figsize=(25,5))

    segments = pysam.AlignmentFile(args.infile, 'rb')
    length = segments.lengths[0]
    print('Reference length: %d' % length)

    # variables
    depths = ( np.zeros(length), np.zeros(length), np.zeros(length) )
    missing_regions = ( [], [], [] )
    inserts = ( np.zeros(args.maxinsert), np.zeros(args.maxinsert), np.zeros(args.maxinsert))
    max_insert = 0
    start_reads = {}

    # calculate depths and insert sizes
    for read in segments:
        try:
            amp_id = read.get_tag('XP')
            if amp_id < 0:
                pool = 2
            else:
                pool = amp_id % 2
        except KeyError:
            pool = 2
        depths[pool][read.reference_start : read.reference_end] += 1
        if read.query_name in start_reads:
            if read.is_reverse:
                insert_size = read.reference_end - start_reads[read.query_name]
            else:
                insert_size = start_reads[read.query_name] - read.reference_start
            del start_reads[read.query_name]
            if 0 < insert_size < args.maxinsert:
                max_insert = max(max_insert, insert_size)
                inserts[pool][insert_size] += 1
            else:
                print('WARN: impossible insert size of %d' % insert_size)
        else:
            if read.is_reverse:
                start_reads[read.query_name] = read.reference_end
            else:
                start_reads[read.query_name] = read.reference_start

    all_inserts = sum(inserts)
    for max_insert_threshold in range(max_insert, 0, -1):
        if all_inserts[max_insert_threshold] > 1:
            break
    max_insert = max_insert_threshold

    if args.pool >= 0:
        analyzed_depths = [ depths[args.pool] ]
    else:
        analyzed_depths = depths

    colors = itertools.cycle(['forestgreen', 'royalblue', 'lightcoral'])
    for y in analyzed_depths:

        color = next(colors)
        x = np.arange(len(y))
        ax.fill_between(x, y, facecolor=color, color=color, alpha=0.5)

    ax.set_yscale("log", base=10)
    ax.set_ylim(ymin=0.7)
    ax.set_title(args.title)

    # check gaps
    y = sum(analyzed_depths)

    # check for holes (depth < 1):
    holes = []
    for i in range(len(y)):
        if y[i] < args.mindepth:
            holes.append(i)

    print('Missing region(s):')
    regions = []
    for region in list(ranges(holes)):
        regions.append(list(region))
        print(region)
    if args.outgap:
        yaml.dump({args.title: regions}, open(args.outgap, 'w'))

    insert_sizes = np.arange(max_insert)
    if args.stat_insert:

        colors = itertools.cycle(['forestgreen', 'royalblue', 'lightcoral'])
        if args.pool >= 0:
            analyzed_inserts = [ inserts[args.pool] ]
        else:
            analyzed_inserts = inserts

        for insert_count in analyzed_inserts:

            insert_count = insert_count[:max_insert]

            # plot data
            axs[1].hist(insert_sizes, len(insert_sizes), weights=insert_count, color=next(colors), alpha=0.5)
            if args.log_insert:
                axs[1].set_yscale("log", base=10)
            axs[1].set_title('Insert Size')

    fig.tight_layout()
    fig.savefig(args.outfile)

    # calculate average depths
    if args.outdepth:
        y = y[ y > args.mindepth ]
        depth = y.mean() if len(y) > 0 else 0
        inserts_freqs = sum( inserts )[:len(insert_sizes)]
        mean_insize = np.average(insert_sizes, weights = inserts_freqs)
        cdf = np.cumsum(inserts_freqs[insert_sizes])
        med_insize = np.searchsorted(cdf, cdf[-1] // 2)
        dev = inserts_freqs * (insert_sizes - mean_insize) ** 2
        stddev = np.sqrt( dev.sum() / (inserts_freqs.sum() - 1) )
        with open(args.outdepth, 'w') as fout:
            fout.write('SAMPLE\tDEPTH\tBASES\tMEAN_IS\tMED_IS\tSD_IS\n')
            fout.write(f'{args.title}\t{int(depth)}\t{len(y)}\t{mean_insize:5.1f}\t{med_insize:5.1f}\t{stddev:5.1f}\n')



def main():
    depthplot( init_argparser().parse_args() )


if __name__ == '__main__':
    main()

