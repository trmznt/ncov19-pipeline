#!/usr/bin/env python3

import argparse
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(prog='combine-report.py')
    p.add_argument('-o', '--outfile', default='full-reports.csv')
    p.add_argument('--outpass', default='qc-pass.tsv')
    p.add_argument('--outfailed', default='qc-failed.tsv')
    p.add_argument('-m', '--metafile')
    p.add_argument('-s', '--statfile', default='stats.tsv')
    p.add_argument('-l', '--lineagefile', default='lineage_report.csv')
    return p.parse_args()


def combine_report(args):

    df_metadata = pd.read_table(args.metafile, sep=',')

    df_stat = pd.read_table(args.statfile, sep='\t')
    df = df_metadata.join(df_stat.set_index('SAMPLE'), on='fn')

    df_lineage = pd.read_table(args.lineagefile, sep=',')
    df = df.join(df_lineage.set_index('taxon'), on='fn')

    df.to_csv(args.outfile, index=False, sep=',')

    # filtering criteria
    # must have less than 1000 N and more than 28000 bp
    # either have ambiguity_score > 0.8 or being predicted by PANGO
    pass_filter = (df['N_BASE'] < 1000) & (df['LENGTH'] > 28000) & ((df['ambiguity_score'] > 0.8) | df['version'].str.startswith('PANGO'))
    df_pass = df[pass_filter]
    df_pass.to_csv(args.outpass, index=False, sep='\t')
    df_failed = df[~ pass_filter]
    df_failed.to_csv(args.outfailed, index=False, sep='\t')
    print(f'Sequence passed QC: {len(df_pass)}')


def main():
    args = parse_args()
    combine_report(args)


if __name__ == '__main__':
    main()

# EOF
