
from seqpy import cout, cerr
from seqpy.core import bioio
from seqpy.cmds import arg_parser

import pandas as pd


def prepare_submission(args):

    out_metadata = args.outprefix + '.csv'
    out_fasta = args.outprefix + '.fas'

    # open metadata file
    if args.metafile.lower().endswith('.csv'):
        separator = ','
    elif args.metafile.lowe().endswith('.tsv'):
        separator = '\t'
    metadata_df = pd.read_table(args.metafile, sep=separator)

    # make sure sequence name is a string (in case the the column is automatically
    # converted to number)
    metadata_df['fn'] = metadata_df['fn'].astype('str')

    # open depth file
    depth_df = pd.read_table(args.depths, sep='\t', header=None)

    # change the first column to string
    depth_df[0] = depth_df[0].astype('str')
    depth_df.set_index(0, inplace=True)

    # open sequence file
    mseq= bioio.load( args.infile )
    mseq_keys = {}
    for i in range(len(mseq)):
        mseq_keys[ mseq[i].label ] = i

    # iterate over metadata rows
    unused = []
    #import IPython; IPython.embed()
    for (i, r) in metadata_df.iterrows():

        if i == 0: continue

        seq_name = r['fn']
        if seq_name not in mseq_keys:
            unused.append( i )
            continue

        # set coverage
        # import IPython; IPython.embed()
        metadata_df.at[i, 'covv_coverage'] = depth_df.loc[seq_name, 1]
        metadata_df.at[i, 'fn'] = out_fasta
        metadata_df.at[i, 'covv_assembly_method'] = 'custom minimap2 + iVar pipeline'

        # set sequence name
        idx = mseq_keys[seq_name]
        mseq[idx].label = r['covv_virus_name']

    # remove unused metadata
    metadata_df.drop( unused, inplace=True )

    # write to new fasta & metadata file
    metadata_df.to_csv(out_metadata, sep=',', index=False)
    bioio.save(mseq, out_fasta)


def init_argparser():

    p = arg_parser('prepare GISAID submission files')

    p.add_argument('--metafile', required=True)
    p.add_argument('--depths', required=True)
    p.add_argument('--outprefix', required=True)
    p.add_argument('infile')

    return p


def main( args ):
    prepare_submission(args)
