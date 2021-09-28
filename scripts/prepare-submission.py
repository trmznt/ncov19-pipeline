
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
    cerr(f'[Reading metadata file {args.metafile}]')
    metadata_df = pd.read_table(args.metafile, sep=separator)

    # make sure sequence name is a string (in case the the column is automatically
    # converted to number)
    metadata_df['fn'] = metadata_df['fn'].astype('str')
    metadata_df['covv_assembly_method'] = metadata_df['covv_assembly_method'].astype('str')
    metadata_df.set_index('fn', drop=False, inplace=True )

    #import IPython; IPython.embed()

    # open infile tsv
    cerr(f'[Reading infile {args.infile}]')
    submission_df = pd.read_table(args.infile, sep='\t')
    submission_df['SAMPLE'] = submission_df['SAMPLE'].astype('str')

    # open sequence file
    cerr(f'[Reading sequence file {args.seqfile}]')
    mseq = bioio.load( args.seqfile )
    mseq_keys = {}
    for i in range(len(mseq)):
        mseq_keys[ mseq[i].label ] = i

    # iterate over submission_df
    used = []
    #import IPython; IPython.embed()
    for (i, s) in submission_df.iterrows():

        sample_id = s['SAMPLE']
        r = metadata_df.loc[sample_id]

        if sample_id not in mseq_keys:
            continue

        cerr(f'[Preparing sample {sample_id}]')
        # set coverage
        # import IPython; IPython.embed()
        metadata_df.at[sample_id, 'covv_coverage'] = s['AVGDEPTH']
        metadata_df.at[sample_id, 'fn'] = out_fasta
        metadata_df.at[sample_id, 'covv_seq_technology'] = args.covv_seq_technology
        metadata_df.at[sample_id, 'covv_assembly_method'] = args.covv_assembly_method

        # set sequence name
        idx = mseq_keys[sample_id]
        mseq[idx].label = r['covv_virus_name']
        mseq[idx].seq = mseq[idx].seq.strip(b'-')
        used.append(sample_id)
        cerr(f'[Finish preparing  sample {sample_id}]')

    # remove unused metadata
    metadata_df = metadata_df.loc[ used ]

    # write to new fasta & metadata file
    metadata_df.to_csv(out_metadata, sep=',', index=False)
    bioio.save(mseq, out_fasta)


def init_argparser():

    p = arg_parser('prepare GISAID submission files')

    p.add_argument('--covv_assembly_method', default='custom minimap2 + iVar pipeline')
    p.add_argument('--covv_seq_technology', default='Illumina NextSeq 550')
    p.add_argument('--metafile', required=True)
    p.add_argument('--seqfile', required=True)
    p.add_argument('--outprefix', required=True)
    p.add_argument('infile')

    return p


def main( args ):
    prepare_submission(args)
