#!/usr/bin/env python3

import sys

def parse_optical_dedup():

    lines = open('logs/optical_dedup.log').read().split('\n')
    for line in lines:
        if line.startswith('Reads In:'):
            reads_in = line.split()[-1]
            continue
        if line.startswith('Reads Out:'):
            reads_out = line.split()[-1]
            break
    return int(reads_in), int(reads_out)

def parse_adapter_trimming():

    lines = open('logs/adapter_trimming.log').read().split('\n')
    for line in lines:
        if line.startswith('Pairs written'):
            pairs_out = line.split()[-2].replace(',','')
            break
    return int(pairs_out)*2

def parse_unique_pairs():

    lines = open('logs/unique_pairs.log').read().split('\n')
    for line in lines:
        if line.startswith('Reads out:'):
            reads_out = line.split()[-1]
            break
    return int(reads_out)

def parse_markdup():

    lines = open('logs/markdup.log').read().split('\n')
    for line in lines:
        if line.startswith('WRITTEN:'):
            reads_out = line.split()[-1]
            break
    return int(reads_out)

def parse_primal_remover():

    reads_out = 0
    lines = open('logs/primal_remover.log').read().split('\n')
    for line in lines:
        if line.startswith('Trimmed pairs:'):
            reads_out += int(line.split()[-1])
            continue
        if line.startswith('Untrimmed pairs:'):
            reads_out += int(line.split()[-1])
            continue
    return reads_out*2

def parse_depth_counter():

    lines = open('maps/depth-counter.txt').read().split('\n')
    t = lines[1].split()
    depth, bases = int(t[1]), int(t[2]),
    mean_is, med_is, sd_is = round(float(t[3])), round(float(t[4])), round(float(t[5]))
    return depth, bases, mean_is, med_is, sd_is

def parse_consensus_fasta():

    lines = open('cons/consensus.fas').read().split('\n')
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

    return len(sequence), sequence.count('N') + sequence.count('n')

def parse_variants():

    lines = open('cons/variants.tsv').read().split('\n')

    snvs = {}
    point_muts = 0
    inframe_gaps = 0
    ooframe_gaps = 0

    for line in lines[1:]:
        tokens = line.split('\t')
        if len(tokens) < 4:
            continue
        key = (tokens[0], tokens[1])
        if key in snvs:
            continue
        snvs[key] = tokens
        if len(tokens[2]) == 1 and len(tokens[3]) == 1:
            point_muts += 1
            continue
        if len(tokens[2]) == 1 and len(tokens[3]) > 1:
            # deletion
            if tokens[3].count('N') % 3 == 0:
                inframe_gaps += 1
            else:
                ooframe_gaps += 1
    return point_muts, inframe_gaps, ooframe_gaps


def main():

    sample_code = sys.argv[1]
    #report_file = sys.argv[2]

    initial_reads, optical_dedup_reads = parse_optical_dedup()
    adapter_trimmed_reads = parse_adapter_trimming()
    properly_mapped_reads = parse_unique_pairs()

    try:
        pcr_dedup_reads = parse_markdup()
    except:
        # pcr dedup is not performed
        pcr_dedup_reads = properly_mapped_reads

    try:
        primer_trimmed_reads = parse_primal_remover()
    except FileNotFoundError:
        # primer trimming is not performed
        primer_trimmed_reads = pcr_dedup_reads

    depth, bases, mean_is, med_is, sd_is = parse_depth_counter()
    seqlength, n_bases = parse_consensus_fasta()
    point_mutations, inframe_gaps, ooframe_gaps = parse_variants()


    # writing report

    #outfile = open(report_file, 'w')
    outfile = sys.stdout
    headers = [ 'SAMPLE', 'RAW', 'OP_DEDUP', 'OP_DEDUP_R', 'ADAPTER', 'ADAPTER_R',
                'PROP_PAIR', 'PROP_PAIR_R', 'PCR_DEDUP', 'PCR_DEDUP_R', 'PRIMAL', 'PRIMAL_R',
                'MEAN_IS', 'MED_IS', 'SD_IS', 'AVGDEPTH', 'BASE', 'LENGTH', 'N_BASE',
                'POINTMUT', 'INFRAME', 'OOFRAME' ]
    outfile.write('%s\n' % '\t'.join( headers ))
    outfile.write(  f'{sample_code}\t'
                    f'{initial_reads}\t'
                    f'{optical_dedup_reads}\t{optical_dedup_reads/initial_reads:5.3f}\t'
                    f'{adapter_trimmed_reads}\t{adapter_trimmed_reads/initial_reads:5.3f}\t'
                    f'{properly_mapped_reads}\t{properly_mapped_reads/initial_reads:5.3f}\t'
                    f'{pcr_dedup_reads}\t{pcr_dedup_reads/initial_reads:5.3f}\t'
                    f'{primer_trimmed_reads}\t{primer_trimmed_reads/initial_reads:5.3f}\t'
                    f'{mean_is}\t{med_is}\t{sd_is}\t'
                    f'{depth}\t{bases}\t{seqlength}\t{n_bases}\t'
                    f'{point_mutations}\t{inframe_gaps}\t{ooframe_gaps}\n'
    )

if __name__ == '__main__':
    main()

