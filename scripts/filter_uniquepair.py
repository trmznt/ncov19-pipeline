#!/usr/bin/env python3

import argparse
import pysam
import bisect
import sys

def filter_unique_pair(infile, outfile, refname, min_match_len, max_nm, outstat):

    inbam = pysam.AlignmentFile(infile, 'rb')
    bam_header = inbam.header.copy().to_dict()

    outbam = pysam.AlignmentFile(outfile, 'wh', header=bam_header)

    statf = open(outstat, 'w') if outstat else None

    current_query_name = None
    current_reads = []
    uniqpair_count = 0
    read_count = 0

    for read in inbam:

        read_count += 1
        if read.reference_name != refname:
            continue

        if current_query_name is None:
            current_query_name = read.query_name
            current_reads.append(read)
            continue

        if current_query_name == read.query_name:
            current_reads.append(read)
            continue

        # read query name is not the same

        # check current reads
        if len(current_reads) == 2:
            uniqpair_count += check_and_write(current_reads, outbam, min_match_len, max_nm, statf)
        current_query_name = read.query_name
        current_reads = [ read ]

    if len(current_reads) == 2:
        uniqpair_count += check_and_write(current_reads, outbam, min_match_len, max_nm, statf)

    return uniqpair_count, read_count


#def write_stats('%d\t%d\t

def check_and_write(reads, outbam, min_match_len, max_nm, statf):

    read1 = reads[0]
    read2 = reads[1]
    if read1.is_reverse:
        read1, read2 = read2, read1
    if (read1.reference_end - read1.reference_start) < min_match_len:
        #print('%d' % (read1.reference_end - read1.reference_start))
        return 0
    if (read2.reference_end - read2.reference_start) < min_match_len:
        #print('%d' % (read2.reference_end - read2.reference_start))
        return 0
    if read1.get_tag('NM')/read1.template_length > max_nm or read2.get_tag('NM')/read2.template_length > max_nm:
        return 0

    if statf:
        statf.write('%d\t%d\t%d\t%d\t%s\n%d\t%d\t%d\t%d\t%s\n' %
                (   read1.query_alignment_start, read1.query_alignment_end, read1.mapping_quality, read1.get_tag('NM'), read1.query_name,
                    read2.query_alignment_start, read2.query_alignment_end, read2.mapping_quality, read2.get_tag('NM'), read2.query_name ))

    outbam.write(read1)
    outbam.write(read2)
    return 1


def init_argparser():

    p = argparse.ArgumentParser(description='Filter unique proper paired reads')
    p.add_argument('-o', '--outfile', required=True)
    p.add_argument('--min_match_len', type=int, default=40)
    p.add_argument('--max_nm', type=float, default=0.13)
    p.add_argument('--refname', default='NC_045512.2')
    p.add_argument('--outstat', default=None)
    p.add_argument('infile')

    return p


def main():

    args = init_argparser().parse_args()

    uniqpair, total_read = filter_unique_pair(args.infile, args.outfile, args.refname, args.min_match_len, args.max_nm, args.outstat)
    print('Found %d unique pairs (%d reads) out of %d reads (%f%%)' %
        (uniqpair, uniqpair*2, total_read, uniqpair*2/total_read*100))

if __name__ == '__main__':

    main()
