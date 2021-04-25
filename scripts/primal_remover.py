#!/usr/bin/env python3

# (c) Hidayat Trimarsanto <anto@eijkman.go.id>

# stand-alone primer sequence remover from amplicon-based sequencing
# some part of the code are copied from align_trim.py of ARTIC network

import argparse
import pysam
import bisect
import sys

from copy import copy
from typing import Iterable, Tuple, TypeVar

def _COUT( text ):
    print( text, file=sys.stdout )

def _CERR( text ):
    print( text, file=sys.stderr )

def cout( text ):
    _COUT(text)

def cerr( text ):
    _CERR( text )

def cexit( text ):
    cerr( text )
    sys.exit(1)

T = TypeVar("T")

def grouped(iterable: Iterable[T], n=2) -> Iterable[Tuple[T, ...]]:
    """s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), ..."""
    return zip(*[iter(iterable)] * n)



# consumesReference lookup for if a CIGAR operation consumes the reference sequence
consumesReference = [True, False, True, True, False, False, False, True]

# consumesQuery lookup for if a CIGAR operation consumes the query sequence
consumesQuery = [True, True, False, False, True, False, False, True]


class PrimerLookup(object):

    def __init__(self, bedfile=None):

        self.amplicon_table = {}
        self.lookup_table = {}
        self.left_keys = {}
        self.right_keys = {}
        # amplicon table contains position or location of primer binding set
        # per amplicon:
        # { 'refseq': [
        #        'amp_no': (left_start, left_end, right_start, right_end),
        #        ...
        #    ]
        # }

        if bedfile:
            self.parse_bedfile(bedfile)

    def get_ampinfo(self, ref, idx):
        return self.lookup_table[ref][idx]

    def parse_bedfile(self, bedfile):
        # bedfile format:
        # REF    START_POS    END_POS        PRIMERCODE    AMPLICON_ID        DIRECTION

        with open(bedfile) as fin:

            for line in fin:

                line = line.strip()
                if not line: continue
                ref, start_pos, end_pos, code, amp_id, direc = line.split('\t')
                start_pos = int(start_pos)
                end_pos = int(end_pos)

                try:
                    amp_ref = self.amplicon_table[ref]
                except KeyError:
                    amp_ref = self.amplicon_table[ref] = {}

                try:

                    amp_info = amp_ref[amp_id]
                    left_start, left_end, right_start, right_end = amp_info
                    if direc == '+':
                        if left_start < 0:
                            left_start = start_pos
                            left_end = end_pos
                        else:
                            left_start = min(left_start, start_pos)
                            left_end = max(left_end, end_pos)
                    else:
                        if right_start < 0:
                            right_start = start_pos
                            right_end = end_pos
                        else:
                            right_start = min(right_start, start_pos)
                            right_end = max(right_end, end_pos)

                    if right_end > 0 and left_start > right_end:
                        cexit("ERR: inconsistent primer position on amplicon %s"
                                % amp_id)
                    amp_ref[amp_id] = (left_start, left_end, right_start, right_end)

                except KeyError:
                    left_start = left_end = right_start = right_end = -1
                    if direc == '+':
                        left_start = start_pos
                        left_end = end_pos
                    else:
                        right_start = start_pos
                        right_end = end_pos

                    amp_ref[amp_id] = ( left_start, left_end, right_start, right_end )

        # create lookup table

        for ref in self.amplicon_table.keys():
            lookup_table =    [ (lf, le, rs, re, amp_id)
                                for (amp_id, (lf, le, rs, re))
                                        in self.amplicon_table[ref].items()
                            ]
            self.lookup_table[ref] = sorted(lookup_table)
            self.left_keys[ref] = [r[0] for r in self.lookup_table[ref]]
            self.right_keys[ref] = [r[3] for r in self.lookup_table[ref]]
            #from pprint import pprint
            #pprint(self.lookup_table[ref])


    def find_left_index(self, segment):
        """ find the first index of amplicon that can contain the segment """
        keys = self.left_keys[segment.reference_name]
        return bisect.bisect_right(keys, segment.reference_start )-1

    def find_right_index(self, segment):
        """ find the last index of amplicon that can contain the segment """
        keys = self.right_keys[segment.reference_name]
        return bisect.bisect_right(keys, segment.reference_end - 1)


def softmask_primer(segment, primer_pos, end, debug):
    """Soft mask an alignment to fit within primer start/end sites.
    Parameters
    ----------
    segment : pysam.AlignedSegment
        The aligned segment to mask
    primer_pos : int
        The position in the reference to soft mask up to (equates to the start/end position of the primer in the reference)
    end : bool
        If True, the segment is being masked from the end (i.e. for the reverse primer)
    debug : bool
        If True, will print soft masking info during trimming
    """
    # get a copy of the cigar tuples to work with
    cigar = copy(segment.cigartuples)

    # get the segment position in the reference (depends on if start or end of the segment is being processed)
    if not end:
        pos = segment.pos
    else:
        pos = segment.reference_end

    # process the CIGAR to determine how much softmasking is required
    eaten = 0
    while 1:

        # chomp CIGAR operations from the start/end of the CIGAR
        try:
            if end:
                flag, length = cigar.pop()
            else:
                flag, length = cigar.pop(0)
            if debug:
                print("Chomped a %s, %s" % (flag, length), file=sys.stderr)
        except IndexError:
            print(
                "Ran out of cigar during soft masking - completely masked read will be ignored", file=sys.stderr)
            break

        # if the CIGAR operation consumes the reference sequence, increment/decrement the position by the CIGAR operation length
        if (consumesReference[flag]):
            if not end:
                pos += length
            else:
                pos -= length

        # if the CIGAR operation consumes the query sequence, increment the number of CIGAR operations eaten by the CIGAR operation length
        if (consumesQuery[flag]):
            eaten += length

        # stop processing the CIGAR if we've gone far enough to mask the primer
        if not end and pos >= primer_pos and flag == 0:
            break
        if end and pos <= primer_pos and flag == 0:
            break

    # calculate how many extra matches are needed in the CIGAR
    extra = abs(pos - primer_pos)
    if debug:
        print("extra %s" % (extra), file=sys.stderr)
    if extra:
        if debug:
            print("Inserted a %s, %s" % (0, extra), file=sys.stderr)
        if end:
            cigar.append((0, extra))
        else:
            cigar.insert(0, (0, extra))
        eaten -= extra

    # softmask the left primer
    if not end:

        # update the position of the leftmost mappinng base
        segment.pos = pos - extra
        if debug:
            print("New pos: %s" % (segment.pos), file=sys.stderr)

        # if proposed softmask leads straight into a deletion, shuffle leftmost mapping base along and ignore the deletion
        if cigar[0][0] == 2:
            if debug:
                print(
                    "softmask created a leading deletion in the CIGAR, shuffling the alignment", file=sys.stderr)
            while 1:
                if cigar[0][0] != 2:
                    break
                _, length = cigar.pop(0)
                segment.pos += length

        # now add the leading softmask
        cigar.insert(0, (4, eaten))

    # softmask the right primer
    else:
        cigar.append((4, eaten))

    # check the new CIGAR and replace the old one
    if cigar[0][1] <= 0 or cigar[-1][1] <= 0:
        raise RuntimeError("invalid cigar operation created - possibly due to INDEL in primer")
    segment.cigartuples = cigar
    return


def softmask(read1, read2, ampinfo):

    # checking both 5'- and 3'- of each reads
    # -----------------===============-----
    # 5-XXXXXXXXXXXXXXXXXXX-3           <read1>
    #     3-XXXXXXXXXXXXXXXXXXX-5       <read2>
    #                  ^^^^^^^^ trimmed both read1 and read2
    #

    ls, le, rs, re, amp_id = ampinfo
    trimmed = False

    # check if read1 needs trimming from 3'
    if rs <= read1.reference_end <= re and rs <= read2.reference_end <= re:
        # trim read1 only if cigar is simple
        if len(read1.cigar) != 1:
            return None
        trim_size = read1.reference_end - rs
        if read1.cigar[0][1] <= trim_size:
            return None
        read1.cigar = [ (0, read1.cigar[0][1] - trim_size), (4, trim_size) ]
        trimmed |= True

    # check if read2 needs trimming from 5'
    if ls <= read2.reference_start <= le and ls <= read1.reference_start <= le:
        # trim read2 only if cigar is simple
        if len(read2.cigar) != 1:
            return None
        trim_size = le - read2.reference_start
        if read2.cigar[0][1] <= trim_size:
            return None
        read2.cigar = [ (4, trim_size), (0, read2.cigar[0][1] - trim_size) ]
        ref_start = read2.reference_start
        read2.reference_start += trim_size
        trimmed |= True

    if read1.reference_start < le:
        #print('8< softmask read1 >8')
        if read1.reference_end <= le+1:
            # whole read is consumed by primer; remove this pair
            return None
        softmask_primer(read1, le+1, False, False)
        trimmed |= True

    if read2.reference_end > rs:
        #print('8< softmask read2 >8')
        if read2.reference_start >= rs-1:
            # whole read is consumed by primer; remove this pair
            return None
        softmask_primer(read2, rs-1, True, False)
        trimmed |= True

    return trimmed


def trim_ligation(read1, read2, primer_lookup, counter = {}):
    """ return False if no trimming, None if pairs are invalid, True if trimmed
    """

    l_idx = primer_lookup.find_left_index(read1)
    r_idx = primer_lookup.find_right_index(read2)
    refname = read1.reference_name
    lookup = primer_lookup.lookup_table[refname]

    if l_idx != r_idx:
        # with adapter ligation, reads must be product of just
        # a single amplicon set, otherwise primer mix has happened
        # XXX: still need to check if we can include the reads

        (ls1, le1, rs1, re1, ampid1) = lookup[l_idx]
        (ls2, le2, rs2, re2, ampid2) = lookup[r_idx]
        read1.set_tag('XP', -1)
        read2.set_tag('XP', -1)
        if read1.reference_start < le1 or read2.reference_end > rs2:
            try:
                softmask_primer(read1, le1, False, False)
                softmask_primer(read2, rs2, True, False)
            except Exception as err:
                cerr('[discarding possible mixup primers amplicons #%d <> #%d]' % (l_idx+1, r_idx+1))
                cerr('[ATT: %d %d %d %d]' % (read1.reference_start, read1.reference_end, read2.reference_start, read2.reference_end))
                cerr( str(err) )
                return None
            return True
        return False

    (ls, le, rs, re, ampid) = lookup[l_idx]

    # count amplicon numbers as well

    try:
        counter[ampid] += 1
    except:
        counter[ampid] = 1

    if read1.reference_start > le or read2.reference_end < rs:
        # read1 and read2 do not start/end at primer binding site
        #cerr('[WARN: reads too short for amplicon #%d: %s bp (expected: %d bp)]' %
        #        (l_idx, read2.reference_end-read1.reference_start, re-ls))
        read1.set_tag('XP', l_idx + 1)
        read2.set_tag('XP', l_idx + 1)

        return False

    softmask_primer(read1, le+1, False, False)
    softmask_primer(read2, rs-1, True, False)
    read1.set_tag('XP', l_idx + 1)
    read2.set_tag('XP', l_idx + 1)

    return True


def trim_tagmentation(read1, read2, primer_lookup, counter=None):
    """ return False if no trimming, None if pairs are invalid, True if trimmed
    """

    l_idx = primer_lookup.find_left_index(read1)
    r_idx = primer_lookup.find_right_index(read2)
    refname = read1.reference_name
    lookup = primer_lookup.lookup_table[refname]
    #print(read1.reference_start, '-', read2.reference_end, ' <> ',
    #        primer_lookup.get_ampinfo(read1.reference_name, l_idx),
    #        primer_lookup.get_ampinfo(read2.reference_name, r_idx) )

    # sanity check with prev and next amplicons
    possible_amplicons = []
    for idx in range(max(0, l_idx-1), min(r_idx+1, len(lookup))):
        ls, le, rs, re, amp_id = lookup[idx]
        if ls <= read1.reference_start and read2.reference_end <= re:
            possible_amplicons.append(idx)

    if len(possible_amplicons) > 1:
        # we have more than 1 amplicon that can yield this pair
        print('WARN: %d <> %d' % (read1.reference_start, read2.reference_end))
        for idx in possible_amplicons:
            print('>>', lookup[idx])
    elif len(possible_amplicons) == 0:
        r_idx = min(len(lookup)-1, r_idx)
        #cerr('ERR: <%d, %d>' % (l_idx, r_idx))
        ls1, le1, rs1, re1, amp_id1 = lookup[l_idx]
        ls2, le2, rs2, re2, amp_id2 = lookup[r_idx]
        read1.set_tag('XP', -1)
        read2.set_tag('XP', -1)
        return softmask(read1, read2, (ls1, le1, rs2, re2, -1))
    else:
        #print('== %d <> %d ==' % (read1.reference_start, read2.reference_end),
        #    lookup[possible_amplicons[0]])
        read1.set_tag('XP', possible_amplicons[0] + 1)
        read2.set_tag('XP', possible_amplicons[0] + 1)
        return softmask(read1, read2, lookup[possible_amplicons[0]])

    return None


def trim_primers(segments, primer_lookup, outfile, trimmer_func=None):

    read_pairs = trimmed_pairs = untrimmed_pairs = invalid_pairs = indel_pairs = rf_pairs = unmapped_pairs = 0
    amplicon_counter = {}

    for (read1, read2) in grouped(segments):
        read_pairs += 1
        if read1.query_name != read2.query_name:
            cexit("ERR: reads were not sorted by name or reads were not in pairs: %s >< %s" %
                    (read1.query_name, read2.query_name))

        if read1.is_unmapped or read2.is_unmapped:
            # unmapped read pairs
            unmapped_pairs += 1
            continue

        # sanity and ordering check
        if read1.reference_end < read1.reference_start or read2.reference_end < read2.reference_start:
            cexit("ERR: found inconsistent read pair direction")

        if read1.is_reverse:
            read1, read2 = read2, read1

        if read1.reference_start > read2.reference_end:
            # possible RF alignment, for now just exit
            cerr("ERR: found possible RF pairs")
            rf_pairs += 1
            continue

        if read1.reference_end > read2.reference_end:

            # read1 passes read2.end, just softmask read1
            # if cigar is complicated, just remove this read pair
            if len(read1.cigar) > 1 or read1.cigar[0][0] != 0:
                cerr("ERR: complicated cigar: %s" % read1.cigarstring)
                cerr('[ATT: %d %d %d %d]' % (read1.reference_start, read1.reference_end, read2.reference_start, read2.reference_end))
                invalid_pairs += 1
                continue

            # modify cigar
            trim_size = read1.reference_end - read2.reference_end
            read1.cigar = [ (0, read1.cigar[0][1] - trim_size), (4, trim_size)]

        if read1.reference_start > read2.reference_start:

            # read2 passes read1.reference_start, just softmask read2
            # only modify cigar if the last cigar is either M or S
            if len(read2.cigar) > 1 or read2.cigar[0][0] != 0:
                cexit("ERR: complicated cigar: %s" % read2.cigarstring)
                cerr('[ATT: %d %d %d %d]' % (read1.reference_start, read1.reference_end, read2.reference_start, read2.reference_end))
                invalid_pairs += 1
                continue

            # modify cigar
            trim_size = read1.reference_start - read2.reference_start
            read2.cigar = [ (4, trim_size), (0, read2.cigar[0][1] - trim_size) ]

        try:
            res = trimmer_func(read1, read2, primer_lookup, amplicon_counter)
        except BaseException as inst:
            raise inst
            print(inst)
            indel_pairs += 1
            continue
        if res is True:
            trimmed_pairs += 1
            #print(' 8< %s >8' % (read1.query_name))
        elif res is False:
            untrimmed_pairs += 1
        elif res is None:
            # discard invalid pairs
            invalid_pairs += 1
            continue
        else:
            raise RuntimeError('unknown trimmer_func() result')
        outfile.write(read1)
        outfile.write(read2)

    print(amplicon_counter)
    return (read_pairs, trimmed_pairs, untrimmed_pairs, invalid_pairs, rf_pairs, indel_pairs, unmapped_pairs, amplicon_counter)


def primal_remover(args):

    primer_lookup = PrimerLookup(bedfile = args.bedfile)

    segments = pysam.AlignmentFile(args.infile, 'rb')
    bam_header = segments.header.copy().to_dict()

    if args.method == 'tagment':
        trim_func = trim_tagmentation
    elif args.method == 'ligate':
        trim_func = trim_ligation
    else:
        cexit('[ERR: unknown trimming method: %s' % args.method)

    outfile = pysam.AlignmentFile(args.outfile, 'wh', header=bam_header)

    read_pairs, trimmed_pairs, untrimmed_pairs, invalid_pairs, rf_pairs, indel_pairs, unmapped_pairs, amplicon_counter = trim_primers(
                            segments, primer_lookup, outfile, trim_func)

    segments.close()
    outfile.close()

    print(amplicon_counter)
    if args.outcount and amplicon_counter:
        with open(args.outcount, 'w') as countfile:
            for k, v in amplicon_counter.items():
                countfile.write('%s\t%d\n' % (k, v))

    if args.logfile:
        with open(args.logfile, 'w') as logfile:
            logfile.write(  'Read pairs: %d\n'
                            'Trimmed pairs: %d\n'
                            'Untrimmed pairs: %d\n'
                            'Invalid pairs: %d\n'
                            'RF pairs: %d\n'
                            'Indel pairs:%d\n'
                            'Unmapped_pairs:%d\n'
                % (read_pairs, trimmed_pairs, untrimmed_pairs, invalid_pairs,
                    rf_pairs, indel_pairs, unmapped_pairs)
            )


def init_argparser(p=None):

    if p is None:
        p = argparse.ArgumentParser(
                description = 'primal_remover -- remove primalscheme primers')

    p.add_argument('-o', '--outfile', default='outfile.bam')
    p.add_argument('-m', '--method', default='tagment')
    p.add_argument('--logfile', default='')
    p.add_argument('--outcount', default='')
    p.add_argument('--bedfile', required=True)
    p.add_argument('infile')

    return p

def main():
    p = init_argparser()
    args = p.parse_args()
    primal_remover(args)


if __name__ == '__main__':
    main()
