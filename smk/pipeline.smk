
import os

try:
    rootdir = os.environ['NCOV19_PIPELINE']
except KeyError:
    raise RuntimeError('NCOV19_PIPELINE environment is not set up!')

refseq = rootdir + "/ref/nc_045512.fasta"
refmap = rootdir + "/ref/nc_045512.mmi"
gff = rootdir + "/ref/NC_045512.2.sorted.gff3"
sample_id = os.path.basename( os.getcwd() )
mode = 'tagment' if config['libprep'].lower() in ['nextera', 'ultraiifs', 'dnaprep'] else 'ligate'
min_depth = config['min_depth']

include: "trim_reads.smk"
include: "consensus.smk"


rule all:
    input:
        "maps/primers_trimmed.bam",
        "maps/depth-counter.txt",
        "maps/trimmed-depths.png",
        "cons/blastn-consensus.txt",
        "cons/consensus.fa",
        "cons/variants.tsv"


rule mapping:
    input:
        read1 = "reads/trimmed.R1.fastq.gz",
        read2 = "reads/trimmed.R2.fastq.gz"
    output:
        bam = "maps/mapped.bam",
        stat = "maps/unique_pairs.txt.gz"
    log:
        log1 = "logs/minimap2.log",
        log2 = "logs/unique_pairs.log"
    shell:
        "minimap2 -ax sr -t 16 --secondary=no {refmap} {input.read1} {input.read2} 2> {log.log1}"
        " | {rootdir}/scripts/filter_uniquepair.py --min_match_len 23 --max_nm 0.29 --outstat {output.stat} -o - - 2> {log.log2}"
        " | samtools fixmate -r -m - {output.bam}"


#rule map_filtering:
#    input:
#        "maps/all.bam"
#    output:
#        temp("maps/mapped.bam")
#    shell:
#        "samtools view -b -f 2 -F 12 {input} > {output}"



#rule fix_mate:
#    input:
#        "maps/all.bam"
#    output:
#        temp("maps/mapped.fixmated.bam")
#    log: "logs/fixmate.log"
#    shell:
#         "samtools fixmate -r -m {input} - 2> {log} | samtools sort -@16 -o {output} -"
#        "samtools sort -n -@16 {input} | samtools fixmate -r -m - - | samtools sort -@16 -o {output} -"

if mode == 'tagment':
    rule pcr_dedup:
        input:
            "maps/mapped.bam"
        output:
            temp("maps/unique_pairs.by-name.bam")
        log: "logs/markdup.log"
        shell:
            "samtools sort -@8 {input} | samtools markdup -r -s - - 2> {log} | samtools sort -n -@8 -o {output} -"


#    rule pair_filtering:
#        input:
#            "maps/mapped.dedup.bam"
#        output:
#            temp("maps/unique_pairs.by-name.sam")
#        log: "logs/unique_pairs.log"
#        shell:
#            "(samtools sort -n -@16 {input} | {rootdir}/scripts/filter_uniquepair.py --min_match_len 23 --max_nm 0.29 --outstat maps/unique_pairs.txt -o {output} - ) && gzip -f maps/unique_pairs.txt"

elif mode == 'ligate':

    rule link_map:
        input:
            "maps/mapped.bam"
        output:
            temp("maps/unique_pairs.by-name.bam")
        shell:
            "ln -s {input} {output}"

else:
    raise RuntimeError("ERR: mode %s is unknown" % mode)


if config['primer_trimmer'].lower() == 'primal_remover':

#    rule pair_compressing:
#        input:
#            "maps/unique_pairs.by-name.sam"
#        output:
#            temp("maps/unique_pairs.by-name.bam")
#        shell:
#            "samtools view -bS -o {output} {input}"


    rule primer_trimming:
        input:
            "maps/unique_pairs.by-name.bam"
        output:
            temp("maps/primers_trimmed.bam"),
        log: "logs/primal_remover.log"
        shell:
            "{rootdir}/scripts/primal_remover.py -m {mode} --bedfile {rootdir}/ref/nCoV-2019-pr.bed --logfile {log}  -o - {input}"
            " | samtools fixmate -m - - | samtools sort -@8 -o {output} -"


#    rule trimmed_compressing:
#        input:
#            "maps/primers_trimmed.sam"
#        output:
#            "maps/primers_trimmed.bam"
#        shell:
#            "samtools fixmate -m {input} - | samtools sort -@16 -o {output} -"


elif config['primer_trimmer'].lower() == 'ivar':

    rule sort_compressing:
        input:
            "maps/unique_pairs.by-name.sam"
        output:
            temp("maps/unique_pairs.by-pos.bam")
        shell:
            "samtools sort -@16 -o {output} {input}"

    rule ivar_primer_trimming:
        input:
            "maps/unique_pairs.by-pos.bam"
        output:
            temp("maps/primers_trimmed-unsorted.bam")
        shell:
            "ivar trim -b {rootdir}/ref/nCoV-2019-ivar.bed -p maps/primers_trimmed-unsorted -i {input} -q 15 -m 33 -s 4"

    rule trimmed_sorting:
        input:
            "maps/primers_trimmed-unsorted.bam"
        output:
            "maps/primers_trimmed.bam"
        shell:
            "samtools sort -@16 -o {output} {input}"

else:
    raise RuntimeError('ERR: primer_trimmer is not defined!')



