
import os

try:
    rootdir = os.environ['NCOV19_PIPELINE']
except KeyError:
    raise RuntimeError('NCOV19_PIPELINE environment is not set up!')

refseq = rootdir + "/ref/nc_045512.fasta"
gff = rootdir + "/ref/GCF_009858895.2_ASM985889v3_genomic.gff"
sample_id = os.path.basename( os.getcwd() )
adapters = {'nextera': ['CTGTCTCTTATACACATCT', 'CTGTCTCTTATACACATCT'] }
mode = 'tagment' if config['libprep'].lower() in ['nextera'] else 'ligate'
min_depth = config['min_depth']

rule all:
    input:
        "cons/blastn-consensus.txt",
        "maps/depth-counter.txt",
        "maps/trimmed-stats.txt.gz",
        "cons/variants.tsv"



rule optical_dedup:
    input:
        read1 = "reads/raw.R1.fastq.gz",
        read2 = "reads/raw.R2.fastq.gz"
    output:
        dedup1 = temp("reads/dedup.R1.fastq.gz"),
        dedup2 = temp("reads/dedup.R2.fastq.gz")
    log: "logs/optical_dedup.log"
    shell:
        "{rootdir}/scripts/bbmap/clumpify.sh in={input.read1} in2={input.read2} out1={output.dedup1} out2={output.dedup2} dedupe optical spany adjacent 2> {log}"


rule adapter_trimming:
    input:
        read1 = "reads/dedup.R1.fastq.gz",
        read2 = "reads/dedup.R2.fastq.gz"
    output:
        trimmed1 = "reads/trimmed.R1.fastq.gz",
        trimmed2 = "reads/trimmed.R2.fastq.gz"
    log: "logs/adapter_trimming.log"
    run:
        shell("cutadapt %s -j 18  -m %d -q %d -O 4 -u 1 -a N%s -U 1 -A N%s -o %s -p %s %s %s 2> %s"
	% ('--nextseq-trim 20' if config['instrument'].lower().startswith('nextseq') else '',
            int(config['read_length'] / 3),
            config['min_read_qual'],
            adapters[config['libprep'].lower()][0],
            adapters[config['libprep'].lower()][1],
            output.trimmed1, output.trimmed2, input.read1, input.read2,
            log
	))


rule mapping:
    input:
        read1 = "reads/trimmed.R1.fastq.gz",
        read2 = "reads/trimmed.R2.fastq.gz"
    output:
        "maps/all.bam"
    log: "logs/minimap2.log"
    shell:
        "minimap2 -ax sr -t 16 {refseq} {input.read1} {input.read2} 2> {log} | samtools view -bS - > {output}"


rule map_filtering:
    input:
        "maps/all.bam"
    output:
        temp("maps/mapped.bam")
    shell:
        "samtools view -b -f 2 -F 12 {input} > {output}"


rule fix_mate:
    input:
        "maps/mapped.bam"
    output:
        temp("maps/mapped.fixmated.bam")
    shell:
        "samtools sort -n -@16 {input} | samtools fixmate -r -m - - | samtools sort -@16 -o {output} -"


rule pcr_dedup:
    input:
        "maps/mapped.fixmated.bam"
    output:
        "maps/mapped.dedup.bam"
    log: "logs/markdup.log"
    shell:
        "samtools markdup -r -s {input} {output} 2> {log}"


rule pair_filtering:
    input:
        "maps/mapped.dedup.bam"
    output:
        temp("maps/unique_pairs.by-name.sam")
    log: "logs/unique_pairs.log"
    shell:
        "(samtools sort -n -@16 {input} | {rootdir}/scripts/filter_uniquepair.py --min_match_len 23 --max_nm 0.29 --outstat maps/unique_pairs.txt -o {output} - ) && gzip maps/unique_pairs.txt"


if config['primer_trimmer'].lower() == 'primal_remover':

    rule pair_compressing:
        input:
            "maps/unique_pairs.by-name.sam"
        output:
            temp("maps/unique_pairs.by-name.bam")
        shell:
            "samtools view -bS -o {output} {input}"


    rule primer_trimming:
        input:
            "maps/unique_pairs.by-name.bam"
        output:
            temp("maps/primers_trimmed.sam")
        log: "logs/primal_remover.log"
        shell:
            "{rootdir}/scripts/primal_remover.py -m {mode} --bedfile {rootdir}/ref/nCoV-2019-pr.bed --logfile {log} -o {output} {input}"


    rule trimmed_compressing:
        input:
            "maps/primers_trimmed.sam"
        output:
            "maps/primers_trimmed.bam"
        shell:
            "samtools fixmate -m {input} - | samtools sort -@16 -o {output} -"

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


rule depth_stat:
    input:
        "maps/primers_trimmed.bam"
    output:
        "maps/trimmed-depths.txt.gz"
    shell:
        "samtools depth {input} | gzip > {output}"

rule depth_counter:
    input:
        "maps/trimmed-depths.txt.gz"
    output:
        "maps/depth-counter.txt"
    shell:
        "{rootdir}/scripts/depth-counter.py {input} {sample_id} > {output}"


rule map_stat:
    input:
        "maps/primers_trimmed.bam"
    output:
        "maps/trimmed-stats.txt.gz"
    shell:
        "samtools stats {input} | gzip > {output}"


rule consensus:
    input:
        "maps/primers_trimmed.bam"
    output:
        "cons/consensus.fa"
    log: "logs/ivar.log"
    shell:
        "samtools mpileup -A -Q 0 {input} | ivar consensus -p cons/consensus -q 17 -t 0 -m {min_depth} 2> {log}"

rule variants:
    input:
        "maps/primers_trimmed.bam"
    output:
        "cons/variants.tsv"
    shell:
        "samtools mpileup -aa -A -d 10000 -B -Q 0 {input} | ivar variants -p cons/variants -q 17 -t 0.5 -r {refseq} -g {gff}"


rule blast_reporting:
    input:
        "cons/consensus.fa"
    output:
        "cons/blastn-consensus.txt"
    shell:
        "blastn -query {input} -subject {refseq} > {output}"

