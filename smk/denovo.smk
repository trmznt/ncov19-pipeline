import os

try:
    rootdir = os.environ['NCOV19_PIPELINE']
except KeyError:
    raise RuntimeError('NCOV19_PIPELINE environment is not set up!')

refseq = rootdir + "/ref/nc_045512.fasta",
sample_id = os.path.basename( os.getcwd() )
adapters = {'nextera': ['CTGTCTCTTATACACATCT', 'CTGTCTCTTATACACATCT'] }

include: "trim_reads.smk"

rule all:
    input:
        "spades/scaffolds.fasta"


rule denovo_assembling:
    input:
        read1 = "reads/trimmed.R1.fastq.gz",
        read2 = "reads/trimmed.R2.fastq.gz"
    output:
        "spades/scaffolds.fasta"
    shell:
        "spades.py --pe1-1 {input.read1} --pe1-2 {input.read2} -o spades --rnaviral"
