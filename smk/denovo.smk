import os

try:
    rootdir = os.environ['NCOV19_PIPELINE']
except KeyError:
    raise RuntimeError('NCOV19_PIPELINE environment is not set up!')

refseq = rootdir + "/ref/nc_045512.fasta",
sample_id = os.path.basename( os.getcwd() )
adapters = {'nextera': ['CTGTCTCTTATACACATCT', 'CTGTCTCTTATACACATCT'] }

rule all:
    input:
        read1 = "reads/trimmed.R1.fastq.gz",
        read2 = "reads/trimmed.R2.fastq.gz"

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

rule denovo_assembling:
    input:
        read1 = "reads/trimmed.R1.fastq.gz",
        read2 = "reads/trimmed.R2.fastq.gz"
    output:
        "spades/scaffold.fasta"