import os

try:
    rootdir = os.environ['NCOV19_PIPELINE']
except KeyError:
    raise RuntimeError('NCOV19_PIPELINE environment is not set up!')

configfile: rootdir + "/smk/config.yaml"

refseq = rootdir + "/ref/nc_045512.fasta"
gff = rootdir + "/ref/NC_045512.2.sorted.gff3"
sample_id = os.path.basename( os.getcwd() )
min_depth = config['min_depth']

include: "merge_bam.smk"
include: "consensus.smk"

rule all:
    input:
	"maps/consensus.bam",
        "maps/depth-counter.txt",
        "maps/consensus-depths.png",
        "cons/blastn-consensus.txt",
        "cons/consensus.fas",
        "cons/variants.tsv",
        "logs/stats.tsv",

rule stats:
    input:
	"maps/depth-counter.txt",
        "cons/consensus.fas",
        "cons/variants.tsv",
    output:
	"logs/stats.tsv"
    shell:
	"{rootdir}/scripts/collect_stats.py {sample_id} > {output}"
