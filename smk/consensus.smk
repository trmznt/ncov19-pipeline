
rule consensus_all:
   input:
       "maps/trimmed-depths.png"
       "cons/consensus.fa"
       "cons/variants.tsv"
       "cons/blastn-consensus.txt" 


rule plot_depths:
    input:
        "maps/primers_trimmed.bam"
    output:
        plotfile = "maps/trimmed-depths.png",
        gapfile = "maps/trimmed-gaps.txt",
        depthfile = "maps/depth-counter.txt"
    shell:
        "{rootdir}/scripts/depthplot.py -t {sample_id} -d {min_depth} --stat_insert --outgap {output.gapfile} --outdepth {output.depthfile} -o {output.plotfile} {input}"


rule consensus:
    input:
        "maps/primers_trimmed.bam"
    output:
        "cons/consensus.fa"
    log: "logs/ivar.log"
    shell:
        "samtools mpileup -d 10000 -Q 0 {input} | ivar consensus -p cons/consensus -q 17 -t 0 -m {min_depth} 2> {log}"

rule variants:
    input:
        "maps/primers_trimmed.bam"
    output:
        "cons/variants.tsv"
    shell:
        "samtools mpileup -a -d 10000 -Q 0 {input} | ivar variants -p cons/variants -q 17 -t 0.5 -r {refseq} -g {gff}"


rule blast_reporting:
    input:
        "cons/consensus.fa"
    output:
        "cons/blastn-consensus.txt"
    shell:
        "blastn -query {input} -subject {refseq} > {output}"

