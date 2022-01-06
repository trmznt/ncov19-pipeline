def get_cons_bams(wildcards):
    reps = sorted(glob_wildcards("runs/run-" + "{rep}/maps").rep)
    return expand( "runs/run-{rep}/maps/consensus.bam", rep = reps)

rule merge_cons:
    """
    Merge consensus bam files for mutiple runs into one for the given sample.
    If the sample has only underwent one run, a symlink will be created.
    """
    input:
	cons = get_cons_bams
    output:
	"maps/consensus.bam"
    run:
        if len(input.cons) > 1:
            shell("samtools merge -@8 {output} {input.cons}")
        else:
            shell("cp {input.cons} {output} && touch -h {output}")
