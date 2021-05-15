
adapters = {
	'nextera': ['CTGTCTCTTATACACATCT', 'CTGTCTCTTATACACATCT'],
	'truseq': ['AGATCGGAAGAGCACACGTCTGAACTCCAGTCA', 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'],
	'ultraii': ['AGATCGGAAGAGCACACGTCTGAACTCCAGTCA', 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'],
	'ultraiifs': ['AGATCGGAAGAGCACACGTCTGAACTCCAGTCA', 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'],
	'dnaprep': ['CTGTCTCTTATACACATCT+ATGTGTATAAGAGACA', 'CTGTCTCTTATACACATCT+ATGTGTATAAGAGACA'],
}
minlen = int(config['minlen']) if 'minlen' in config else (int(config['read_length']) / 3)
maxlen = int(config['maxlen']) if 'maxlen' in config else 0

def adapter_arguments(libprep):
    arguments = []
    adapter_1, adapter_2 = adapters[libprep]
    arguments += ['-a %s' % x for x in adapter_1.split('+')]
    arguments += ['-A %s' % x for x in adapter_2.split('+')]
    return ' '.join(arguments)

def is_nextseq():
    return config['instrument'].lower().startswith('nextseq')


rule optical_dedup:
    input:
        read1 = "reads/raw.R1.fastq.gz",
        read2 = "reads/raw.R2.fastq.gz"
    output:
        dedup1 = temp("reads/dedup.R1.fastq.gz"),
        dedup2 = temp("reads/dedup.R2.fastq.gz")
    log: "logs/optical_dedup.log"
    shell:
        "clumpify.sh in={input.read1} in2={input.read2} out1={output.dedup1} out2={output.dedup2} dedupe optical %s 2> {log}"
        % ('spany adjacent' if is_nextseq() else '')


rule adapter_trimming:
    input:
        read1 = "reads/dedup.R1.fastq.gz",
        read2 = "reads/dedup.R2.fastq.gz"
    output:
        trimmed1 = "reads/trimmed.R1.fastq.gz",
        trimmed2 = "reads/trimmed.R2.fastq.gz"
    log: "logs/adapter_trimming.log"
    run:
        shell("cutadapt %s -j 18 %s -m %d -q %d -O 3 %s -o %s -p %s %s %s > %s"
	% ('--nextseq-trim 20' if is_nextseq() else '',
            ('--length %d' % maxlen) if maxlen > 0 else '',
            minlen,
            config['min_read_qual'],
            adapter_arguments(config['libprep'].lower()),
            output.trimmed1, output.trimmed2, input.read1, input.read2,
            log
	))

