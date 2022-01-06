#!/bin/sh
#wrapper to run snakemake within a sample directory to merge consensus bams of the sample that underwent multi
ple runs
#the layout directory follows:
#scripts/
#smk/
#data
#   /ngs-XXX/
#       /config.yaml
#       /XXXX/
#       /XXXX/sample1
#       /XXXX/sample2
#
#run this wrapper within data/NGS-XXXX/XXXX

_script="$(readlink -f ${BASH_SOURCE[0]})"

## Delete last component from $_script ##
_mydir="$(dirname $_script)"
export NCOV19_PIPELINE="$(dirname $_mydir)"

snakemake -s ${NCOV19_PIPELINE}/smk/pipeline_merge.smk --configfile ../../config.yaml -j1 $1
