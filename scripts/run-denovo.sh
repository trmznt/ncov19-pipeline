#!/bin/sh
# wrapper to run snakemake within sample directory
# the layout directory follows:
# scripts/
# smk/
# data/
#     /ngs-XXX/
#             /config.yaml
#             /XXXX/
#             /XXXX/sample1
#             /XXXX/sample2
#
# run this wrapper within data/ngs-XXX/XXXX

_script="$(readlink -f ${BASH_SOURCE[0]})"

## Delete last component from $_script ##
_mydir="$(dirname $_script)"
export NCOV19_PIPELINE="$(dirname $_mydir)"

snakemake -s ${NCOV19_PIPELINE}/smk/denovo.smk --configfile ../../config.yaml -j1 $1

