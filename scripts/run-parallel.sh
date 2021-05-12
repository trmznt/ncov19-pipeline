#!/bin/sh
_script="$(readlink -f ${BASH_SOURCE[0]})"

## Delete last component from $_script ##
_mydir="$(dirname $_script)"

parallel --eta --sshlogin 8/cu-4,8/cu-5 --workdir $PWD/{} ${_mydir}/run-pipeline.sh all ::: `ls`

