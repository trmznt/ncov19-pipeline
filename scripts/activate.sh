
_script="$(readlink -f ${BASH_SOURCE[0]})"

## Delete last component from $_script ##
_mydir="$(dirname $_script)"

export NCOV19_PIPELINE="$(dirname $_mydir)"

