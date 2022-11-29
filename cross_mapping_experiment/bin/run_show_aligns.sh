#! /usr/bin/bash
# Chunyu Zhao 2020-10-21
set -e

if [ $# -ne 2 ]; then
    echo "Usage: $0 DELTA_FILE CSTRING"
    exit 1
fi

delta="$1"
mfs="$2"
sarr=($mfs)

ref=${sarr[0]}
qry=${sarr[1]}

workdir=`dirname "${delta}"`
mkdir -p ${workdir}/${ref}

show-aligns -r -m 0 ${delta} ${ref} ${qry} | sed '/^$/d' > ${workdir}/${ref}/${qry}.aln
