#! /usr/bin/bash
# Chunyu Zhao 2020-10-21
set -e

if [ $# -ne 2 ]; then
    echo "Usage: $0 COORDS_FILE OUT_FILE"
    exit 1
fi

coords_file="$1"
outfile="$2"


awk 'BEGIN {FS=" \\|"} ; {print $7}' ${coords_file} | sed -e 's/^ *//' | awk '{print $1, $2}' | grep UHGG | sort -u -k1,2 > ${outfile}
