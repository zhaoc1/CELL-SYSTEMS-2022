#!/bin/bash
set -e
set -x

if [ $# -ne 1 ]; then
    echo "Usage: $0 SPECIES"
    exit 1
fi

species="$1"

proj_dir="/mnt/cz/20220602_close_species/${species}"

reads_dir="${proj_dir}/reads"


cd $reads_dir

grep "^@" */2_trimmomatic/cov_20_paired_1.fastq | cut -d'-' -f1 | sort | uniq -c > ${proj_dir}/R1.readcounts
