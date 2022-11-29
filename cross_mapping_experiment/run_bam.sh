#!/bin/bash
set -e
set -x

if [ $# -ne 2 ]; then
    echo "Usage: $0 SPECIES BIN_ON"
    exit 1
fi

species="$1"
on_bin="$2"

proj_dir="/mnt/cz/20220602_close_species/${species}"

#midas_outdir="${proj_dir}/${on_bin}/4_midas/${off_bin}"
midas_outdir="${proj_dir}/${on_bin}/4_midas_local"

outdir="/mnt/cz/20220602_close_species/bamcounts_bowtie2_local/${species}"
mkdir -p ${outdir}

ls $midas_outdir/*/cov_20X/temp/snps/*/*.bam | xargs -Ixx bash -c \
  'samtools view xx | cut -f1 | cut -d'-' -f1 | sort | uniq -c | awk "$0"' '{print "xx", $0}' > $outdir/$on_bin.bamcounts
