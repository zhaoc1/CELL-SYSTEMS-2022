#!/bin/bash
set -e
set -x

if [ $# -ne 3 ]; then
    echo "Usage: $0 SPECIES BIN_ON BIN_OFF"
    exit 1
fi

species="$1"
on_bin="$2"
off_bin="$3"

proj_dir="/mnt/cz/20220602_close_species/${species}"
data_dir="${proj_dir}/${on_bin}/2_reads/${off_bin}"

################### READS
r1="${data_dir}/cov_20_1.fastq.gz"
r2="${data_dir}/cov_20_2.fastq.gz"

####################BT2_INDEX
bt2_dir="${proj_dir}/${on_bin}/3_bt2db/${off_bin}"
repgenome="${bt2_dir}/repgenomes.fa"
repspecies="${bt2_dir}/repgenome.species"

bwa index ${repgenome}

midas_outdir="${proj_dir}/${on_bin}/4_midas_bwa/${off_bin}"
mkdir -p ${midas_outdir}

bwa mem -M -t 8 ${repgenome} ${r1} ${r2} | samtools view --threads 8 -b - | \
  samtools sort --threads 8 -o ${midas_outdir}/repgenomes.bam
samtools index ${midas_outdir}/repgenomes.bam

bam_outdir="${proj_dir}/bwa_bamcounts/${on_bin}"
mkdir -p ${bam_outdir}

samtools view ${midas_outdir}/repgenomes.bam | cut -f1,3 | awk -F'\t' '{sub(/-.+$/,"",$1)}1' OFS='\t' | sort | uniq -c | awk '{print $1, $2, $3}' OFS='\t' > ${bam_outdir}/${off_bin}.bamcounts

echo "${bam_outdir}/${off_bin}.bamcounts"
echo "${species}:${on_bin}:${off_bin} DONE"
