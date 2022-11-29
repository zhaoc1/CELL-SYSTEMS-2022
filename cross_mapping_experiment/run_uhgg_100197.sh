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

data_dir="${proj_dir}/${on_bin}/2_reads/${on_bin}"
midasdb_dir="/mnt/cz/midas_uhgg0"


r1="${data_dir}/cov_20_1.fastq.gz"
r2="${data_dir}/cov_20_2.fastq.gz"


bt2_dir="/mnt/cz/midas_uhgg0/bowtie2_indexes"
repspecies="${bt2_dir}/repgenomes.species"

echo $repspecies


midas_outdir="${proj_dir}/${on_bin}/uhgg_all"
mkdir -p ${midas_outdir}
python3 -m iggtools midas_run_snps \
    --sample_name cov_20X -1 ${r1} -2 ${r2} \
    --midasdb_name uhgg --midasdb_dir ${midasdb_dir} \
    --prebuilt_bowtie2_indexes ${bt2_dir}/repgenomes \
    --prebuilt_bowtie2_species ${repspecies} \
    --advanced --analysis_ready --fragment_length 1000 --select_threshold=-1 \
    --num_cores 32 ${midas_outdir}

echo "${species}:${on_bin}:${off_bin} DONE"
