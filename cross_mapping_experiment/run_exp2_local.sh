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
midasdb_dir="${proj_dir}/${on_bin}/1_dbs/${off_bin}"

mkdir -p ${data_dir}

################### READS
r1="${data_dir}/cov_20_1.fastq.gz"
r2="${data_dir}/cov_20_2.fastq.gz"

######################MIDAS
bt2_dir="${proj_dir}/${on_bin}/3_bt2db/${off_bin}"
midas_outdir="${proj_dir}/${on_bin}/4_midas_local/${off_bin}"
repspecies="${bt2_dir}/repgenome.species"

mkdir -p ${midas_outdir}

python3 -m iggtools midas_run_snps \
    --sample_name cov_20X -1 ${r1} -2 ${r2} \
    --midasdb_name uhgg --midasdb_dir ${midasdb_dir} \
    --prebuilt_bowtie2_indexes ${bt2_dir}/repgenomes \
    --prebuilt_bowtie2_species ${repspecies} \
    --aln_mode local \
    --advanced --analysis_ready --fragment_length 1000 --select_threshold=-1 \
    --num_cores 8 ${midas_outdir}


echo "${species}:${on_bin}:${off_bin} DONE"
