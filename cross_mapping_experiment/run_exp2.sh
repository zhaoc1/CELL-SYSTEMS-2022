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

on1="${proj_dir}/reads/${on_bin}/2_trimmomatic/cov_20_paired_1.fastq"
on2="${proj_dir}/reads/${on_bin}/2_trimmomatic/cov_20_paired_2.fastq"

off1="${proj_dir}/reads/${off_bin}/2_trimmomatic/cov_20_paired_1.fastq"
off2="${proj_dir}/reads/${off_bin}/2_trimmomatic/cov_20_paired_2.fastq"


if [ "$on_bin" == "$off_bin" ]; then
    echo "sim"
    cat $on1 | gzip -c > $r1
    cat $on2 | gzip -c > $r2
else
    echo "mngs"
    cat $on1 $off1 | gzip -c > $r1
    cat $on2 $off2 | gzip -c > $r2
fi


####################BT2_INDEX
bt2_dir="${proj_dir}/${on_bin}/3_bt2db/${off_bin}"
mkdir -p ${bt2_dir}

toc="${midasdb_dir}/genomes.tsv"
repspecies="${bt2_dir}/repgenome.species"
tail -n +2 ${toc} | awk '{print $2}' > ${repspecies}

python -m iggtools build_bowtie2db \
    --midasdb_name uhgg --midasdb_dir ${midasdb_dir} --num_cores 8 \
    --bt2_indexes_name repgenomes --species_list ${repspecies} \
    --bt2_indexes_dir ${bt2_dir}

######################MIDAS

midas_outdir="${proj_dir}/${on_bin}/4_midas/${off_bin}"
mkdir -p ${midas_outdir}
python3 -m iggtools midas_run_snps \
    --sample_name cov_20X -1 ${r1} -2 ${r2} \
    --midasdb_name uhgg --midasdb_dir ${midasdb_dir} \
    --prebuilt_bowtie2_indexes ${bt2_dir}/repgenomes \
    --prebuilt_bowtie2_species ${repspecies} \
    --advanced --analysis_ready --fragment_length 1000 --select_threshold=-1 \
    --num_cores 8 ${midas_outdir}


echo "${species}:${on_bin}:${off_bin} DONE"
