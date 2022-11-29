#!/bin/bash
set -e
set -x

if [ $# -ne 1 ]; then
    echo "Usage: $0 SPECIES"
    exit 1
fi

species="$1"

proj_dir="/mnt/cz/20220602_close_species/${species}"

ls $proj_dir | grep "^ani" | grep -v bamcounts | xargs -Ixx bash -c \
  "ls $proj_dir/xx/1_dbs | xargs -Iyy -P 8 bash -c 'bash /mnt/cz/cross_mapping_simulation/run_exp2.sh ${species} xx yy' "

## bwa
ls $proj_dir | grep "^ani" | grep -v bamcounts | xargs -Ixx bash -c \
  "ls $proj_dir/xx/1_dbs | xargs -Iyy -P 8 bash -c 'bash /mnt/cz/cross_mapping_simulation/run_exp2_bwa.sh ${species} xx yy' "

## exp1
ls $proj_dir | grep "^ani" | grep -v bamcounts | xargs -Ixx bash -c \
  "ls $proj_dir/xx/1_dbs | xargs -Iyy -P 4 bash -c 'bash /mnt/cz/cross_mapping_simulation/run_exp1.sh ${species} xx yy' "
