#!/usr/bin/env bash
#$ -cwd
#$ -l h_vmem=16G
#$ -l h_rt=8:00:00
#$ -j y
#$ -R y
#$ -P gscid
#$ -o /dev/null
#$ -e /dev/null

source /broad/software/scripts/useuse
use Anaconda3

source activate midas
source midas_env.sh

if [ $# -lt 1 ]; then 
    echo "usage: ${0} sample_list.txt"
    exit 1
fi

SAMPLE=($(sed -n ${SGE_TASK_ID}p "$1"))

set -euxo pipefail

sample_name=${SAMPLE[0]%/*}
sample_name=${sample_name#samples/}
dirname=${SAMPLE[0]%/*}
fname=${SAMPLE[0]##*/}
timepoint=${fname%.1.fq.gz}

mkdir -p log/midas/${sample_name}

# Redirect all stdout/stderr to a file
exec > log/midas/${sample_name}/${timepoint}.log
exec 2>&1

TMPDIR=/broad/hptmp/ldijk

R1=${SAMPLE[0]}
R2=${SAMPLE[1]}

OUTPFX="midas/${sample_name}/${timepoint}"

run_midas.py species ${OUTPFX} \
    -1 ${SAMPLE[0]} \
    -2 ${SAMPLE[1]}

run_midas.py genes ${OUTPFX} \
    -1 ${SAMPLE[0]} \
    -2 ${SAMPLE[1]} \
    --species_id Escherichia_coli_58110 
    
run_midas.py snps ${OUTPFX} \
    -1 ${SAMPLE[0]} \
    -2 ${SAMPLE[1]} \
    --species_id Escherichia_coli_58110 
