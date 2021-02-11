#!/usr/bin/env bash
#$ -cwd
#$ -l h_vmem=64G
#$ -l h_rt=4:00:00
#$ -j y
#$ -R y
#$ -P gscid
#$ -e /dev/null
#$ -o /dev/null

source /broad/software/scripts/useuse
use Anaconda3

source activate strainge

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

mkdir -p log/straingst/${sample_name}

# Redirect all stdout/stderr to a file
exec > log/straingst/${sample_name}/${timepoint}.log
exec 2>&1

TMPDIR=/broad/hptmp/ldijk
STRAINGST_DB=/gsap/garage-bacterial/Users/Bruce/newbench/db/ecoli_chrom_db_90.hdf5

R1=${SAMPLE[0]}
R2=${SAMPLE[1]}

KMERSET="straingst/${sample_name}/${timepoint}.hdf5"
RESULTS_TSV="${KMERSET%.hdf5}.tsv"

mkdir -p "${KMERSET%/*}"

straingst kmerize -k 23 -o "${KMERSET}" "${R1}" "${R2}"
straingst run "${STRAINGST_DB}" "${KMERSET}" -o "${RESULTS_TSV}"

rm ${KMERSET}
