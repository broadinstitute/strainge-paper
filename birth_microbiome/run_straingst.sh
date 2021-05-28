#!/usr/bin/env bash
#$ -cwd
#$ -l h_vmem=32G
#$ -l h_rt=2:00:00
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

sample_name=${SAMPLE[0]%_1.fastq.gz}
sample_name=${sample_name#data/}
dirname=${SAMPLE[0]%/*}
fname=${SAMPLE[0]##*/}

mkdir -p "log/${sample_name}"

# Redirect all stdout/stderr to a file
exec > log/${sample_name}/straingst.log
exec 2>&1

TMPDIR=/broad/hptmp/ldijk
STRAINGST_DB=/gsap/archive-bacterial/Projects/StrainGE/refseq/enterococcus/strainge_db/enterococcus.hdf5

R1=${SAMPLE[0]}
R2=${SAMPLE[1]}

KMERSET="straingst/kmersets/${sample_name}.hdf5"
RESULTS_TSV="straingst/${sample_name}.tsv"

mkdir -p "${KMERSET%/*}"
mkdir -p "${RESULTS_TSV%/*}"

if [[ ! -s "${KMERSET}" ]]; then
    straingst kmerize -k 23 -o "${KMERSET}" "${R1}" "${R2}"
fi

straingst run "${STRAINGST_DB}" "${KMERSET}" -o "${RESULTS_TSV}"
