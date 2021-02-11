#!/usr/bin/env bash
#$ -cwd
#$ -l h_vmem=32G
#$ -l h_rt=8:00:00
#$ -j y
#$ -R y
#$ -P gscid
#$ -e /dev/null
#$ -o /dev/null

source /broad/software/scripts/useuse
use Anaconda3

source activate strainphlan

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

mkdir -p log/strainphlan/${sample_name}

# Redirect all stdout/stderr to a file
exec > log/strainphlan/${sample_name}/${timepoint}.log
exec 2>&1

outdir="strainphlan/${sample_name}"
mkdir -p "${outdir}"

metaphlan2.py "${SAMPLE[0]},${SAMPLE[1]}" \
    "${outdir}/${timepoint}_profile.txt" \
    --bowtie2out "${outdir}/${timepoint}_bowtie2.txt" \
    --samout "${outdir}/${timepoint}.sam.bz2" \
    --input_type fastq \
    --nproc 2

sample2markers.py --ifn_samples "${outdir}/${timepoint}.sam.bz2" \
    --input_type sam --output_dir "${outdir}"
