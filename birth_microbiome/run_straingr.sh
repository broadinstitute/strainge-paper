#!/usr/bin/env bash
#$ -cwd
#$ -l h_vmem=32G
#$ -l h_rt=24:00:00
#$ -j y
#$ -R y
#$ -P gscid
#$ -e /dev/null
#$ -o /dev/null

source /broad/software/scripts/useuse
use Anaconda3

source activate strainge

if [ $# -lt 2 ]; then 
    echo "usage: ${0} sample_list.txt individual_list.txt"
    exit 1
fi

SAMPLE=($(sed -n ${SGE_TASK_ID}p "$1"))

set -euxo pipefail

sample_name=${SAMPLE[0]%_1.fastq.gz}
sample_name=${sample_name#data/}
dirname=${SAMPLE[0]%/*}
fname=${SAMPLE[0]##*/}
individual=$(sed -n ${SGE_TASK_ID}p "$2")

mkdir -p "log/${sample_name}"

# Redirect all stdout/stderr to a file
exec > "log/${sample_name}/straingr.log"
exec 2>&1

TMPDIR=/broad/hptmp/ldijk

R1=${SAMPLE[0]}
R2=${SAMPLE[1]}

REF="straingr/concat_refs/${individual}.fa"
BAM="straingr/${individual}/${sample_name}.bam"

mkdir -p "${BAM%/*}"

if [[ ! -f "${BAM}" ]]; then
    bwa mem -I 300 "${REF}" "${R1}" "${R2}" \
        | samtools sort -O bam -o "${BAM}"

    samtools index "${BAM}"
fi

HDF5=${BAM%.bam}.hdf5
VCF=${HDF5%.hdf5}.vcf
SUMMARY=${HDF5%.hdf5}.tsv

straingr call "${REF}" "${BAM}" \
    --hdf5-out "${HDF5}" --vcf "${VCF}" --summary "${SUMMARY}" -t all
