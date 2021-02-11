#!/usr/bin/env bash
#$ -cwd
#$ -l h_vmem=32G
#$ -l h_rt=16:00:00
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

mkdir -p log/straingr/${sample_name}

# Redirect all stdout/stderr to a file
exec > log/straingr/${sample_name}/${timepoint}.log
exec 2>&1

TMPDIR=/broad/hptmp/ldijk

R1=${SAMPLE[0]}
R2=${SAMPLE[1]}
REF="straingr/${sample_name}/concat_refs.fa"
SORTED_BAM="straingr/${sample_name}/${timepoint}.sorted.bam"
METRICS="straingr/${sample_name}/${timepoint}.duplicate_metrics.txt"
BAM="straingr/${sample_name}/${timepoint}.bam"

mkdir -p "${BAM%/*}"

if [[ ! -f "${BAM}" ]]; then
    bwa mem -I 300 "${REF}" "${R1}" "${R2}" \
        | samtools sort  -O bam -o "${SORTED_BAM}"

    picard -Xmx16G MarkDuplicates \
        I="${SORTED_BAM}" O="${BAM}" M="${METRICS}"

    samtools index "${BAM}"
    rm "${SORTED_BAM}"
fi

HDF5=${BAM%.bam}.hdf5
VCF=${BAM%.bam}.vcf
SUMMARY=${BAM%.bam}.tsv

straingr call "${REF}" "${BAM}" \
    --hdf5-out "${HDF5}" --vcf "${VCF}" --summary "${SUMMARY}"
