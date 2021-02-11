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

mkdir -p "log/${sample_name}"

# Redirect all stdout/stderr to a file
exec > log/${sample_name}/straingr.log
exec 2>&1

TMPDIR=/broad/hptmp/ldijk

R1=${SAMPLE[0]}
R2=${SAMPLE[1]}
REF="straingr/${sample_name}/concat_refs.fa"
SORTED_BAM="straingr/${sample_name}/alignments.sorted.bam"
METRICS="straingr/${sample_name}/alignments.duplicate_metrics.txt"
BAM="straingr/${sample_name}/alignments.bam"

PREPARE_REF_TPL="/gsap/archive-bacterial/Projects/StrainGE/strainge-benchmarks/strainge_db/ecoli_chrom_db/{ref}.fa.gz"

mkdir -p "${BAM%/*}"

if [[ ! -f "${REF}" ]]; then
    orig_refs_file="mutated_refs/${sample_name%/*}/orig_refs.txt"
    refs=($(awk '{print $1}' < ${orig_refs_file}))

    straingr prepare-ref -r ${refs[@]} -p "${PREPARE_REF_TPL}" -o "${REF}"
    bwa index "${REF}"
fi

if [[ ! -f "${BAM}" ]]; then
    bwa mem -I 300 "${REF}" "${R1}" "${R2}" \
        | samtools sort  -O bam -o "${SORTED_BAM}"

    picard -Xmx30G MarkDuplicates \
        I="${SORTED_BAM}" O="${BAM}" M="${METRICS}"

    samtools index "${BAM}"
    rm "${SORTED_BAM}"
fi

HDF5=${BAM%/*}/straingr.hdf5
VCF=${HDF5%.hdf5}.vcf
SUMMARY=${HDF5%.hdf5}.tsv

straingr call "${REF}" "${BAM}" \
    --hdf5-out "${HDF5}" --vcf "${VCF}" --summary "${SUMMARY}" -t all
