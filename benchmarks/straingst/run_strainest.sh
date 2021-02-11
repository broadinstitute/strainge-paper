#!/usr/bin/env bash
#$ -cwd
#$ -l h_vmem=16G
#$ -l h_rt=36:00:00
#$ -j y
#$ -R y
#$ -P gscid
#$ -e /dev/null
#$ -o /dev/null

source /broad/software/scripts/useuse
use Anaconda3

source activate strainest


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
exec > log/${sample_name}/strainest2.log
exec 2>&1

TMPDIR=/broad/hptmp/ldijk

STRAINST_DB_PREFIX=strainest_db/straingst_like_db/MA
STRAINST_DB_SNV=strainest_db/straingst_like_db/snv.txt

R1=${SAMPLE[0]}
R2=${SAMPLE[1]}

TRIMMED_R1=${R1%.1.fq.gz}.trimmed.1.fq
TRIMMED_R2=${R2%.2.fq.gz}.trimmed.2.fq
SINGLES=${R1%.1.fq.gz}.singles.fq

output_dir=strainest2/${sample_name}
mkdir -p ${output_dir}

START=$(date --rfc-3339=seconds)

if [[ ! -f "${output_dir}/abund.txt" ]]; then
    if [[ ! -f "${output_dir}/alignments.bam" ]]; then
        # Trim reads
        sickle pe -f ${R1} -r ${R2} -q 20 -t sanger \
            -o ${TRIMMED_R1} -p ${TRIMMED_R2} -s ${SINGLES}

        # Align reads to database
        bowtie2 --very-fast --no-unal -x ${STRAINST_DB_PREFIX} \
            -1 ${TRIMMED_R1} -2 ${TRIMMED_R2} \
            | samtools sort -O bam -o ${output_dir}/alignments.bam

        samtools index ${output_dir}/alignments.bam
    fi

    coverage=${dirname##*/}

    if [[ "${coverage}" == "10x" ]]; then
        min_depth=6
    else
        min_depth=1
    fi

    # Run StrainEst
    strainest est ${STRAINST_DB_SNV} ${output_dir}/alignments.bam ${output_dir} \
        --min-depth-absolute ${min_depth}

    rm ${TRIMMED_R1} ${TRIMMED_R2} ${SINGLES} ${output_dir}/alignments.bam
fi

END=$(date --rfc-3339=seconds)
