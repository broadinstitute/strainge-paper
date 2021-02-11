#!/usr/bin/env bash
#$ -cwd
#$ -l h_vmem=16G
#$ -l h_rt=4:00:00
#$ -j y
#$ -R y
#$ -P gscid
#$ -e /dev/null
#$ -o /dev/null

source /broad/software/scripts/useuse
use Anaconda3

source activate bib


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
exec > log/${sample_name}/bib.log
exec 2>&1

TMPDIR=/broad/hptmp/ldijk
BITSEQ=./bib/BitSeq

R1=${SAMPLE[0]}
R2=${SAMPLE[1]}

START=$(date --rfc-3339=seconds)

# Trim reads
# Align reads to database
output_dir=bib/${sample_name}
mkdir -p ${output_dir}

CORE_FASTA_PREFIX=./bib/core_alignment_gapless
bowtie2 --very-fast --no-unal -a -x ${CORE_FASTA_PREFIX} \
    -1 ${R1} -2 ${R2} \
    | samtools sort -n -O bam -o ${output_dir}/alignments.bam

#samtools index ${output_dir}/alignments.bam

$BITSEQ/parseAlignment ${output_dir}/alignments.bam --trSeqFile ./bib/core_alignment_gapless.fasta \
    -o ${output_dir}/alignment_info.prob --trInfoFile ${output_dir}/genome_info.tr --uniform --verbose

$BITSEQ/estimateVBExpression -o ${output_dir}/abundance.txt ${output_dir}/alignment_info.prob -t ${output_dir}/genome_info.tr

END=$(date --rfc-3339=seconds)
