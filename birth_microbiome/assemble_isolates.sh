#!/usr/bin/env bash
#$ -cwd
#$ -l os=RedHat7
#$ -l h_vmem=32G
#$ -pe smp 4
#$ -j y
#$ -R y
#$ -q gscid
#$ -e /dev/null
#$ -o /dev/null

source /broad/software/scripts/useuse
use Anaconda3

source activate strainge

if [ $# -lt 1 ]; then 
    echo "usage: ${0} isolate_list.txt"
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
exec > log/${sample_name}/assembly.log
exec 2>&1

export TMPDIR=/broad/hptmp/ldijk

R1=${SAMPLE[0]}
R2=${SAMPLE[1]}

trimmed_reads_pfx=data/isolates_trimmed/${sample_name#isolates/}
TRIMMED_R1="${trimmed_reads_pfx}.trimmed_1.fastq.gz"
TRIMMED_R2="${trimmed_reads_pfx}.trimmed_2.fastq.gz"
UNPAIRED_R1="${trimmed_reads_pfx}.unpaired_1.fastq.gz"
UNPAIRED_R2="${trimmed_reads_pfx}.unpaired_2.fastq.gz"

if [[ ! -s "${TRIMMED_R1}" ]]; then
    trimmomatic PE ${R1} ${R2} \
        "${TRIMMED_R1}" "${UNPAIRED_R1}" \
        "${TRIMMED_R2}" "${UNPAIRED_R2}" \
        ILLUMINACLIP:adapters_illumina_pe.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
fi

OUTDIR=assemblies/${sample_name}
mkdir -p "${OUTDIR}"

spades.py --isolate -1 "${TRIMMED_R1}" -2 "${TRIMMED_R2}" -t 4 -o "${OUTDIR}"
