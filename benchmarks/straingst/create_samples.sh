#!/usr/bin/env bash
#$ -cwd
#$ -l h_vmem=16G
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
    echo "usage: ${0} ref_list.txt"
    exit 1
fi

OUTDIR="samples/pure_sim/10x/"
mkdir -p "${OUTDIR}"

REF=($(sed -n ${SGE_TASK_ID}p "$1"))

set -euxo pipefail

# Redirect all stdout/stderr to a file
mkdir -p "log/"
exec > log/sample${SGE_TASK_ID}.log
exec 2>&1

art_illumina -ss HS25 -na -i "${REF[1]}" -f 10 \
    -l 100 -m 200 -s 10 -o "${OUTDIR}/sample${SGE_TASK_ID}."
