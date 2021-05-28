#!/usr/bin/env bash
#$ -cwd
#$ -l h_vmem=16G
#$ -l h_rt=0:30:00
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

REF=$(sed -n ${SGE_TASK_ID}p "$1")

DIRNAME=${REF%/*}
DIRNAME=${DIRNAME/\/gaemr\/work//}
OUTDIR=${DIRNAME}/mob_suite

set -euxo pipefail

# Redirect all stdout/stderr to a file
exec > "${DIRNAME}/mob_suite.log"
exec 2>&1

mob_recon --infile "${REF}" --outdir "${OUTDIR}" --force
