#!/usr/bin/env bash
#$ -cwd
#$ -l h_vmem=16G
#$ -l h_rt=4:00:00
#$ -j y
#$ -R y
#$ -P gscid
#$ -o download.log

source /broad/software/scripts/useuse

source activate strainge

if [ $# -lt 1 ]; then 
    echo "usage: ${0} download.cmds"
    exit 1
fi

set -x

CMD=($(sed -n ${SGE_TASK_ID}p "$1"))
eval "${CMD[@]}"
