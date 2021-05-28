#!/usr/bin/env bash
#$ -cwd
#$ -l h_vmem=16G
#$ -l h_rt=0:30:00
#$ -j y
#$ -R y
#$ -P gscid
#$ -o log/strainge_compare_complete.log

source /broad/software/scripts/useuse
use Anaconda3

source activate strainge

if [ $# -lt 1 ]; then 
    echo "usage: ${0} compare_list.txt"
    exit 1
fi

CMD=($(sed -n ${SGE_TASK_ID}p "$1"))
eval "${CMD[@]}"
