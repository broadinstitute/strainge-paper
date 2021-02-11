#! /bin/bash

#$ -N straingr_compare
#$ -q gscid
#$ -cwd
#$ -pe smp 16
#$ -l h_vmem=256g
#$ -R y
#$ -e straingr_compare.err
#$ -o straingr_compare.out

export TMPDIR=/broad/hptmp/tstraub

source /broad/software/scripts/useuse
use Anaconda3
source activate strainge

~tstraub/bin/parallel -j 16 -a compare.sh
