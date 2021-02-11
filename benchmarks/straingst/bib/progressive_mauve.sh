#!/usr/bin/env bash
#$ -cwd
#$ -l h_vmem=64G
#$ -l h_rt=48:00:00
#$ -j y
#$ -R y
#$ -P gscid
#$ -e /dev/null
#$ -o mauve.log

source /broad/software/scripts/useuse
use Anaconda3

source activate bib

export TMPDIR=/broad/hptmp/ldijk

progressiveMauve --output=full_alignment.xmfa rep_genomes/*.fa
