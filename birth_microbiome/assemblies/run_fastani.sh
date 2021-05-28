#!/usr/bin/env bash
#$ -cwd
#$ -l h_vmem=32G
#$ -pe smp 16
#$ -j y
#$ -R y
#$ -q gscid
#$ -o fastani.log

source /broad/software/scripts/useuse
use Anaconda3

source activate strainge

fastANI --ql reflist_chrom.txt --rl reflist_chrom.txt -t 16 -o fastani_chrom.output --matrix
