#!/usr/bin/env bash
#$ -cwd
#$ -l h_vmem=16G
#$ -l h_rt=24:00:00
#$ -j y
#$ -R y
#$ -P gscid
#$ -e /dev/null
#$ -o /dev/null

source /broad/software/scripts/useuse
use Anaconda3

source activate midas
source midas_env.sh

set -euxo pipefail

mkdir -p "log/"

# Redirect all stdout/stderr to a file
exec >> "log/midas_strain_tracking.log"
exec 2>&1

TMPDIR=/broad/hptmp/ldijk

COVS=(0.1 0.5 1 10)

for cov in ${COVS[@]}; do
    strain_tracking.py id_markers --indir all_vs_all/midas/${cov}x/Escherichia_coli_58110 \
        --out all_vs_all/midas/${cov}x/Escherichia_coli_58110.markers \
        --samples $(cat midas_training_samples.txt) --min_reads 1 || true;

    strain_tracking.py track_markers --indir all_vs_all/midas/${cov}x/Escherichia_coli_58110 \
        --out all_vs_all/midas/${cov}x/Escherichia_coli_58110.marker_sharing \
        --markers all_vs_all/midas/${cov}x/Escherichia_coli_58110.markers --min_reads 1 || true;
done


