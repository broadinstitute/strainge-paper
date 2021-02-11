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
exec >> "log/midas_all_vs_all.log"
exec 2>&1

TMPDIR=/broad/hptmp/ldijk

COVS=(0.1 0.5 1 10)

for cov in ${COVS[@]}; do

    outdir="all_vs_all/midas/${cov}x"
    sample_dir="${outdir}/samples"
    mkdir -p "${sample_dir}"

    for f in midas/equal/equal*/${cov}x/*; do
        sample_set=${f#midas/equal/}
        sample_set=${sample_set%%/*}

        fname=${f##*/}
        fname=${fname/equal/${sample_set}}

        ln -s $(readlink -e $f) ${sample_dir}/${fname}
    done

    merge_midas.py species -i ${sample_dir} -t dir "${outdir}"
    merge_midas.py genes -i ${sample_dir} -t dir "${outdir}" \
        --species_id Escherichia_coli_58110
    merge_midas.py snps -i ${sample_dir} -t dir "${outdir}" --all_samples \
        --species_id Escherichia_coli_58110

done
