#!/usr/bin/env bash
#$ -cwd
#$ -l h_vmem=16G
#$ -l h_rt=2:00:00
#$ -j y
#$ -R y
#$ -P gscid
#$ -e /dev/null
#$ -o /dev/null

source /broad/software/scripts/useuse
use Anaconda3

source activate strainphlan

set -uxo pipefail

mkdir -p "log/"

# Redirect all stdout/stderr to a file
exec > "log/strainphlan_all_vs_all.log"
exec 2>&1

TMPDIR=/broad/hptmp/ldijk

COVS=(0.1 0.5 1 10)

for cov in ${COVS[@]}; do
    outdir="all_vs_all/strainphlan/${cov}x"
    sample_dir="${outdir}/samples"
    mkdir -p "${sample_dir}"

    for f in strainphlan/equal/equal*/${cov}x/*.markers; do
        sample_set=${f#strainphlan/equal/}
        sample_set=${sample_set%%/*}

        fname=${f##*/}
        fname=${fname/equal/${sample_set}}

        ln -s $(readlink -e $f) ${sample_dir}/${fname}
    done

    # Relaxed parameters
    strainphlan.py --ifn_samples ${sample_dir}/*.markers --output_dir "${outdir}" \
        --mpa_pkl /cil/shed/apps/internal/conda/envs/strainphlan/bin/metaphlan_databases/mpa_v295_CHOCOPhlAn_201901.pkl \
        --print_clades_only > "${outdir}/clades.txt"

    strainphlan.py --ifn_markers "strainphlan/db_markers/s__Escherichia_coli.markers.fasta" \
        --ifn_samples ${sample_dir}/*.markers --clades s__Escherichia_coli \
        --mpa_pkl /cil/shed/apps/internal/conda/envs/strainphlan/bin/metaphlan_databases/mpa_v295_CHOCOPhlAn_201901.pkl \
        --output_dir "${outdir}"

    distmat -nucmethod 2 -outfile "${outdir}/s__Escherichia_coli.distmat" "${outdir}/s__Escherichia_coli.fasta" || true;
done
