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

source activate strainphlan

set -uxo pipefail

if [ $# -lt 1 ]; then 
    echo "usage: ${0} subdir [orig_refs]"
    exit 1
fi

mkdir -p "log/$1"

# Redirect all stdout/stderr to a file
exec > "log/$1/strainphlan_collect.log"
exec 2>&1

TMPDIR=/broad/hptmp/ldijk

if [[ $# -gt 1 ]]; then
    orig_refs=($(cat "$2"))
    ifn_ref_genomes="--ifn_ref_genomes ${orig_refs[@]}"
else
    ifn_ref_genomes=""
fi

indir="$1"

outdir1="${indir}/default"
outdir2="${indir}/relaxed"

mkdir -p "${outdir1}"
mkdir -p "${outdir2}"

strainphlan.py --ifn_samples ${indir}/*.markers --output_dir "${outdir1}" \
    --mpa_pkl /cil/shed/apps/internal/conda/envs/strainphlan/bin/metaphlan_databases/mpa_v295_CHOCOPhlAn_201901.pkl \
    --print_clades_only > "${outdir1}/clades.txt"

strainphlan.py --ifn_markers "strainphlan/db_markers/s__Escherichia_coli.markers.fasta" \
    --ifn_samples ${indir}/*.markers --clades s__Escherichia_coli \
    ${ifn_ref_genomes} \
    --mpa_pkl /cil/shed/apps/internal/conda/envs/strainphlan/bin/metaphlan_databases/mpa_v295_CHOCOPhlAn_201901.pkl \
    --output_dir "${outdir1}"

# Don't exit script on fail
distmat -nucmethod 2 -outfile "${outdir1}/s__Escherichia_coli.distmat" "${outdir1}/s__Escherichia_coli.fasta" || true

# Relaxed parameters
strainphlan.py --ifn_samples ${indir}/*.markers --output_dir "${outdir2}" \
    --mpa_pkl /cil/shed/apps/internal/conda/envs/strainphlan/bin/metaphlan_databases/mpa_v295_CHOCOPhlAn_201901.pkl \
    --relaxed_parameters \
    --print_clades_only > "${outdir2}/clades.txt"

strainphlan.py --ifn_markers "strainphlan/db_markers/s__Escherichia_coli.markers.fasta" \
    --ifn_samples ${indir}/*.markers --clades s__Escherichia_coli \
    ${ifn_ref_genomes} \
    --mpa_pkl /cil/shed/apps/internal/conda/envs/strainphlan/bin/metaphlan_databases/mpa_v295_CHOCOPhlAn_201901.pkl \
    --relaxed_parameters \
    --output_dir "${outdir2}"

distmat -nucmethod 2 -outfile "${outdir2}/s__Escherichia_coli.distmat" "${outdir2}/s__Escherichia_coli.fasta"
