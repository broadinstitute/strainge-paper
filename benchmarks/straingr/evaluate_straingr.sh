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
    echo "usage: ${0} sample_list.txt"
    exit 1
fi

SAMPLE=($(sed -n ${SGE_TASK_ID}p "$1"))

set -euxo pipefail

sample_name=${SAMPLE[0]%/*}
sample_name=${sample_name#samples/}
dirname=${SAMPLE[0]%/*}
fname=${SAMPLE[0]##*/}

mkdir -p "log/${sample_name}"

# Redirect all stdout/stderr to a file
exec > log/${sample_name}/evaluation.log
exec 2>&1

TMPDIR=/broad/hptmp/ldijk

REF="straingr/${sample_name}/concat_refs.fa"
VCF="straingr/${sample_name}/straingr.vcf"

sample_without_cov="${sample_name%/*}"
orig_refs_file="mutated_refs/${sample_without_cov}/orig_refs.txt"

# Evaluate each ref separately
while IFS= read -r line; do
    refdata=($line)
    refname=${refdata[0]}
    refpath=$(readlink -e "${refdata[1]}")

    # Get uncompressed path to ref
    refdir="${refpath%/*}"
    ref_fname="${refpath##*/}"
    uncompressed_path="${refdir}/uncompressed/${ref_fname%.gz}"

    if [[ ! -f "${uncompressed_path}" ]]; then
        mkdir -p "${refdir}/uncompressed"
        gunzip -c ${refpath} > ${uncompressed_path}
    fi

    truth_vcf="mutated_refs/${sample_without_cov}/${refname}.mutations.vcf"
    query_vcf="straingr/${sample_name}/straingr.vcf"

    ./run_sompy.py "${refname}" "${uncompressed_path}" "${truth_vcf}" "${query_vcf}"

done < ${orig_refs_file};

