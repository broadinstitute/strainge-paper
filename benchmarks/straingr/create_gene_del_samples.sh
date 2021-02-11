#!/usr/bin/env bash
#$ -cwd
#$ -l h_vmem=16G
#$ -l h_rt=12:00:00
#$ -j y
#$ -R y
#$ -P gscid
#$ -e /dev/null
#$ -o /dev/null

source /broad/software/scripts/useuse
use Anaconda3

source activate strainge

mutated_ref=("mutated_refs/gene_del/ref${SGE_TASK_ID}/"*.fa)
sample_name="gene_del/ref${SGE_TASK_ID}"
OUTPUT_DIR="samples/${sample_name}"

mkdir -p "log/${sample_name}"

exec > "log/${sample_name}/create_sample.log"
exec 2>&1

set -x
shopt -s extglob

rand_dir=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 8 | head -n 1)
TMPDIR="/broad/hptmp/ldijk/${rand_dir}"
mkdir -p "${TMPDIR}"

coverages=(0.1 0.5 1 10)
for cov in ${coverages[@]}; do
    sample_out_dir="${OUTPUT_DIR}/${cov}x"
    mkdir -p "${sample_out_dir}"

    tmp_dir="${TMPDIR}/${cov}x"
    mkdir -p "${tmp_dir}"

    # Simulate reads from each ref at fixed coverage
    ref_name=${mutated_ref[0]##*/}
    art_illumina -ss HS25 -na -i ${mutated_ref} -f ${cov} \
        -l 100 -m 200 -s 10 -o ${tmp_dir}/${ref_name%fa}

    gen_titration_sample.py \
        -R1 ${tmp_dir}/*.1.fq -R2 ${tmp_dir}/*.2.fq \
        -B1 SRS014613.1.fq.gz -B2 SRS014613.2.fq.gz \
        --depth 3000000000 -o "${sample_out_dir}/sample"
done

