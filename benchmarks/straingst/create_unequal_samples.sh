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

if [ $# -lt 1 ]; then
    echo "usage: ${0} ref_list.txt"
    exit 1
fi


REF_FILE=$(sed -n ${SGE_TASK_ID}p "$1")
sample_dir=${REF_FILE%/*}
sample_name=${sample_dir#samples/}

mkdir -p "log/${sample_name}"

exec >> "log/${sample_name}/create_sample.log"
exec 2>&1

set -x

rand_dir=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 8 | head -n 1)
TMPDIR="/broad/hptmp/ldijk/${rand_dir}"
mkdir -p "${TMPDIR}"

ref_paths=( $(awk '{ print $2; }' < "${REF_FILE}") )

if [[ ${#ref_paths[@]} == 2 ]]; then
    # First 10x:1x
    sample_out_dir="${sample_dir}/10x_1x"
    mkdir -p "${sample_out_dir}"
    tmp_dir="${TMPDIR}/10x_1x"
    mkdir -p "${tmp_dir}"

    # Simulate reads from each ref at fixed coverage
    ref_name1=${ref_paths[0]##*/}
    art_illumina -ss HS25 -na -i ${ref_paths[0]} -f 10 \
        -l 100 -m 200 -s 10 -o ${tmp_dir}/${ref_name1%fna}

    ref_name2=${ref_paths[1]##*/}
    art_illumina -ss HS25 -na -i ${ref_paths[1]} -f 1 \
        -l 100 -m 200 -s 10 -o ${tmp_dir}/${ref_name2%fna}

    cat ${tmp_dir}/*.1.fq > ${TMPDIR}/sample.1.fq
    cat ${tmp_dir}/*.2.fq > ${TMPDIR}/sample.2.fq

    gen_titration_sample.py \
        -R1 $TMPDIR/sample.1.fq -R2 $TMPDIR/sample.2.fq \
        -B1 SRS014613.1.fq.gz -B2 SRS014613.2.fq.gz \
        --depth 3000000000 -o "${sample_out_dir}/sample"

    # Then 1x:0.5x
    sample_out_dir="${sample_dir}/1x_0.5x"
    mkdir -p "${sample_out_dir}"
    tmp_dir="${TMPDIR}/1x_0.5x"
    mkdir -p "${tmp_dir}"

    # Simulate reads from each ref at fixed coverage
    ref_name1=${ref_paths[0]##*/}
    art_illumina -ss HS25 -na -i ${ref_paths[0]} -f 1 \
        -l 100 -m 200 -s 10 -o ${tmp_dir}/${ref_name1%fna}

    ref_name2=${ref_paths[1]##*/}
    art_illumina -ss HS25 -na -i ${ref_paths[1]} -f 0.5 \
        -l 100 -m 200 -s 10 -o ${tmp_dir}/${ref_name2%fna}

    cat ${tmp_dir}/*.1.fq > ${TMPDIR}/sample.1.fq
    cat ${tmp_dir}/*.2.fq > ${TMPDIR}/sample.2.fq

    gen_titration_sample.py \
        -R1 $TMPDIR/sample.1.fq -R2 $TMPDIR/sample.2.fq \
        -B1 SRS014613.1.fq.gz -B2 SRS014613.2.fq.gz \
        --depth 3000000000 -o "${sample_out_dir}/sample"

else
    # First 10x:1x:0.5x
    sample_out_dir="${sample_dir}/10x_1x_0.5x"
    mkdir -p "${sample_out_dir}"
    tmp_dir="${TMPDIR}/10x_1x_0.5x"
    mkdir -p "${tmp_dir}"

    # Simulate reads from each ref at fixed coverage
    ref_name1=${ref_paths[0]##*/}
    art_illumina -ss HS25 -na -i ${ref_paths[0]} -f 10 \
        -l 100 -m 200 -s 10 -o ${tmp_dir}/${ref_name1%fna}

    ref_name2=${ref_paths[1]##*/}
    art_illumina -ss HS25 -na -i ${ref_paths[1]} -f 1 \
        -l 100 -m 200 -s 10 -o ${tmp_dir}/${ref_name2%fna}

    ref_name3=${ref_paths[2]##*/}
    art_illumina -ss HS25 -na -i ${ref_paths[2]} -f 0.5 \
        -l 100 -m 200 -s 10 -o ${tmp_dir}/${ref_name3%fna}

    cat ${tmp_dir}/*.1.fq > ${TMPDIR}/sample.1.fq
    cat ${tmp_dir}/*.2.fq > ${TMPDIR}/sample.2.fq

    gen_titration_sample.py \
        -R1 $TMPDIR/sample.1.fq -R2 $TMPDIR/sample.2.fq \
        -B1 SRS014613.1.fq.gz -B2 SRS014613.2.fq.gz \
        --depth 3000000000 -o "${sample_out_dir}/sample"

    # Then 1x:0.5x:0.1x
    sample_out_dir="${sample_dir}/1x_0.5x_0.1x"
    mkdir -p "${sample_out_dir}"
    tmp_dir="${TMPDIR}/1x_0.5x_0.1x"
    mkdir -p "${tmp_dir}"

    # Simulate reads from each ref at fixed coverage
    ref_name1=${ref_paths[0]##*/}
    art_illumina -ss HS25 -na -i ${ref_paths[0]} -f 1 \
        -l 100 -m 200 -s 10 -o ${tmp_dir}/${ref_name1%fna}

    ref_name2=${ref_paths[1]##*/}
    art_illumina -ss HS25 -na -i ${ref_paths[1]} -f 0.5 \
        -l 100 -m 200 -s 10 -o ${tmp_dir}/${ref_name2%fna}

    ref_name3=${ref_paths[2]##*/}
    art_illumina -ss HS25 -na -i ${ref_paths[2]} -f 0.1 \
        -l 100 -m 200 -s 10 -o ${tmp_dir}/${ref_name3%fna}

    cat ${tmp_dir}/*.1.fq > ${TMPDIR}/sample.1.fq
    cat ${tmp_dir}/*.2.fq > ${TMPDIR}/sample.2.fq

    gen_titration_sample.py \
        -R1 $TMPDIR/sample.1.fq -R2 $TMPDIR/sample.2.fq \
        -B1 SRS014613.1.fq.gz -B2 SRS014613.2.fq.gz \
        --depth 3000000000 -o "${sample_out_dir}/sample"

fi
