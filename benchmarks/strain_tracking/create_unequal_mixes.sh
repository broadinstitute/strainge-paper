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

t=$SGE_TASK_ID

mkdir -p log/unequal/

exec >> log/unequal/timepoint${t}.log
exec 2>&1

set -x

rand_dir=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 8 | head -n 1)
TMPDIR="/broad/hptmp/ldijk/${rand_dir}"
mkdir -p "${TMPDIR}"

if [ $# -lt 2 ]; then
    echo "usage: ${0} ref_paths.txt out_dir"
    exit 1
fi

ref_paths=( $(cat ${1}) )

coverages=(0.1 0.5 1 10)
output_dir="${2}"
mkdir -p "${output_dir}"

case $t in 
    # A few single strain samples
    1)
        # Ref 1
        art_illumina -ss HS25 -na -i ${ref_paths[0]} -f 1 \
            -l 100 -m 200 -s 10 -o $TMPDIR/unequal.t${t}.

        echo -e "${ref_paths[0]}\t1" > "${output_dir}/unequal.t${t}.metadata"
        ;;
    2)
        # Ref 1
        art_illumina -ss HS25 -na -i ${ref_paths[0]} -f 10 \
            -l 100 -m 200 -s 10 -o $TMPDIR/unequal.t${t}.

        echo -e "${ref_paths[0]}\t10" > "${output_dir}/unequal.t${t}.metadata"
        ;;
    3)
        # Ref 2
        art_illumina -ss HS25 -na -i ${ref_paths[1]} -f 1 \
            -l 100 -m 200 -s 10 -o $TMPDIR/unequal.t${t}.

        echo -e "${ref_paths[1]}\t1" > "${output_dir}/unequal.t${t}.metadata"
        ;;
    4)
        # Ref 2
        art_illumina -ss HS25 -na -i ${ref_paths[1]} -f 10 \
            -l 100 -m 200 -s 10 -o $TMPDIR/unequal.t${t}.

        echo -e "${ref_paths[1]}\t10" > "${output_dir}/unequal.t${t}.metadata"
        ;;
    [5-9])
        rand_int=$(shuf -i 0-3 -n 1)
        cov1=${coverages[$rand_int]}
        rand_int=$(shuf -i 0-3 -n 1)
        cov2=${coverages[$rand_int]}
        # Both
        art_illumina -ss HS25 -na -i ${ref_paths[0]} -f ${cov1} \
            -l 100 -m 200 -s 10 -o $TMPDIR/unequal.t${t}.r1.
        art_illumina -ss HS25 -na -i ${ref_paths[1]} -f ${cov2} \
            -l 100 -m 200 -s 10 -o $TMPDIR/unequal.t${t}.r2.

        cat $TMPDIR/unequal.t${t}.r{1,2}.1.fq > $TMPDIR/unequal.t${t}.1.fq
        cat $TMPDIR/unequal.t${t}.r{1,2}.2.fq > $TMPDIR/unequal.t${t}.2.fq

        echo -e "${ref_paths[0]}\t${cov1}" > "${output_dir}/unequal.t${t}.metadata"
        echo -e "${ref_paths[1]}\t${cov2}" >> "${output_dir}/unequal.t${t}.metadata"
        ;;
    10)
        rand_int=$(shuf -i 0-3 -n 1)
        cov1=${coverages[$rand_int]}
        rand_int=$(shuf -i 0-3 -n 1)
        cov2=${coverages[$rand_int]}
        # Both
        art_illumina -ss HS25 -na -i ${ref_paths[0]} -f ${cov1} \
            -l 100 -m 200 -s 10 -o $TMPDIR/unequal.t${t}.r1.
        art_illumina -ss HS25 -na -i ${ref_paths[1]} -f ${cov2} \
            -l 100 -m 200 -s 10 -o $TMPDIR/unequal.t${t}.r2.

        cat $TMPDIR/unequal.t${t}.r{1,2}.1.fq > $TMPDIR/unequal.t${t}.1.fq
        cat $TMPDIR/unequal.t${t}.r{1,2}.2.fq > $TMPDIR/unequal.t${t}.2.fq

        echo -e "${ref_paths[0]}\t${cov1}" > "${output_dir}/unequal.t${t}.metadata"
        echo -e "${ref_paths[1]}\t${cov2}" >> "${output_dir}/unequal.t${t}.metadata"
        ;;
    *)
        echo "Invalid timepoint ${t}"
        exit 1
        ;;
esac

gen_titration_sample.py \
    -R1 $TMPDIR/unequal.t${t}.1.fq -R2 $TMPDIR/unequal.t${t}.2.fq \
    -B1 ../benchmark_data_gen/SRS014613.1.fq.gz -B2 ../benchmark_data_gen/SRS014613.2.fq.gz \
    --depth 3000000000 -o "${output_dir}/unequal.t${t}"
