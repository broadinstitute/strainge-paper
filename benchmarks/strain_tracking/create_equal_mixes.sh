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

if [ $# -lt 2 ]; then
    echo "usage: ${0} ref_paths.txt out_dir"
    exit 1
fi

t=$SGE_TASK_ID
sample_name=${1%/*}
sample_name=${sample_name#mutated_refs/}

mkdir -p "log/${sample_name}"

exec >> "log/${sample_name}/timepoint${t}.log"
exec 2>&1

set -x

rand_dir=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 8 | head -n 1)
TMPDIR="/broad/hptmp/ldijk/${rand_dir}"
mkdir -p "${TMPDIR}"

echo "Generating mixes using refs $1, output dir: $2"

ref_paths=( $(cat ${1}) )

if [[ ${t} -gt 6 ]]; then
    base_ref_ix=1
    t_alt=$((t - 6))
else
    base_ref_ix=0
    t_alt=$t
fi

ref1=$base_ref_ix
ref2=$(( base_ref_ix + 2))

coverages=(0.1 0.5 1 10)
for cov in ${coverages[@]}; do
    output_dir="${2}/${cov}x"
    mkdir -p "${output_dir}"

    case $t_alt in 
        [1-2])
            # Ref 1
            art_illumina -ss HS25 -na -i ${ref_paths[$ref1]} -f ${cov} \
                -l 100 -m 200 -s 10 -o $TMPDIR/equal.t${t}.
            ;;
        [3-4])
            # Ref 2
            art_illumina -ss HS25 -na -i ${ref_paths[$ref2]} -f ${cov} \
                -l 100 -m 200 -s 10 -o $TMPDIR/equal.t${t}.
            ;;
        [5-6])
            # Both
            art_illumina -ss HS25 -na -i ${ref_paths[$ref1]} -f ${cov} \
                -l 100 -m 200 -s 10 -o $TMPDIR/equal.t${t}.r1.
            art_illumina -ss HS25 -na -i ${ref_paths[$ref2]} -f ${cov} \
                -l 100 -m 200 -s 10 -o $TMPDIR/equal.t${t}.r2.

            cat $TMPDIR/equal.t${t}.r{1,2}.1.fq > $TMPDIR/equal.t${t}.1.fq
            cat $TMPDIR/equal.t${t}.r{1,2}.2.fq > $TMPDIR/equal.t${t}.2.fq
            ;;
        *)
            echo "Invalid timepoint ${t_alt}"
            exit 1
            ;;
    esac

    gen_titration_sample.py \
        -R1 $TMPDIR/equal.t${t}.1.fq -R2 $TMPDIR/equal.t${t}.2.fq \
        -B1 ../benchmark_data_gen/SRS014613.1.fq.gz -B2 ../benchmark_data_gen/SRS014613.2.fq.gz \
        --depth 3000000000 -o "${output_dir}/equal.t${t}"
done
