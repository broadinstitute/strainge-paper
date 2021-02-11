#!/usr/bin/env bash

for f in mutated_refs/simple/*/mutated_refs.txt; do
    sample_name=${f%/*}
    sample_name=${sample_name#mutated_refs/}

    outdir=samples/${sample_name}

    qsub -t 1-8 ./create_single.sh "${f}" "${outdir}"
done

# for f in mutated_refs/equal/*/mutated_refs.txt; do
#     sample_name=${f%/*}
#     sample_name=${sample_name#mutated_refs/}
# 
#     outdir=samples/${sample_name}
# 
#     qsub -t 6 ./create_equal_mixes.sh "${f}" "${outdir}"
# done
