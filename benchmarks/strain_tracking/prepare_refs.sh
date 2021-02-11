#!/usr/bin/env bash

sample_sets=($(find straingst/ -mindepth 3 -maxdepth 3 -type d))
for set in "${sample_sets[@]}"; do
  straingr_dir=straingr/${set#straingst/}
  mkdir -p "${straingr_dir}"

  straingr prepare-ref -S /gsap/garage-bacterial/Users/Bruce/newbench/db/similarities.tsv \
    -s "${set}"/*.tsv -p "../strainge_db/ecoli_chrom_db/{ref}.fa.gz" \
    -o "${straingr_dir}/concat_refs.fa"

  bwa index "${straingr_dir}/concat_refs.fa"
done
