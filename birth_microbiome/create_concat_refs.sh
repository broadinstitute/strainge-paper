#!/usr/bin/env bash

for f in $(\ls -1 straingst/metagenomes/*.tsv); do
    sample=${f##*/}
    sample=${sample%.tsv}

    lineno=($(grep -n "${sample}" metagenomes_all.txt))
    lineno=${lineno[0]%%:*}
    individual=$(sed -n ${lineno}p "individual_list.txt")
    echo $individual;

    ref="straingr/concat_refs/${individual}.fa"

    # straingr prepare-ref -S /gsap/archive-bacterial/Projects/StrainGE/refseq/r_gnavus/db_${DB}/similarities.tsv \
    #     -p "/gsap/archive-bacterial/Projects/StrainGE/refseq/r_gnavus/db_${DB}/{ref}.fa" \
    #     -o straingr/${DB}/concat_refs/${patient}.fa \
    #     -s straingst/${DB}/metagenomes/${patient}_*.tsv 2> log/prepare_ref_${patient}_${DB}.log

    # bwa index straingr/${DB}/concat_refs/${patient}.fa
done
