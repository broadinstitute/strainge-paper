#!/usr/bin/env python3
"""
This script will generate for each reference genome a BED file with features
spanning its entire contigs. This is used for evaluation with Illumina's
som.py, where we can specify this bed file to only evaluate those regions.

In this way, we can evaluate SNP calls for a single reference genome, when that
particular reference is part of a concatenated reference for StrainGR.
"""

from pathlib import Path

import skbio

from strainge.io.utils import open_compressed

REF_DIR = Path("mutated_refs")


def generate_contig_bed(fasta, output):
    with open_compressed(fasta) as f:
        for r in skbio.io.read(f, "fasta"):
            output.write("\t".join([r.metadata['id'], "0", str(len(r))])
                         + "\n")


for fname in REF_DIR.glob("**/orig_refs.txt"):
    with fname.open() as f:
        for line in f:
            if not line:
                continue

            ref_name, path = line.split()
            opath = fname.with_name(f"{ref_name}.contigs.bed")
            with opath.open('w') as o:
                generate_contig_bed(path, o)
