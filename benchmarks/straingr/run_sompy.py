#!/usr/bin/env python3
"""
Evaluate the variant calls made by StrainGR

This script determines the contig names to only evaluate one particular
reference.
"""

import sys
import shlex
import subprocess
from pathlib import Path

import skbio

from strainge.io.utils import open_compressed

refname = sys.argv[1]
ref = sys.argv[2]
truth_vcf = sys.argv[3]
query_vcf = sys.argv[4]
callable_bed = str(Path(query_vcf).with_suffix(".callable.bed"))
output_prefix = str(Path(query_vcf).with_name(refname))

with open_compressed(ref) as f:
    contigs = [r.metadata['id'] for r in skbio.io.read(f, "fasta")]

location = ",".join(contigs)

subprocess.run(["samtools", "faidx", ref], check=True)

ref = shlex.quote(ref)
callable_bed = shlex.quote(callable_bed)
truth_vcf = shlex.quote(truth_vcf)
query_vcf = shlex.quote(query_vcf)
output_prefix = shlex.quote(output_prefix)

subprocess.run(f"""
source activate happy;

som.py -r {ref} -l {location} -T {callable_bed} \\
    -o {output_prefix} \\
    <(bcftools filter -i'TYPE="snp"' {truth_vcf}) \\
    <(bcftools filter -i'TYPE~"snp"' {query_vcf})
""", shell=True, check=True, executable="/bin/bash")
