#!/usr/bin/env python3
import sys
import random
import subprocess
from pathlib import Path

import pandas

NUM_SAMPLES = 20

random.seed(43)

references_meta = pandas.read_csv(
    '../strainge_db/ecoli_db/references_meta.tsv', sep='\t', comment='#',
    index_col=0, header=None, names=['ref_id', 'path', 'accesson'])

ref_outdir = Path("mutated_refs/gene_del")


for i in range(NUM_SAMPLES):
    out_dir = ref_outdir / f"ref{i+1}"
    out_dir.mkdir(parents=True, exist_ok=True)

    ref = random.choice(list(references_meta.index))

    fasta = Path(references_meta.loc[ref, 'path']).resolve()
    if fasta.suffixes[-1] == '.gz':
        gff = fasta.with_suffix("").with_suffix(".gff.gz")
    else:
        gff = fasta.with_suffix(".gff.gz")

    out_prefix = str(out_dir / f"{ref}")
    seed = random.randrange(2**32 - 1)
    print(f"ref{i+1}", "Seed for mutations:", seed, file=sys.stderr)
    subprocess.run([
        "mutate_ref.py", str(fasta), str(gff), "-s", "100", "-o", out_prefix,
        "-d", str(7.5), "--seed", str(seed), "-m", "5000"
    ])

    orig_ref_link = out_dir / ("orig_ref" + "".join(fasta.suffixes))
    orig_ref_link.symlink_to(fasta)
