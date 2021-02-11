#!/usr/bin/env python3
import random
import subprocess
from pathlib import Path

import pandas

NUM_SAMPLES = 10

random.seed(42)

references_meta = pandas.read_csv(
    '../strainge_db/ecoli_chrom_db/references_meta.tsv', sep='\t', comment='#',
    index_col=0, header=None, names=['ref_id', 'path', 'accesson'])

ref_outdir = Path("mutated_refs")


for stype in ['simple', 'equal', 'unequal']:
    for i in range(NUM_SAMPLES):
        out_dir = ref_outdir / f"{stype}/{stype}{i+1}"
        out_dir.mkdir(parents=True, exist_ok=True)

        refs = random.sample(list(references_meta.index), 2)

        with (out_dir / "orig_refs.txt").open("w") as f:
            print(references_meta.loc[refs, 'path'], sep='\n', file=f)

        for ref in refs:
            path = references_meta.loc[ref, 'path']

            # Create derived reference A
            out_prefix = out_dir / f"{ref}-a"
            subprocess.run([
                "mutate_ref.py", path, "-s", "99.9", "-o", out_prefix,
                "--ins-weight", "0", "--del-weight", "0"
            ])

            # Create derived reference B
            out_prefix = out_dir / f"{ref}-b"
            subprocess.run([
                "mutate_ref.py", path, "-s", "99.9", "-o", out_prefix,
                "--ins-weight", "0", "--del-weight", "0"
            ])

        with (out_dir / "mutated_refs.txt").open("w") as f:
            for ref in refs:
                path_a = out_dir / f"{ref}-a.fa"
                path_b = out_dir / f"{ref}-b.fa"
                print(path_a.resolve(), path_b.resolve(), sep='\n', file=f)
