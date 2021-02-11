#!/usr/bin/env python3
import sys
import random
import subprocess
from pathlib import Path

import pandas

NUM_SAMPLES = 20

random.seed(42)

references_meta = pandas.read_csv(
    '../strainge_db/ecoli_chrom_db/references_meta.tsv', sep='\t', comment='#',
    index_col=0, header=None, names=['ref_id', 'path', 'accesson'])

ref_outdir = Path("mutated_refs")

sample_types = {
    'single': 1,
    'equal2': 2,
    'equal3': 3,
    'equal4': 4,
    'unequal': 2
}

for stype, num_refs in sample_types.items():
    for i in range(NUM_SAMPLES):
        out_dir = ref_outdir / f"{stype}/{stype}{i+1}"
        out_dir.mkdir(parents=True, exist_ok=True)

        refs = random.sample(list(references_meta.index), num_refs)

        with (out_dir / "orig_refs.txt").open("w") as f:
            references_meta.loc[refs, 'path'].to_csv(f, sep='\t', header=False)

        if len(sys.argv) > 2 and sys.argv[1] == 'fix':
            continue

        for ref in refs:
            path = references_meta.loc[ref, 'path']

            out_prefix = out_dir / f"{ref}"
            subprocess.run([
                "mutate_ref.py", path, "-s", "99.9", "-o", out_prefix,
                "--ins-weight", "0", "--del-weight", "0"
            ])

        with (out_dir / "mutated_refs.txt").open("w") as f:
            for ref in refs:
                ref_path = out_dir / f"{ref}.fa"
                print(ref_path.resolve(), file=f)
