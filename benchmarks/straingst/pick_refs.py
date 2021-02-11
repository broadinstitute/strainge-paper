#!/usr/bin/env python3
import random
import subprocess
from pathlib import Path

import pandas

NUM_SAMPLES = 15

random.seed(41)

references_meta = pandas.read_csv(
    '../strainge_db/ecoli_chrom_db/references_meta.tsv', sep='\t', comment='#',
    index_col=0, header=None, names=['ref_id', 'path', 'accesson'])

ref_outdir = Path("samples")

sample_types = {
    'single': 1,
    'equal2': 2,
    'equal3': 3,
    'equal4': 4,
    'unequal2': 2,
    'unequal3': 3,
}

for stype, num_refs in sample_types.items():
    for i in range(NUM_SAMPLES):
        out_dir = ref_outdir / f"{stype}/sample{i+1}"
        out_dir.mkdir(parents=True, exist_ok=True)

        refs = random.sample(list(references_meta.index), num_refs)
        uncompressed_refs = []

        for ref in refs:
            path = Path(references_meta.loc[ref, 'path'])
            uncompressed_path = (path.parent / "uncompressed" /
                                 path.with_suffix("").name)  # remove .gz
            if not uncompressed_path.is_file():
                uncompressed_path.parent.mkdir(exist_ok=True, parents=True)

                subprocess.run(f'gunzip -c "{path}" > "{uncompressed_path}"',
                               shell=True, check=True)

            uncompressed_refs.append([ref, uncompressed_path])

        with (out_dir / "refs.txt").open("w") as f:
            for name, path in uncompressed_refs:
                print(name, path, sep='\t', file=f)
