import json
import subprocess
from pathlib import Path
from collections import defaultdict

from strainge.io.utils import parse_straingst

STRAINGR_DIR = Path('straingr/equal')
STRAINGST_DIR = Path('straingst/equal')
REF_DIR = Path('/gsap/archive-bacterial/Projects/StrainGE/strainge-benchmarks/'
               'strainge_db/ecoli_chrom_db/')
OUTPUT_DIR = Path("all_vs_all/strainge/")

samples_per_ref = defaultdict(lambda: defaultdict(set))

for sample_set_dir in sorted(STRAINGST_DIR.glob('**/*x')):
    if not sample_set_dir.is_dir():
        continue

    print("Collecting samples from", sample_set_dir)

    sample_set = sample_set_dir.parts[-2]
    cov = sample_set_dir.parts[-1]
    straingr_dir = STRAINGR_DIR / sample_set / cov

    for sample in sample_set_dir.glob("equal*.tsv"):
        orig_straingr_hdf5 = straingr_dir / (sample.stem + ".hdf5")
        new_hdf5 = OUTPUT_DIR / cov / "samples" / (sample.stem.replace("equal", sample_set) + ".hdf5")
        new_hdf5.parent.mkdir(exist_ok=True, parents=True)
        if not new_hdf5.is_file():
            new_hdf5.symlink_to(orig_straingr_hdf5.resolve())

        with sample.open() as f:
            for r in parse_straingst(f):
                samples_per_ref[cov][r['strain']].add(new_hdf5)

for cov, strains in samples_per_ref.items():
    output_dir = OUTPUT_DIR / cov

    for ref, samples in strains.items():
        print(ref, samples)
        ref_path = REF_DIR / (ref + ".fa.gz")
        output = output_dir / (ref + ".dist.tsv")

        if output.is_file():
            print(output, "already exists, skipping.")
            continue

        p = subprocess.run(['straingr', 'dist', '-r', ref_path, "-d", "jc", "-c", "0.5",
                            "-o", output, *samples])
        p.check_returncode()
