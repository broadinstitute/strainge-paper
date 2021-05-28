import subprocess
from pathlib import Path
from collections import defaultdict
import pandas

accession_meta = pandas.read_csv("data/samples_ena.txt", sep='\t', index_col='secondary_sample_accession')
sample_meta = pandas.read_csv("data/samples_meta.txt", sep='\t', index_col='Accession')

accession_meta = accession_meta.join(sample_meta).reset_index().set_index('run_accession')

samples_per_individual = defaultdict(list)

print("Collecting samples...")

for fname in Path("straingst/metagenomes").glob("*.tsv"):
    sample = fname.stem
    individual = accession_meta.loc[sample, 'Individual']

    samples_per_individual[individual].append(fname)

print("Done.")

for individual, samples in samples_per_individual.items():
    print("Building ref for", individual)
    refpath = f"straingr/concat_refs/{individual}.fa"

    with open(f"straingr/concat_refs/{individual}.log", "w") as f:
        p = subprocess.run([
            "straingr", "prepare-ref",
            "-S", "/gsap/archive-bacterial/Projects/StrainGE/refseq/enterococcus/strainge_db/similarities.tsv",
            "-p" "/gsap/archive-bacterial/Projects/StrainGE/refseq/enterococcus/strainge_db/{ref}.fa.gz",
            "-o", refpath,
            "-s", *samples
        ], text=True, capture_output=True)

        f.write(p.stdout)
        f.write("\n\nStd err:\n")
        f.write(p.stderr)
        f.write("\n")

        p = subprocess.run(["bwa", "index", refpath], text=True, capture_output=True)

        f.write(p.stdout)
        f.write("\n\nStd err:\n")
        f.write(p.stderr)
