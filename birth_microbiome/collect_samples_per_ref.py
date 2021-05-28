from pathlib import Path
from collections import defaultdict

import pandas

from strainge.io.utils import parse_straingst

accession_meta = pandas.read_csv("data/samples_ena.txt", sep='\t', index_col='secondary_sample_accession')
sample_meta = pandas.read_csv("data/samples_meta.txt", sep='\t', index_col='Accession')

accession_meta = accession_meta.join(sample_meta).reset_index().set_index('run_accession')

print("Checking StrainGST results...")

samples_per_ref = defaultdict(set)

for straingst_tsv in Path().glob("straingst/metagenomes/*.tsv"):
    sample = straingst_tsv.stem
    individual = accession_meta.loc[sample, 'Individual']

    collapsed_refs = pandas.read_csv(
        f"straingr/concat_refs/{individual}.collapsed.tsv", sep='\t', index_col=0)

    with straingst_tsv.open() as f:
        for r in parse_straingst(f):
            collapsed_ref = collapsed_refs.loc[r['strain'], 'post_cluster']
            samples_per_ref[collapsed_ref].add((individual, sample))

for ref, samples in samples_per_ref.items():
    sample_list_txt = Path(f"strainge_compare/{ref}/sample_list.txt")
    sample_list_txt.parent.mkdir(parents=True, exist_ok=True)

    with sample_list_txt.open('w') as f:
        for individual, sample in sorted(samples):
            fname = Path(f"straingr/{individual}/metagenomes/{sample}.hdf5")

            if fname.is_file():
                print(fname, file=f)
            else:
                print("No straingr results for", fname)
