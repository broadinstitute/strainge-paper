from pathlib import Path
import pandas

accession_meta = pandas.read_csv("data/samples_ena.txt", sep='\t', index_col='secondary_sample_accession')
sample_meta = pandas.read_csv("data/samples_meta.txt", sep='\t', index_col='Accession')

accession_meta = accession_meta.join(sample_meta).reset_index().set_index('run_accession')


with open("metagenomes_all.txt") as f:
    for line in f:
        f1, f2 = line.strip().split()

        sample = Path(f1).name.replace("_1.fastq.gz", "")

        print(accession_meta.loc[sample, 'Individual'])
