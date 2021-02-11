import random
import subprocess
from pathlib import Path


NUM_IN_DB = 100
NUM_OUT_DB = 100
DB_DIR = Path("/gsap/archive-bacterial/Projects/StrainGE/strainge-benchmarks/"
              "strainge_db/ecoli_chrom_db/")


all_strains = {
    f.name.replace(".fa.gz", ""): f for f in DB_DIR.glob("*.fa.gz")
}

with (DB_DIR / "clusters.tsv").open() as f:
    in_db = {line.strip().split()[0] for line in f}

out_db = all_strains.keys() - in_db

in_db = list(in_db)
out_db = list(out_db)

for strain in random.sample(in_db, NUM_IN_DB):
    # Make sure ref is uncompressed for ART
    strain_path = all_strains[strain].resolve()
    uncompressed_path = (strain_path.parent / "uncompressed" /
                         strain_path.with_suffix("").name)  # remove .gz
    if not uncompressed_path.is_file():
        uncompressed_path.parent.mkdir(exist_ok=True, parents=True)

        subprocess.run(f'gunzip -c "{strain_path}" > "{uncompressed_path}"',
                       shell=True, check=True)

    print(strain, uncompressed_path, sep='\t')

for strain in random.sample(out_db, NUM_OUT_DB):
    # Make sure ref is uncompressed for ART
    strain_path = all_strains[strain].resolve()
    uncompressed_path = (strain_path.parent / "uncompressed" /
                         strain_path.with_suffix(""))  # remove .gz
    if not uncompressed_path.is_file():
        subprocess.run(f'gunzip -c "{strain_path}" > "{uncompressed_path}"',
                       shell=True, check=True)

    print(strain, uncompressed_path, sep='\t')
