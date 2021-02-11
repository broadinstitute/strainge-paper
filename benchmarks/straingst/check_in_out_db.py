from pathlib import Path
import pandas


in_db = set()
with open("/gsap/archive-bacterial/Projects/StrainGE/straingst-benchmarks/db/keep90.txt") as f:
    for line in f:
        in_db.add(Path(line.strip()).stem)

all_refs = []
for refs in Path("samples/").glob("**/refs.txt"):
    refs_df = pandas.read_csv(refs, sep='\t', names=['ref', 'path'])
    refs_df['sample_type'] = refs.parts[1]
    refs_df['in_db'] = refs_df['ref'].map(lambda r: r in in_db)

    all_refs.append(refs_df)

all_refs_df = pandas.concat(all_refs)

print(all_refs_df.groupby(['sample_type', 'in_db']).count())
print(all_refs_df.groupby('in_db').count())

