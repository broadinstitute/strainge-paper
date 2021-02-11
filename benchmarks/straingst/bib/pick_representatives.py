import subprocess
from pathlib import Path

import numpy
import pandas
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

STRAINEST_PATH = Path("../strainest_db")


with open(STRAINEST_PATH / "dist_matrix.phylip") as f:
    fiter = iter(f)
    num_entries = int(next(fiter))
    reference_names = []

    matrix = numpy.zeros((num_entries, num_entries))
    for i, line in enumerate(fiter):
        parts = line.strip().split('\t')
        reference_names.append(parts[0].replace(".sketch", ""))
        row = numpy.array([float(v) for v in parts[1:]])
        if len(row) > 0:
            matrix[i, 0:len(row)] = row

    index_upper = numpy.triu_indices(num_entries)
    matrix[index_upper] = matrix.T[index_upper]

dmatrix = pandas.DataFrame(matrix, index=reference_names,
                           columns=reference_names)
d = squareform(dmatrix.values)
Z = linkage(d, method='complete')

clusters = fcluster(Z, 20, criterion='maxclust')

for cluster in numpy.unique(clusters):
    ix = dmatrix.loc[clusters == cluster].index

    c_dmatrix = dmatrix.loc[ix, ix]
    mean_dist = c_dmatrix.mean(axis=1)

    representative = mean_dist.idxmin()
    print(representative)

    orig_path = Path(f"../../strainge_db/ecoli_db/{representative}.fa.gz")
    target = Path(f"rep_genomes/{representative}.fa")

    subprocess.run(f'gunzip -c "{orig_path}" > "{target}"', shell=True)
