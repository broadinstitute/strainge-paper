import sys
import itertools
from pathlib import Path

import numpy
import pandas
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster


similarities = pandas.read_csv("straingst_kmersim.tsv", sep='\t',
                               index_col=['kmerset1', 'kmerset2']).sort_index()
similarities['ani'] = pandas.to_numeric(similarities['ani'])

species_repr = "Esch_coli_K-12_substr_MG1655_GCF_001544635.1"
ix = similarities.loc[species_repr, 'ani'] < 0.9
too_distant = similarities.loc[species_repr, 'ani'][ix].index

print("These strains are too distant from Esch_coli_K-12_substr_MG1655:",
      too_distant, file=sys.stderr)

similarities.drop(index=too_distant, level=1, inplace=True)
similarities.drop(index=too_distant, level=0, inplace=True)

labels = set(similarities.index.unique(level=0)).union(
             similarities.index.unique(level=1))
to_keep = set()
with open("straingst_like_db/reflist.txt") as f:
    for line in f:
        to_keep.add(Path(line).stem)
drop = labels - to_keep

similarities.drop(index=drop, level=1, inplace=True)
similarities.drop(index=drop, level=0, inplace=True)

similarities = similarities.reset_index().set_index(
    ['kmerset1', 'kmerset2']).sort_index()

labels = list(set(similarities.index.unique(level=0)).union(
    set(similarities.index.unique(level=1))))
labels.sort()
label_ix = {label: i for i, label in enumerate(labels)}

dmatrix = numpy.zeros((len(labels), len(labels)))

for kmerset1, kmerset2 in itertools.combinations(labels, 2):
    ix1 = label_ix[kmerset1]
    ix2 = label_ix[kmerset2]
    if (kmerset1, kmerset2) in similarities.index:
        dmatrix[ix1, ix2] = 1-similarities.loc[(kmerset1, kmerset2), 'ani']
        dmatrix[ix2, ix1] = 1-similarities.loc[(kmerset1, kmerset2), 'ani']
    else:
        dmatrix[ix1, ix2] = 1-similarities.loc[(kmerset2, kmerset1), 'ani']
        dmatrix[ix2, ix1] = 1-similarities.loc[(kmerset2, kmerset1), 'ani']

Z = linkage(squareform(dmatrix), method='complete')

clusters = fcluster(Z, 10, criterion='maxclust')

labels = numpy.array(labels)

for cluster in numpy.unique(clusters):
    strains_in_cluster = labels[clusters == cluster]

    cluster_strains = similarities.loc[(strains_in_cluster,
                                        strains_in_cluster), 'ani'].copy()
    if len(strains_in_cluster) == 1:
        print(strains_in_cluster[0])
    else:
        mean_dist = cluster_strains.groupby(level=1).mean()

        rep = mean_dist.idxmax()

        print(rep)
