import math
import itertools

import numpy
import pandas
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from sklearn.metrics import silhouette_score

EQUAL_A_SEED = 'equal.t1'
EQUAL_B_SEED = 'equal.t3'
EQUAL_MIX_SEED = 'equal.t5'

TRUTH_EQUAL_SMALL_SET = list(f'equal.t{i}' for i in range(1, 7))
TRUTH_EQUAL_COMPLETE_SET = list(f'equal.t{i}' for i in range(1, 13))
TRUTH_EQUAL_A = {f'equal.t{i}' for i in [1, 2, 5, 6]}
TRUTH_EQUAL_B = {f'equal.t{i}' for i in [3, 4, 5, 6]}
TRUTH_EQUAL_SAME_A = {f'equal.t{i}' for i in [1, 2]}
TRUTH_EQUAL_SAME_B = {f'equal.t{i}' for i in [3, 4]}
TRUTH_EQUAL_SAME_MIX = {f'equal.t{i}' for i in [5, 6]}

SINGLE_STRAIN_SAMPLES = {f'equal.t{i}' for i in [1, 2, 3, 4, 7, 8, 9, 10]}
TRUTH_SNPSET1_SINGLE_A = {f'equal.t{i}' for i in [1, 2]}
TRUTH_SNPSET2_SINGLE_A = {f'equal.t{i}' for i in [7, 8]}
TRUTH_SNPSET1_SINGLE_B = {f'equal.t{i}' for i in [3, 4]}
TRUTH_SNPSET2_SINGLE_B = {f'equal.t{i}' for i in [9, 10]}

TRUTH_SNPSET1_ALL_A = {f'equal.t{i}' for i in [1, 2, 5, 6]}
TRUTH_SNPSET1_ALL_B = {f'equal.t{i}' for i in [3, 4, 5, 6]}
TRUTH_SNPSET2_ALL_A = {f'equal.t{i}' for i in [7, 8, 11, 12]}
TRUTH_SNPSET2_ALL_B = {f'equal.t{i}' for i in [9, 10, 11, 12]}

SNPSET1_A_SEED = 'equal.t1'
SNPSET1_B_SEED = 'equal.t3'
SNPSET2_A_SEED = 'equal.t7'
SNPSET2_B_SEED = 'equal.t9'

SET_LIST = range(1, 11)
COVERAGES = [0.1, 0.5, 1, 10]

PER_SAMPLE_TRUTH = {}

for s in range(1, 11):
    PER_SAMPLE_TRUTH.update({
        f'equal{s}.t1': {f'equal{s}.t{i}' for i in [1, 2, 5, 6]},
        f'equal{s}.t2': {f'equal{s}.t{i}' for i in [1, 2, 5, 6]},
        f'equal{s}.t3': {f'equal{s}.t{i}' for i in [3, 4, 5, 6]},
        f'equal{s}.t4': {f'equal{s}.t{i}' for i in [3, 4, 5, 6]},
        f'equal{s}.t5': {f'equal{s}.t{i}' for i in [1, 2, 3, 4, 5, 6]},
        f'equal{s}.t6': {f'equal{s}.t{i}' for i in [1, 2, 3, 4, 5, 6]},
        f'equal{s}.t7': {f'equal{s}.t{i}' for i in [7, 8, 11, 12]},
        f'equal{s}.t8': {f'equal{s}.t{i}' for i in [7, 8, 11, 12]},
        f'equal{s}.t9': {f'equal{s}.t{i}' for i in [9, 10, 11, 12]},
        f'equal{s}.t10': {f'equal{s}.t{i}' for i in [9, 10, 11, 12]},
        f'equal{s}.t11': {f'equal{s}.t{i}' for i in [7, 8, 9, 10, 11, 12]},
        f'equal{s}.t12': {f'equal{s}.t{i}' for i in [7, 8, 9, 10, 11, 12]},
    })


SAMPLE_PAIRS = list(itertools.combinations(PER_SAMPLE_TRUTH.keys(), 2))
TRUTH_LABELS = numpy.array([s2 in PER_SAMPLE_TRUTH[s1] for s1, s2 in SAMPLE_PAIRS])


def mcc(tp, fp, tn, fn):
    if min([tp+fp, tp+fn, tn+fp, tn+fn]) == 0:
        return (tp * tn) - (fp * fn)
    else:
        return (((tp * tn) - (fp * fn)) /
                math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))


def automatic_cluster(dmatrix: pandas.DataFrame, max_clusters=10):
    silhouette_scores = []
    assigned_clusters = []
    ks = range(2, max_clusters + 1)

    # Special case when we only have two elements. Clustering doesn't make
    # sense in that case. Only assign the same cluster to both elements if
    # the distance is 0.
    if len(dmatrix) == 2:
        if dmatrix.iloc[0, 1] == 0.0:
            return 1, numpy.array([1, 1]), None
        else:
            return 2, numpy.array([1, 2]), None

    for k in ks:
        Z = linkage(squareform(dmatrix), method='average')
        clusters = fcluster(Z, k, criterion='maxclust')

        num_clusters = len(numpy.unique(clusters))
        if num_clusters == 1:
            # Only one cluster created, even though we set maxclust > 1 in
            # fcluster. Force the point the farthest away from all other
            # points to a separate cluster.
            farthest = dmatrix.mean(axis=1).idxmax()
            ix = list(dmatrix.index).index(farthest)
            clusters[ix] = 2

        silhouette = silhouette_score(dmatrix, clusters, 'precomputed')
        silhouette_scores.append(silhouette)
        assigned_clusters.append((clusters, Z))

    ix = numpy.argmax(silhouette_scores)

    return (ks[ix], *assigned_clusters[ix])


def calculate_mcc(cluster: set, truth: set, all_samples: set):

    tp = len(truth & cluster)
    fp = len(cluster - truth)
    fn = len(truth - cluster)
    tn = len((all_samples - truth) - cluster)

    return tp, fp, fn, tn, mcc(tp, fp, tn, fn)
