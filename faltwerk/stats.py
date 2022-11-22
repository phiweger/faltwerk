from itertools import product
import random
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal

from hdbscan import HDBSCAN
import numpy as np

from esda import fdr
from esda.getisord import G_Local
from esda.moran import Moran_Local
from libpysal.weights import DistanceBand
import markov_clustering as mc
import networkx as nx

from faltwerk.geometry import get_alpha_carbon_atoms, euclidean_distance


# ------------------------------------------------------------------------------
# Spatial autocorrelation

def find_hotspots(model, features, method: Literal['getis_ord', 'moran'] = 'getis_ord', angstrom=8, false_discovery_rate=0.05, test_two_sided=False):
    '''
    Getis-Ord statistic for spatial association.

    "The analysis of Spatial Association by Use of Distance Statistics",
    Getis & Ord, Geographical Analysis, 1992

    Gi, other than Gi_star will exclude the position in the
    "center" of the current fn call (i != j).

    - https://pysal.org/esda/generated/esda.Moran_Local.html
    - https://pysal.org/esda/generated/esda.G_Local.html

    On calculation of p-values from permutations see pysal docs or look for
    "p_z_sim" in:

    - https://squidpy.readthedocs.io/en/latest/_modules/squidpy/gr/_ppatterns.html

    Usage:

    >>> from faltwerk.stats import find_hotspots
    >>> # Assuming the following residues have been found by some analysis, eg
    >>> # positive selection:
    >>> m = Fold('some.pdb')
    >>> res = [45, 47, 112, 123, 124]
    >>> selection = [1 if i in res else 0 for i in range(len(m))]
    >>> hotspots = find_hotspots(m, selection)

    Optional arguments:

    - method to test for local association -- "getis_ord" (default) or "moran"
    - angstrom -- maximum distance of what constitutes a residue's neighbor (default 8)
    - false_discovery_rate -- max. FDR to adjust for (default 0.05)
    - test_two_sided -- Test two-sided? (default False, ie one-sided)
    '''
    ca = list(get_alpha_carbon_atoms(model, only_coords=True))
    dist = DistanceBand(ca, angstrom, p=2, binary=True)
    # <star> .. include the present observation, see Getis and Ord, 1992

    # Cast to float otherwise G_Local() etc. cause problems when np is compiled
    # with numba to C bc/ it does not allow mixed types.
    features = np.array(features, dtype='float')

    if method == 'getis_ord':
        local = G_Local(features, dist, 'B', permutations=1000, star=True)
    
    elif method == 'moran':
        local = Moran_Local(features, dist, 'B', permutations=1000)

    else:
        raise ValueError('Method not implemented')

    if not test_two_sided:
        ps = local.p_z_sim
    else:
        ps = local.p_z_sim * 2

    FDR = fdr(ps, false_discovery_rate)
    # For Getis-Ord, we could use:
    # ps = local.p_norm
    # FDR = fdr(local.p_norm, false_discovery_rate)
    # ... but the resulting FDR p-values seem near identical
    return [1 if i < FDR else 0 for i in ps]


# ------------------------------------------------------------------------------
# Point pattern analysis
# https://geographicdata.science/book/notebooks/08_point_pattern_analysis.html#ripley-s-alphabet-of-functions

def cluster(fold, mask, method='HDBSCAN', angstrom=8, **kwargs):
    '''
    Instead of testing for local association, group residues based on density.

    Usage:

    >>> from faltwerk.stats import cluster
    >>> from faltwerk.geometry import get_alpha_carbon_atoms
    >>> cluster(model, mask, min_cluster_size=2)
    >>> cluster(model, hotspots)

    Optional arguments:

    - method -- "HDBSCAN" (default) or "MCL"
    - angstrom -- only for MCL clustering, what constitutes a neighbor (default 8)
    - kwargs -- only for HDBSCAN, passed to the clustering fn
    '''
    assert all([i in set([0, 1, False, True]) for i in set(mask)]), \
    'Ambiguous mask, only (0, 1) or boolean allowed'

    assert method in ['HDBSCAN', 'MCL'], f'Method {method} not implemented'

    points = list(get_alpha_carbon_atoms(fold, only_coords=True))
    X = [i for i, j in zip(points, mask) if j]
    
    if method == 'HDBSCAN':
        clusterer = HDBSCAN(**kwargs)

        yhat = clusterer.fit_predict(X)
        # HDBSCAN: -1 means no cluster assignment, others start w/ 0;
        # We want cluster labels > 0, and define 0 as no cluster assignment
        yhat = [i+1 if i > -1 else i for i in yhat] 
        yrev = list(yhat.copy()[::-1])

    elif method == 'MCL':
        yhat =  mcl(fold, mask, angstrom, **kwargs)
        yrev = yhat.copy()
        yrev.reverse()  # inplace operation

    l = []
    for i in mask:
        if i == 0:
            l.append(0)
        else:
            l.append(yrev.pop())
    return l


def mcl(fold, mask, angstrom=8, inflation=1.4, **kwargs):
    '''
    Helper fn for Markov chain clustering
    '''
    coords = list(enumerate(get_alpha_carbon_atoms(fold, only_coords=True)))
    candidates = {k: v for (k, v), m in zip(coords, mask) if m}
    
    edges = []
    for p1, p2 in product(candidates.keys(), candidates.keys()):
        if p1 != p2:
            dist = euclidean_distance(candidates[p1], candidates[p2])
            if dist <= angstrom:
                edges.append([p1, p2])

    G = nx.from_edgelist(edges)
    # G.edges

    matrix = nx.to_scipy_sparse_matrix(G)
    # Run MCL with default parameters
    result = mc.run_mcl(matrix, inflation=inflation)           
    clusters = mc.get_clusters(result)
    # print(clusters)
    # Returns list of tuples:
    # [(0, 1, 2, 3, ...), (8, 9, 10, 11, ...), (25, 26, 27, 28, ...)]

    # Vis?
    # mc.draw_graph(matrix, clusters, node_size=50, with_labels=True, edge_color="silver")

    cnt = 1  # cluster labels > 0, we define 0 as no cluster assignment
    yhat = []
    for c in clusters:
        for i in c:
            # Add cluster labels
            yhat.append(cnt)
        cnt += 1

    return yhat

