from typing import Literal

from hdbscan import HDBSCAN
import numpy as np

from esda import fdr
from esda.getisord import G_Local
from esda.moran import Moran_Local
from libpysal.weights import DistanceBand

from faltwerk.geometry import get_alpha_carbon_atoms 


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
    '''
    ca = list(get_alpha_carbon_atoms(model, only_coords=True))
    dist = DistanceBand(ca, angstrom, p=2, binary=True)
    # <star> .. include the present observation, see Getis and Ord, 1992
    
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


# TODO: MCL
def cluster(fold, mask, *args, **kwargs):
    '''
    from faltwerk.stats import cluster
    from faltwerk.geometry import get_alpha_carbon_atoms

    mask = [1 if i < 0.05 else 0 for i in d['meme']['positive']['scores']]
    cluster(model, mask, min_cluster_size=2)
    '''
    points = list(get_alpha_carbon_atoms(fold, only_coords=True))
    X = [i for i, j in zip(points, mask) if j]
    clusterer = HDBSCAN(*args, **kwargs)
    return clusterer.fit_predict(X)



# 

'''
TODO:

random forest, predict whether a residue is part of a "positive selection" hot spot from its features, then rank variable importance to see which features 
matter most. This assumes that a single objective acted on the protein, which
is unlikely but let's try this anyway.

==?

enrichment: given spatial clusters, are they enriched in any feature?

- distance to ligand
- solvent accessibility
- interface


find spatial clusters > segment with HDBSCAN > for each cluster, predict whether a residue is in or out, the rank features by importance to find most
likely explanatory factor for this cluster. Stop at any step (clusters too small, ...) -- or logistic regression (Bayes); something simple -- https://stats.stackexchange.com/questions/192310/is-random-forest-suitable-for-very-small-data-sets
'''


'''
TODO: We could feed the clusters to the Getis-Ord statistic or calculate eg
the (adjusted) Rand score.

Or feed the significant areas to itself now using other features (interface,
...)
'''





