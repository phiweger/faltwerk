from pathlib import Path

from libpysal.weights import DistanceBand

from faltwerk.models import Fold
from faltwerk.geometry import get_alpha_carbon_atoms, is_close


def test_distance_band():
    '''
    https://stackoverflow.com/questions/32527861/python-unit-test-that-uses-an-external-data-file
    '''
    here = Path(__file__).parent
    rel = 'data/1AAY_alphafold/test_676a7_unrelaxed_rank_1_model_2.pdb'
    p = here.parent / rel
    model = Fold(p)
    
    ca = list(get_alpha_carbon_atoms(model, only_coords=True))
    dist = DistanceBand(ca, 8, p=2, binary=True)

    x = {}
    for i in range(len(model)):
        x[i] = list(is_close(i, model, 8, 'alpha_carbons'))

    for i in range(len(model)):
        foo = [ix for ix, u in enumerate(x[i]) if u and ix != i] 
        bar = list(dist.neighbors[i])
        assert foo == bar
