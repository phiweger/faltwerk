from pathlib import Path

import numpy as np

import faltwerk
from faltwerk.learn import tokenize
from faltwerk.models import Fold


def test_tokenize():
    fp = 'data/alphafold2/transferrin/test_08df6_unrelaxed_rank_1_model_3.pdb'
    fp = Path(faltwerk.__file__).parent.parent / fp
    model = Fold(fp)
    assert all(tokenize(model)[:5] == np.array([203, 293, 193,  13, 353]))
    return None