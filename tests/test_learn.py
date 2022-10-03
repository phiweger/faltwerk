from pathlib import Path

import numpy as np

import faltwerk
from faltwerk.geometry import get_foldseek_vae_states
from faltwerk.learn import tokenize_states_and_seq, tokens
from faltwerk.models import Fold


def test_tokenize():
    fp = 'data/alphafold2/transferrin/test_08df6_unrelaxed_rank_1_model_3.pdb'
    fp = Path(faltwerk.__file__).parent.parent / fp
    assert all(tokenize_states_and_seq(fp)[:5] == np.array([203, 293, 193,  13, 353]))

    model = Fold(fp)
    raw = model.sequence[0] + get_foldseek_vae_states(model)[0]
    assert tokens[raw] == 203

    return None