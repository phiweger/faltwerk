from itertools import product

import numpy as np

from faltwerk.geometry import get_foldseek_vae_states

# Assign each residue-state pair a separate token out of 400; reserve 0 for
# mask; eg AA > 1, AB > 2, ... // foldseek states use the amino acid alphabet.
alphabet = 'ACDEFGHIKLMNPQRSTVWY'
tokens = {i+j: n+1 for n, (i, j) in enumerate(product(alphabet, alphabet))}
assert len(tokens) == len(alphabet)**2


def tokenize(model):
    '''
    TODO: Don't initialize lookup table every time

    Turn the model's states and aa sequence into tokens
    '''
    z = zip(model.sequence, get_foldseek_vae_states(model))
    return np.array([tokens[i+j] for i, j in z])



