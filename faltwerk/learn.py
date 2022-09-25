'''
Below I experiment, not intended for use!
'''


from itertools import product
from multiprocess import Pool, cpu_count
# https://stackoverflow.com/questions/41385708/multiprocessing-example-giving-attributeerror
# https://stackoverflow.com/questions/5666576/show-the-progress-of-a-python-multiprocessing-pool-imap-unordered-call

import numpy as np
from tqdm import tqdm

from faltwerk.geometry import (
    get_foldseek_vae_states,
    get_foldseek_vae_states_from_path,
    )
from faltwerk.models import Fold


# Assign each residue-state pair a separate token out of 400; reserve 0 for
# mask; eg AA > 1, AB > 2, ... // foldseek states use the amino acid alphabet.
alphabet = 'ACDEFGHIKLMNPQRSTVWY'
tokens = {i+j: n+3 for n, (i, j) in enumerate(product(alphabet, alphabet))}
assert len(tokens) == len(alphabet)**2


def tokenize_states_and_seq(fp):
    '''
    Turn the model's states and aa sequence into tokens
    '''
    model = Fold(fp)
    z = zip(model.sequence, get_foldseek_vae_states(model))
    x = [bos_token_id] + [tokens[i+j] for i, j in z] + [eos_token_id]
    return np.array(x)


def process_many(files, fn=tokenize_states_and_seq, threads=cpu_count()):
    with Pool(threads) as pool:
        results = []
        for i in tqdm(pool.imap_unordered(fn, files), total=len(files)):
            results.append(i)
    return results



class FoldDataset():
    pass


class FoldTokenizer():
    '''
    https://huggingface.co/docs/transformers/model_doc/roberta#transformers.RobertaConfig
    '''
    def __init__(self):
        self.bos_token_id = 0  # .. beginning of sentence 
        self.pad_token_id = 1
        self.eos_token_id = 2  # .. end of sequence

        alphabet = 'ACDEFGHIKLMNPQRSTVWY'
        pr = enumerate(product(alphabet, alphabet))
        self.tokens = {i+j: n+3 for n, (i, j) in pr}
        # assert len(tokens) == len(alphabet)**2

    def __call__(self, fp):
        model = Fold(fp)
        z = zip(model.sequence, get_foldseek_vae_states(model))
        x = [self.bos_token_id] + [self.tokens[i+j] for i, j in z] + [self.eos_token_id]
        return np.array(x)

