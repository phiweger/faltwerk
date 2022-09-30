from collections import defaultdict, Counter
from copy import deepcopy
import io
import json
from pathlib import Path
import pkg_resources
import re
try:
    from typing import Union
except ImportError:
    from typing_extensions import Union
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd

from faltwerk.io import read_pdb, save_pdb, load_scores, stream
from faltwerk.parsers import HMMERStandardOutput
from faltwerk.utils import (
    align_structures,
    filter_aa,
    get_sequence,
    search_domains,
    )


# TODO: Complex(fp, one_fold_per='chain')
# or
# folds = cx.split_into_folds()
# subset chain, to stream, load into fold (from_stream=True)
class Complex():
    '''
    Datastructure to deal with a multiple chain
    
    TODO: A Complex object should be made up of two or more Fold objects
    '''
    def __init__(self, fp, scores=None, reindex: bool=False, reindex_start_position: int=1):
        self.path = Path(fp)

        # TODO: This will only read in the first sequence! Thus
        # cx.sequence does not hold all sequences
        self.structure = read_pdb(self.path, strict=False, reindex=reindex, reindex_start_position=reindex_start_position)
        
        self.chains = {}
        self.annotation = {}

        original_ix = []
        for chain in self.structure.get_chains():
            self.chains[chain.id] = chain
            original_ix.extend([i for i in range(len(chain))])

        # Positions for the overall complex are sorted by name. So if there
        # are chains (B, C, D) then 0 marks the first position in B.
        absolute_pos, relative_pos, names = [], [], []
        cnt = 0
        for chain in [v for (k, v) in sorted(self.chains.items())]:
            # print(chain)
            # aa = filter_aa(chain)
            
            # TODO: i+1 bc/ 3Dmol.js and PDB count positions from 1?
            # We go with 0-based here bc/ pythonic and intuitive.
            p, q = zip(*[(i, chain.id) for i in range(len(chain))])
            
            for i in range(len(chain)):
                absolute_pos.append(cnt)
                cnt += 1

            relative_pos.extend(p)
            names.extend(q)

        assert original_ix == relative_pos, 'Unsorted PDB file'
        self.len = cnt
        # Provide different indices for chains, which can be useful in plotting
        self.annotate_('positions', absolute_pos)
        self.annotate_('chain_positions', relative_pos)
        self.annotate_('chains', names)

        s = sorted(set(self.annotation['chains']))
        d = {i: j for i, j in zip(s, range(len(s)))}
        self.annotate_('chains_numbered', [d[k] for k in self.annotation['chains']])
        self.chains_numbered = d

        # In complexes, there is the index position, and the relative one of the
        # chain.
        self.scores = {}
        if scores:
            sc = load_scores(scores)
            for k, v in sc.items():
                try:
                    self.annotate_(k, v)
                except TypeError:
                    # Skip single values like "max_pae"
                    continue
        return None

    def get_chains(self):
        v = self.chains.values()
        for chain in v:
            sorted_chains = [x for _, x in sorted(
                zip([len(i) for i in v], v),
                key=lambda x: x[0],
                reverse=True
                )]
            return sorted_chains

    def annotate_(self, key, values, check=True):
        if check:
            assert len(values) == len(self), 'No 1:1 mapping of labels to positions'
        self.annotation[key] = values
        return None

    def __len__(self):
        return self.len

    def to_stream(self):
        return stream(self.structure)

    def __getitem__(self, ix):
        '''
        cx[1]
        # <Chain id=C>
        '''
        k = {v: k for k, v in self.chains_numbered.items()}[ix]
        return self.chains[k]


class Fold():
    '''
    Datastructure to deal with a single chain
    '''
    def __init__(self, fp, quiet=True, annotate=True, strict=True):
        self.path = Path(fp)

        if not quiet:
            print(f'Loading structure in {self.path.name}')
        self.structure = read_pdb(self.path, strict=strict)
        # Reads in ONE sequence; of more chains are present, treat as complex
        self.sequence = get_sequence(next(self.structure.get_chains()))
        # structure > model > chain > residue > atom
        self.transformed = False
        self.annotation = {}

        if annotate:
            aa = filter_aa(self.structure)
            self.annotate_('position', [i+1 for i in range(len(aa))])
        return None

    def __repr__(self):
        return 'Fold loaded from ' + self.path.name

    def __iter__(self):
        yield from self.structure

    def __len__(self):
        '''Number of amino acids in the sequence'''
        # return len(list(self.structure.get_residues()))
        return len(self.sequence)

    def align_to(self, ref, mode=0, minscore=0.5):
        tmscore, rot, tra = align_structures(ref.path, self.path, mode=mode, minscore=minscore)
        cp = deepcopy(self)
        cp.structure.transform(rot, tra)
        cp.transformed = True
        return tmscore, cp

    def rename_chains_(self, renames: dict) -> None:
        '''
        - https://stackoverflow.com/questions/70246451/how-do-i-change-the-chain-name-of-a-pdb-file
        - modifies in place, "_" suffix convention like in pytorch
        '''
        for model in self.structure:
            for chain in model:
                old_name = chain.get_id()
                new_name = renames.get(old_name)
                if new_name:
                    print(f'Renaming chain {old_name} to {new_name}')
                    chain.id = new_name
                else:
                    print(f'Keeping chain name {old_name}, no new name found')
        return None

    def add_scores(self, fp):
        with open(fp, 'r') as file:
            scores = json.load(file)
            self.annotate_('plddt', scores['plddt'])
        return None

    def to_stream(self):
        stream = io.StringIO()
        return save_pdb(self.structure, stream).getvalue()

    def annotate_(self, key, values, check=True):
        if check:
            assert len(values) == len(self), f'No 1:1 mapping of labels to positions for key "{key}"'
        self.annotation[key] = values
        return None

    def annotate_many_(self, many):
        for k, v in many.items():
            self.annotate_(k, v)
        return None

    def select(self):
        '''
        Assume one model and one chain?

        structure > model > chain > residue > atom
        
        currently, there are no eg atom objects implemented in py3Dmol:
        https://github.com/3dmol/3Dmol.js/issues/498

        select: 
        0:A::CA
        ::343-368,85-89,734:
        l = ','.join([1, 2, 3, 4, 5, 66, 67, 88, 89, 90])
        
        f'::{l}:'
        if unique, 343-368 should work too
        
        check type, if string decompose, if list assume/ check everything 
        else unique
        '''

        # Maybe this should just be a PDB to json fn from utils
        # -- also, use existing code? pdb-tools:
        # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6343223/        
        pass

    # def __iter__(self, item):
    #     self.fold.select('residues')[item]


class AlphaFold():
    '''
    fold = AlphaFold(...)
    fold.best -> structure, wrapped in Fold object (keeps track of filepath ...)
    '''
    def __init__(self, indir, workdir=None):
        self.models = {}
        self.workdir = workdir

        for n, i in enumerate(self.read_alphafold(indir, workdir)):
            self.models[n] = i

        self.best = self.models[0]
        
        return None

    def read_alphafold(self, filedir, outdir=None):
        '''
        Load structures and (quality) scores, align them, and put prepare data
        for visualization.
        '''
        files = Path(filedir).glob('*.pdb')
        d = {}
    
        # Load scores
        for i in files:
            model = int(re.match(r'.*model_(\d).pdb', i.name).group(1))
            # print(f'Loading model {model}')
            fp = str(i.resolve())
            fp = fp.replace(f'{model}.pdb', f'{model}_scores.json')
            
            with open(fp, 'r') as file:
                scores = json.load(file)
           
            fold = Fold(i)
            fold.annotate_('plddt', scores['plddt'])
    
            v = np.mean(scores['plddt'])
            d[fold] = v
    
        # Rank models by pLDDT, best is reference
        try:
            ref, *queries = [i for i, j in sorted(d.items(), key=lambda x: x[1], reverse=True)]
        except ValueError:
            raise FileNotFoundError()

        print(f'Best model (pLDDT): {ref.path.name}')
        print(f'Align remaining models to best and rename')
        # Align into the same space
        rest = []
        for qry, chain in zip(queries, 'BCDE'):
            _, trx = qry.align_to(ref)
            trx.rename_chains_({'A': chain})
            rest.append(trx)
    
        if outdir:
            outdir = Path(outdir)
            
            if not outdir.exists():
                outdir.mkdir(parents=True)
    
            name = ref.path.name.replace('.pdb', '.reference.pdb')
            _ = save_pdb(ref.structure, outdir / name)
    
            for qry in rest:
                name = qry.path.name.replace('.pdb', '.transform.pdb')
                _ = save_pdb(qry.structure, outdir / name)
    
        return [ref] + rest


class Binding():
    '''
    f = Fold(...)
    b = Binding(f, fp_interactions)
    b.predict_binding_(fp_hmms)

    - get sequeence from pdb file
    - hmmsearch
    - parse result
    - align seqs
    - expose (domain, ligand pairs)

    d
    # prints all interactions, same as
    d.interactions
    # Could be a dict
    # ('PFxxxx', 'ZN'): [...] .. array of binding frequencies
    # d.fold.sequence .. get original sequence
    # type(d.fold) is Fold

    interaction = d.get_interaction(['PFxxxx', 'ZN'])
    interaction = d[('PFxxxx', 'ZN')]

    view = plot_interaction(interaction)
    '''
    def __init__(self, fold, option='confident'):
        self.fold = fold
        self.options = {
            'confident': 'InteracDome_v0.3-confident.tsv',
            'representable': 'InteracDome_v0.3-representable.tsv',
            'non-redundant': 'InteracDome_v0.3-representableNR.tsv',
            }
        self.interactions = self.read_interactions(self.options[option])
        self.domains = None
        self.ids = None

    def read_interactions(self, fn):
        fp = pkg_resources.resource_filename('faltwerk', f'/data/ligands/{fn}')
        # fp = Path(__file__).parents[1] / f'data/ligands/{fn}'
        df = pd.read_csv(fp, sep='\t', comment='#')
        return df

    def predict_binding_(self, hmms):
        self.domains = search_domains(self.fold, hmms)
        self.ids = set(self.domains['acc'])
        # TODO: filter?

        # Note: By using a set we assume there are no two domain: ligand pairs
        # for the same ligand, eg zinc binding to the start AND end of a domain.
        ligands = defaultdict(set)
        self._pfam_map = {}
        '''
        PF00096_zf-C2H2 NUCACID_
        PF00096_zf-C2H2 DNABASE_
        PF13894_zf-C2H2_4 DNABACKBONE_
        PF13894_zf-C2H2_4 ZN
        '''
        for _, i in self.interactions.iterrows():
            for j in set(self.domains['acc']):
                # InteractDome does not preserve Pfam versions, they used v31
                if j.split('.')[0] in i["pfam_id"]:
                    # print(i.pfam_id, i.ligand_type)
                    ligands[j].add(i.ligand_type)
                    self._pfam_map[j] = i.pfam_id
        
        # "Turn off" defaultdict
        self.ligands = dict(ligands)

    def get_binding(self, domain, ligand):
        '''
        Returns binding frequencies for protein sequence of fold
        '''
        # Get binding frequencies for (domain, ligand) pair
        if not isinstance(self.domains, pd.DataFrame):
            raise ValueError('Please annotate domains first using "predict_binding_()"')

        n = self.interactions
        bf = list(map(float, n[(n['pfam_id'] == self._pfam_map[domain]) & (n['ligand_type'] == ligand)]['binding_frequencies'].item().split(',')))
        # Match query sequence residues to domain states
        df = self.domains
        # We can have multiple domain hits in one protein, even of the same
        # domain, think zinc finger (PDB 1AAY).
        result = np.zeros(len(self.fold))
        for _, i in df[df['acc'] == domain].iterrows():
            # We reverse the vector so when we pop()
            u = bf[::-1]
            v = []
            for residue, state in zip(i.sequence_align, i.match_state_align):
                if residue == '-':
                    # Domain state not present in residue (inserted "-"s)
                    _ = u.pop()
                    continue

                if state != '.':
                    v.append(u.pop())
                else:
                    v.append(0.)

            n_gaps = Counter(i.sequence_align)['-']
            assert len(i.sequence_align) - n_gaps == len(v)
            # d[(i.acc, j.ligand_type)] = u
            # print(i.ali_start, i.ali_stop)
            # print(len(v))
            # print(len(result))
            # print(len(result[i.ali_start-1:i.ali_stop]))
            result[i.ali_start-1:i.ali_stop] = v

        return result

    def get_domain(self, domain):
        '''
        Returns vector of residue positions with 1 where domain present and
        0 otherwise.

        for i in b.ids:
            _ = b.get_domain(i)
        '''
        
        if not isinstance(self.domains, pd.DataFrame):
            raise ValueError('Please annotate domains first using "predict_binding_()"')

        # q = 'PF00405.16'
        ranges = []
        for i in b.domains[b.domains['acc'] == q].itertuples():
            ranges.append([i.ali_start, i.ali_stop])
    
        result = np.zeros(len(self.fold))
        for r in ranges:
            for i in range(*r):
                result[i] = 1

        return result


        


        

'''


with screed.open("rcsb_pdb_1AAY.fasta") as file:
    for line in file:
        seq = line.sequence
        ln = len(seq)


x = HMMERStandardOutput("result.txt")

result = np.zeros(ln)

# d = {}
for _, i in x.dom_hits.iterrows():
    sub = df[[i.acc.split('.')[0] in name for name in df['pfam_id']]]

    for _, j in sub.iterrows():
        if not j.ligand_type == 'ZN' or not i.acc == 'PF00096.25':
            continue

        v = list(map(float, j.binding_frequencies.split(',')))[::-1]
        # We reverse the vector so when we pop() we get the first, then 2nd ...
        print(seq[i.ali_start-1:i.ali_stop])
    
        print(i.sequence_align)
        print(i.match_state_align)


        u = []
        for residue, state in zip(i.sequence_align, i.match_state_align):
            if state != '.':
                u.append(v.pop())
            else:
                u.append(0.)
        assert len(i.sequence_align) == len(u)
        # d[(i.acc, j.ligand_type)] = u
        result[i.ali_start-1:i.ali_stop] = u


with open('zinc.csv', 'w+') as out:
    out.write(','.join(map(str, result)) + '\n')

'''






