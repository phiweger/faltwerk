from collections import defaultdict, Counter
from copy import deepcopy
import io
import json
from pathlib import Path
import pkg_resources
import re
import subprocess
import tempfile
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


class Complex():
    '''
    Core object to manipulate multi-component protein structures.

    Many proteins interact in protein-protein complexes, and there is an
    increasing number of models such as AlphaFold v2 that can be used to
    predict such binding patterns. The ``Complex`` object collects methods
    to manipulate multi-component structures, selecting individual components,
    annotating binding sites, etc.

    Basic usage:

    >>> from faltwerk import Complex, Layout
    >>> cx = Complex('af2_multichain_prediction.pdb', reindex=1)
    
    Remove and select chains:
    
    >>> cx - 'AB'  # rm chains, same as cx - ['A', 'B']
    >>> cx * 'C'   # select
    >>> cx.chains['C']
    >>> cx[2]

    Visualise:

    >>> Layout(cx).geom_ribbon().render()
    '''
    def __init__(self, fp, scores=None, reindex: bool=False, reindex_start_position: int=1):
        '''
        Create a ``faltwerk.Complex`` object.

        Optional arguments:

        - scores (default None) -- add quality scores in ColabFold format [json]
        - reindex (default False) -- reindex residues (if not numbered 1:n)
        - reindex_start_position (default 1) -- where to start residue index
        '''
        self.path = Path(fp)
        
        # "The Structure Object" -- https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ
        self.structure = read_pdb(
            self.path,
            strict=False,
            reindex=reindex,
            reindex_start_position=reindex_start_position)
        
        self.chains = {}
        self.annotation = {}

        original_ix = []
        for chain in self.structure.get_chains():
            self.chains[chain.id] = chain
            original_ix.extend([i for i in range(len(chain))])

        self.sequences = {i: get_sequence(ch) for i, ch in self.chains.items()}

        # Positions for the overall complex are sorted by name. So if there
        # are chains (B, C, D) then 0 marks the first position in B.
        absolute_pos, relative_pos, names = [], [], []
        cnt = 0
        for chain in [v for (k, v) in sorted(self.chains.items())]:
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
    
    def delete_chains(self, chains):
        return delete_chains(self, chains)
    
    def select_chains(self, chains):
        return delete_chains(self, chains, invert=True)
    
    def __sub__(self, chains):
        return delete_chains(self, chains)
    
    def __mul__(self, chains):
        return delete_chains(self, chains, invert=True)


class Fold():
    '''
    Core object to manipulate single-component protein structures.

    The ``Fold`` object is the basis for many protein-centric operations in the
    ``faltwerk`` library.

    Basic usage:

    >>> from faltwerk import Fold
    >>> model = Fold('af2_prediction.pdb')
    '''
    def __init__(self, fp, quiet=True, annotate=True, strict=True, reindex: bool=False, reindex_start_position: int=1):
        '''
        Create a ``faltwerk.Fold`` object.

        Optional arguments:

        - strict (default True) -- check that single chain
        - annotate (default True) -- add a track that annotates positions N > C 
        '''
        self.path = Path(fp)

        if not quiet:
            print(f'Loading structure in {self.path.name}')
        self.structure = read_pdb(self.path, strict=strict, reindex=reindex, reindex_start_position=reindex_start_position)
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
        '''
        Align a structure to another one using ``foldseek``. Returns the
        Tm score (> 0.5 is good) and a copy of the query structure:

        >>> qry = Fold(...)
        >>> ref = Fold(...)
        >>> tm_score, cp_qry = qry.align_to(ref)

        Optional arguments (passed to ``foldseek``, see
        https://github.com/steineggerlab/foldseek):

        - mode (default 0)
        - minscore (default 0.5)
        '''
        tmscore, rot, tra = align_structures(ref.path, self.path, mode=mode, minscore=minscore)
        cp = deepcopy(self)
        cp.structure.transform(rot, tra)
        cp.transformed = True
        return tmscore, cp

    def rename_chains_(self, renames: dict) -> None:
        '''
        Rename chains inplace (note the pytorch-style "_" fn name suffix). AF2
        will for example run multiple models, to load them, each has to have a
        unique chain name.

        >>> model = Fold(...)
        >>> model.rename_chains_({'A': 'foo', 'B': 'bar'})  # inplace operation
        '''
        # https://stackoverflow.com/questions/70246451/how-do-i-change-the-chain-name-of-a-pdb-file
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
        '''
        Add ColabFold formatted quality scores.
        '''
        with open(fp, 'r') as file:
            scores = json.load(file)
            self.annotate_('plddt', scores['plddt'])
        return None

    def to_stream(self):
        stream = io.StringIO()
        return save_pdb(self.structure, stream).getvalue()

    def annotate_(self, key, values, check=True):
        '''
        Annotate a structure with a track of some feature (solvent access,
        selected sites, surface probability, ...), with one value for each
        residue. Inplace operation.

        Optional arguments:

        - check (default True) -- length features == length protein?

        >>> from faltwerk.io import load_bfactor_column
        >>> features = load_bfactor_column('dMASIF.pdb')
        >>> model = Fold(...)
        >>> model.annotate_('surface', features)
        '''
        if check:
            assert len(values) == len(self), f'No 1:1 mapping of labels to positions for key "{key}"'
        self.annotation[key] = values
        return None

    def annotate_many_(self, many: dict):
        '''
        Batch annotate features inplace.

        >>> model.annotate_many_({'foo': arr1, 'bar': arr2})
        '''
        for k, v in many.items():
            self.annotate_(k, v)
        return None


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
    Object to predict ligand binding based on mapped Pfam domains, using the
    method first explored in the "InteracDome":

    - https://protdomain.princeton.edu/interacdome
    - https://merenlab.org/2020/07/22/interacdome/
    - https://academic.oup.com/nar/article/47/2/582/5232439

    The basic idea there was to take all PDB structures with interacting
    ligands, mark the interface between the two on the (linear) sequence, and
    then map the (linear) sequence to the protein sequence underlying a
    protein structure, marking the likely interacting residues on this query
    structure. There can be multiple ligands for each Pfam domain.

    Because this approach was "trained" on Pfam v31, and all coordinates are
    based on this version, please only use v31 as reference.

    There are three sets of contact data to make predictions from:

    ``confident`` .. "correspond to domain-ligand interactions that had 
    nonredundant instances across three or more distinct PDB entries and 
    achieved a cross-validated precision of at least 0.5. [Kobren and Singh] 
    recommend using this collection to annotate potential ligand-binding 
    positions in protein sequences"
    
    ``representable``/ ``nonredundant`` .. "correspond to domain-ligand 
    interactions that had nonredundant instances across three or more distinct 
    PDB structures. [Kobren and Singh] recommend using this collection to learn 
    more about domain binding properties" (``nonredundant`` is the nonredundant
    version of ``representable``)

    `Update`: For an end-to-end approach to the same problem, see DiffDock 
    -- https://github.com/gcorso/DiffDock

    Usage:

    >>> from faltwerk import Fold, Binding, Layout
    >>> model = Fold('faltwerk/data/zinc_finger/1mey.pdb)
    >>> hmms = 'path/to/pfam_v31/Pfam-A.hmm'
    >>> b = Binding(model, option='non-redundant')
    >>> b.predict_binding_(hmms)  # inplace operation
    >>> # Explore: b.domains, b.ligands
    >>> model.annotate_('binding', b.get_binding('PF13912.5', 'ZN'))
    >>> Layout(model).geom_surface('binding').render()
    >>> # Marked in yellow are two symmetric iron binding sites
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

        >>> for i in b.ids:
        >>>    _ = b.get_domain(i)
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


def delete_chains(model, chains, invert=False):
    '''
    # http://www.bonvinlab.org/pdb-tools/
    python pdb_delchain.py -A,B 1CTF.pdb

    When invert, instead of deleting the chains, keep them and remove all others.
    '''    
    with tempfile.TemporaryDirectory() as p:
    # p = tmp.name
        fp = f'{p}/fold.pdb'
        save_pdb(model.structure, fp)

        if invert:
            chains = [i for i in model.chains.keys() if i not in chains]
        
        remaining = len(list(model.structure.get_chains())) - len(chains)

        steps = [
            f'pdb_delchain -{",".join(chains)} {fp} > {p}/deleted.pdb',
        ]
        command = '; '.join(steps)
        log = subprocess.run(command, capture_output=True, shell=True)
        assert log.returncode == 0, log.stderr
        
        if remaining == 1:
            new_model = Fold(f'{p}/deleted.pdb', reindex=True)    
        else:
            new_model = Complex(f'{p}/deleted.pdb', reindex=True)
        return new_model




