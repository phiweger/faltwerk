from collections import defaultdict
# import shutil
import subprocess
import tempfile
try:
    from typing import Union, List
except ImportError:
    from typing_extensions import Union, List

from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from libpysal.weights import DistanceBand
import numpy as np
import screed

from faltwerk.models import Fold


def get_alpha_carbon_atoms(x: Union[Atom, Residue], only_coords=False):
    '''
    alpha carbon: https://foldit.fandom.com/wiki/Alpha_carbon
    '''
    if type(x) == Fold:
        residues = x.structure.get_residues()
    else:
        residues = x.get_residues()
    for res in residues:
        l = []
        for i in res.get_atoms():
            # There is only one alpha carbon per amino acid.
            if i.get_id() == 'CA':
                if only_coords:
                    yield i.coord
                else:
                    yield i


def get_coordinate(x: Union[Atom, Residue]):
    '''
    Calculating center of mass is much slower than looking up the coordinate
    of a carbon atom.
    '''
    if type(x) == Residue:
        return x.center_of_mass()
    elif type(x) == Atom:
        return x.coord
    else:
        raise ValueError('Unsupported type')


def euclidean_distance(a, b):
    '''
    https://stackoverflow.com/questions/1401712/how-can-the-euclidean-distance-be-calculated-with-numpy
    '''
    return np.linalg.norm(a-b)


def is_close(pos, fold, radius, coordinates='alpha_carbons'):
    '''
    fold = Fold('test_676a7_unrelaxed_rank_1_model_2.pdb')
    list(is_close(1, fold, 10))
    # [True, True, True, True, True, False, False, ...]

    coordinates .. center_of_mass, alpha_carbon

    TODO: pseudo single-atom representation of side chains:

    > Specifically, we defined this distance according to the sites' side chain
    center of masses. A consequence of approximating DTL with respect to the
    closest ligand-binding sites is that by definition, any ligand-binding
    residue has a DTL of 0. -- Kiefl et al.,
    https://www.biorxiv.org/content/10.1101/2022.03.02.482602v1.full.pdf

    - https://pymolwiki.org/index.php/Sidechaincenters
    - https://bioinformatics.stackexchange.com/questions/18162/pymol-python-script-for-selecting-a-residues-sidechain-and-calculating-its-cent
    '''

    if coordinates == 'center_of_mass':
        chain = list(fold.structure.get_residues())
    
    elif coordinates == 'alpha_carbons':
        chain = list(get_alpha_carbon_atoms(fold))
    
    a = get_coordinate(chain[pos])

    for i in chain:
        b = get_coordinate(i)
        # https://stackoverflow.com/questions/1401712/how-can-the-euclidean-distance-be-calculated-with-numpy
        # dist = np.linalg.norm(a - b)
        dist = euclidean_distance(a, b)
        if dist < radius:
            yield True
        else:
            yield False


def get_foldseek_vae_states(fold):
    '''
    https://github.com/steineggerlab/foldseek/issues/15
    '''
    tmp = tempfile.TemporaryDirectory()
    p = tmp.name
    fp = fold.path.resolve().__str__()

    steps = [
        f'foldseek createdb {fp} {p}/db',
        f'foldseek lndb {p}/db_h {p}/db_ss_h',
        f'foldseek convert2fasta {p}/db_ss {p}/db_ss.fasta',
    ]
    command = '; '.join(steps)
    log = subprocess.run(command, capture_output=True, shell=True)
    assert log.returncode == 0, log.stderr

    # shutil.copyfile(f'{p}/db_ss.fasta', outfile)
    with screed.open(f'{p}/db_ss.fasta') as file:
       return next(file).sequence


def get_foldseek_vae_states_from_path(fp):
    '''
    Turn the model's states and aa sequence into tokens
    '''
    model = Fold(fp)
    return get_foldseek_vae_states(model)


def distance_to_closest_active_site(fold, binding_frequencies, threshold=0.5):
    '''
    Usage:

    f = Fold('serine_hydroxymethyltransferase.pdb')
    b = Binding(f, 'confident')
    b.predict_binding_(pfam)
    bf = b.get_binding('PF00464.18', 'SER')
    distance_to_closest_active_site(f, bf, .5)
    '''
    residues = list(fold.structure.get_residues())
    bf = binding_frequencies

    assert len(residues) == len(bf) 
    active = [r for r, f in zip(residues, bf) if f >= threshold]

    l = []
    for r in residues:
        rm = r.center_of_mass()
        d = np.min([euclidean_distance(rm, a.center_of_mass()) for a in active])
        l.append(float(d))
    
    return l


def get_complex_interface(chains: List, angstrom=10, map_to_chains=False):
    '''
    from foldspace.models import Complex
    from foldspace.geometry import get_complex_interface

    cx = Complex('/path/to/model.pdb')
    interface = get_complex_interface(cx, 10)

    map_to_chains .. remember which interface atom belongs to which chain(s)
    '''
    coords = []
    labels = []
    ix_residues = []
    pairwise = []

    for chain in chains:
        # TODO: Add a more general iterator, for both complexes and single folds
        for atom in chain.get_atoms():
            if atom.get_id() == 'CA':
                coords.append(get_coordinate(atom))
                labels.append(chain.id)
                # For each chain, residue numbering starts at 1
                ix_residues.append(atom.parent.id[1])
                
    lu = {i: j for i, j in enumerate(labels)}
    # lu .. lookup: {712: 'C', 713: 'C', 714: ...
    lu_res = {i: j for i, j in enumerate(zip(labels, ix_residues))}  
    # {0: ('B', 1), 1: ('B', 2), 2: ...

    dist = DistanceBand(coords, angstrom, p=2, binary=True)
    # {0: [1, 2], 1: [0, 2, 3, 4], 2: ...
    # assert len(dist.neighbors.keys()) == len(cx)
    
    interface = set()
    for k, v in dist.neighbors.items():
        origin = lu[k]
        for i in v:
            if lu[i] != origin:
                interface.add(k)
                pairwise.append([lu_res[k], lu_res[i]])
    
    #import pdb
    #pdb.set_trace()

    if map_to_chains:
        return [(k, lu_res[k]) for k in interface], pairwise
    else:
        return interface


def get_interface(A: Chain, B: Chain, angstrom=10):
    '''
    Get all residue positions that are likely binding partners btw/ 2 chains.

    Example:

    get_interface(cx.chains['B'], cx.chains['D'])
    '''
    d = defaultdict(set)
    _, pairwise = get_complex_interface(
        [A, B], map_to_chains=True, angstrom=angstrom)
    for (k1, v1), (k2, v2) in pairwise:
        d[k1].add(v1)
    return dict(d)  # return dict not defaultdict to avoid unexpect behav


def distance_to_positions(model, positions):
    '''
    import altair as alt

    from foldspace.models import Complex
    from foldspace.geometry import get_complex_interface, distance_to_positions

    cx = Complex('/path/to/model.pdb')
    interface = get_complex_interface(cx, 10)
    distance_to_interface = distance_to_positions(model, interface)
    df['distance_to_interface'] = distance_to_interface
    
    alt.Chart(df).mark_boxplot(extent=1.5).encode(
        x='selection:O',
        y='distance_to_interface:Q',
    )
    '''
    atoms = list(get_alpha_carbon_atoms(model, only_coords=True))
    l = []
    for atom in atoms:
        # Restrict positions to only those in the protein structure, ignore the rest
        dist = np.min([euclidean_distance(atom, atoms[p]) for p in positions if p in range(len(model))])
        l.append(dist)
    return l

