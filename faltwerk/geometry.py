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
    Given a (Biopython) Atom or Residue object, return a generator of backbone
    carbon atoms.
    
    - https://foldit.fandom.com/wiki/Alpha_carbon
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
    Calculate the center of mass for an atom or residue. Note: Calculating the
    center of mass is much slower than looking up the coordinate of a carbon
    atom.
    '''
    if type(x) == Residue:
        return x.center_of_mass()
    elif type(x) == Atom:
        return x.coord
    else:
        raise ValueError('Unsupported type')


def euclidean_distance(a, b):
    '''
    Calculate the Euclidean disance between two points.

    https://stackoverflow.com/questions/1401712/how-can-the-euclidean-distance-be-calculated-with-numpy
    '''
    return np.linalg.norm(a-b)


def is_close(pos, fold, radius, coordinates='alpha_carbons'):
    '''
    Find all neighbors within a specified distance of a residue.

    Usage:

    >>> m = Fold('test.pdb')
    >>> list(is_close(1, fold, 10))
    >>> # [True, True, True, True, True, False, False, ...]

    Optional arguments:

    - radius -- Angstrom radius around residue of interest
    - coordinates -- "center_of_mass" or "alpha_carbon" (default)
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
    ``foldseek`` is a program to search similar protein structures. Methodically
    it is really clever as it maps each residue to a latent state that
    represents local geometry, and then searches homologous proteins based on
    this VAE alphabet.

    This function is for exploration only, maybe those states can be useful when
    training neural nets.

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
    Same as ``get_foldseek_vae_states`` but from path.
    '''
    model = Fold(fp)
    return get_foldseek_vae_states(model)


def distance_to_closest_active_site(fold, binding_frequencies, threshold=0.5):
    '''
    Calculate the (Euclidean) distance of each residue to the closest active
    site (or really any annotation, for that matter).

    Usage:

    >>> m = Fold('serine_hydroxymethyltransferase.pdb')
    >>> b = Binding(m, 'confident')
    >>> b.predict_binding_(hmms)
    >>> bf = b.get_binding('PF00464.18', 'SER')
    >>> distance_to_closest_active_site(m, bf, .5)

    Optional arguments:

    - threshold -- cut-off what counts as binding site (default 0.5)
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
    Given a Complex object, ie a protein structure complex, return all
    binding interfaces.

    Usage:

    >>> from faltwerk import Complex
    >>> from foldspace.geometry import get_complex_interface
    >>> cx = Complex('/path/to/model.pdb')
    >>> interface = get_complex_interface(cx, 10)

    Optional arguments:

    - map_to_chains -- remember which interface atom belongs to which chain(s) (default True)
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

    Usage:

    >>> get_interface(cx.chains['B'], cx.chains['D'])

    Optional arguments:

    - angstrom -- maximum distance of what counts as interface (default 10)
    '''
    d = defaultdict(set)
    _, pairwise = get_complex_interface(
        [A, B], map_to_chains=True, angstrom=angstrom)
    for (k1, v1), (k2, v2) in pairwise:
        d[k1].add(v1)
    return dict(d)  # return dict not defaultdict to avoid unexpect behav


def distance_to_positions(model, positions):
    '''
    Sometimes it is useful to analyse the distance of residues to some feature
    such as a binding interface.

    See eg:

    https://www.biorxiv.org/content/10.1101/2022.03.02.482602v1.full.pdf

    Usage:

    >>> import altair as alt
    >>> from faltwerk import Complex
    >>> from faltwerk.geometry import get_complex_interface, distance_to_positions
    >>> cx = Complex('/path/to/model.pdb')
    >>> interface = get_complex_interface(cx, 10)
    >>> distance_to_interface = distance_to_positions(model, interface)

    >>> # Create a dataframe with residue numbers and whether the sites are
    >>> # under positive selection or not.
    >>> df['distance_to_interface'] = distance_to_interface

    >>> alt.Chart(df).mark_boxplot(extent=1.5).encode(
    >>>    x='selection:O',
    >>>    y='distance_to_interface:Q',
    >>> )
    '''
    atoms = list(get_alpha_carbon_atoms(model, only_coords=True))
    l = []
    for atom in atoms:
        # Restrict positions to only those in the protein structure, ignore the rest
        dist = np.min([euclidean_distance(atom, atoms[p]) for p in positions if p in range(len(model))])
        l.append(dist)
    return l

