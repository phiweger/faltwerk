import collections
from itertools import combinations
from math import log, e
from pathlib import Path
import subprocess
try:
    from typing import Union
except ImportError:
    from typing_extensions import Union
import tempfile

from Bio import SeqUtils
from Bio.PDB import PDBIO
from Bio.PDB.Structure import Structure
import numpy as np

from faltwerk.parsers import HMMERStandardOutput


def chunks(l, n):
    '''
    Yield successive n-sized chunks from l (stackoverflow, 312443).

    a = [1, 2, 3, 4]
    list(chunks(a, 2))
    # [[1, 2], [3, 4]]

    Returns empty list if list empty.

    For overlapping chunks, see windows()
    '''
    for i in range(0, len(l), n):
        yield l[i:i + n]


def rmdir(directory):
    '''
    DEPRECATED, using tempfile ".cleanup()" fn
    https://stackoverflow.com/questions/13118029/deleting-folders-in-python-recursively/49782093#49782093
    '''
    directory = Path(directory)
    for item in directory.iterdir():
        if item.is_dir():
            try:
                rmdir(item)
            except NotADirectoryError:
                # NotADirectoryError: [Errno 20] Not a directory:
                # '.../tmp/latest' -- a symlink (?) created by foldseek
                item.unlink()
        else:
            item.unlink()
    directory.rmdir()


'''
OMG getting this right took me forever, so some more comments on how to align
two protein structures.

Foldseek follows the same format as TM-align; checking both programs on the same
input results in nearly identical results -- examples from:

https://github.com/steineggerlab/foldseek/tree/master/data

TMalign d1asha.pdb d1jl7a.pdb -m matrix.txt

 **************************************************************************
 *                        TM-align (Version 20170708)                     *
 * An algorithm for protein structure alignment and comparison            *
 * Based on statistics:                                                   *
 *       0.0 < TM-score < 0.30, random structural similarity              *
 *       0.5 < TM-score < 1.00, in about the same fold                    *
 * Reference: Y Zhang and J Skolnick, Nucl Acids Res 33, 2302-9 (2005)    *
 * Please email your comments and suggestions to: zhng@umich.edu          *
 **************************************************************************

Name of Chain_1: d1asha.pdb
Name of Chain_2: d1jl7a.pdb
Length of Chain_1:  147 residues
Length of Chain_2:  147 residues

Aligned length=  138, RMSD=   2.54, Seq_ID=n_identical/n_aligned= 0.138
TM-score= 0.75831 (if normalized by length of Chain_1)
TM-score= 0.75831 (if normalized by length of Chain_2)
(You should use TM-score normalized by length of the reference protein)

(":" denotes aligned residue pairs of d < 5.0 A, "." denotes other aligned residues)
--ANKTRELCMKSLEHAKVDTSNEARQDGIDLYKHMFENYPPLRKYFKSREEYTAEDVQNDPFFAKQGQKILLACHVLCATYDDRETFNAYTRELLDRHARDH-VHMPPEVWTDFWKLFEEYLGKKTT--LDEPTKQAWHEIGREFAKE-IN--K-
  ::::::::::::::::.. . ..::::::::::::::::::::::: :: ::  ::  : ::::::::::::::::::::::::::::::::::::::::. ::::::::::::::::::::::::  ::::::::::::::::::: ::  .
GLSAAQRQVVASTWKDIAGA-D-NGAGVGKECLSKFISAHPEMAAVFG-FS-GA--SD--P-GVAELGAKVLAQIGVAVSHLGDEGKMVAEMKAVGVRHKGYGNKHIKAEYFEPLGASLLSAMEHRIGGKMNAAAKDAWAAAYGDISGALISGLQS


The (rotation, translation) matrix looks like this:

cat matrix.txt
 -------- Rotation matrix to rotate Chain_1 to Chain_2 ------
 m          t(m)         u(m,1)         u(m,2)         u(m,3)
 1    -20.5720310015   0.3041357799   0.6084949658   0.7329633715
 2    -23.4973527840   0.9222606725   0.0046469628  -0.3865406287
 3     42.2311133600  -0.2386140802   0.7935441275  -0.5597776687
 Code for rotating Chain_1 from (x,y,z) to (X,Y,Z):
    do i=1,L
      X(i)=t(1)+u(1,1)*x(i)+u(1,2)*y(i)+u(1,3)*z(i)
      Y(i)=t(2)+u(2,1)*x(i)+u(2,2)*y(i)+u(2,3)*z(i)
      Z(i)=t(3)+u(3,1)*x(i)+u(3,2)*y(i)+u(3,3)*z(i)
    enddo


Note how it says "Rotation matrix to rotate Chain_1 to Chain_2"; which protein
is projected onto which is important, bc/ we have to apply the rotation matrix
to __CHAIN 1__ here.

Applying the transformation is relatively straightforward:

structure.transform?? self.coord = np.dot(self.coord, rot) + tran

Or

def manual(coord, rot, tra):
    x, y, z = coord
    X = tra[0] + x*rot[0, 0] + y*rot[0, 1] + z*rot[0, 2] 
    Y = tra[1] + x*rot[1, 0] + y*rot[1, 1] + z*rot[1, 2] 
    Z = tra[2] + x*rot[2, 0] + y*rot[2, 1] + z*rot[2, 2]
    return X, Y, Z

Note that we use the PDB method .transform() and thus have to tanspose .T the
matrix; the align_structures() fn will return the transposed matrix.
'''
def align_structures(query: Structure, target: Structure, mode=0, minscore=0.5):
    '''
    PyMOL just wraps programs as well:

    - https://pymolwiki.org/index.php/TMalign
    - https://raw.githubusercontent.com/Pymol-Scripts/Pymol-script-repo/master/tmalign.py

    https://janakiev.com/blog/python-shell-commands/

    https://stackoverflow.com/questions/17742789/running-multiple-bash-commands-with-subprocess

    conda create -n foldseek -c conda-forge -c bioconda foldseek
    conda activate foldseek
    # wget https://mmseqs.com/foldseek/foldseek-osx-universal.tar.gz; tar xvzf foldseek-osx-universal.tar.gz; export PATH=$(pwd)/foldseek/bin/:$PATH

    To superimpose two structures we have to first find out how to transform 
    one into the other:

    https://github.com/steineggerlab/foldseek/issues/13

    mode:

    --cov-mode INT

    0: coverage of query and target
    1: coverage of target
    2: coverage of query
    3: target seq. length has to be at least x% of query length
    4: query seq. length has to be at least x% of target length
    5: short seq. needs to be at least x% of the other seq. length [0]

    --tmscore-threshold [0.5]
    accept alignments with a tmsore > thr [0.0,1.0] [0.500]
    '''
    assert is_tool('foldseek'), '"foldseek" appears not to be installed'
    tmp = tempfile.TemporaryDirectory()
    p = tmp.name
    print(f'Aligning {Path(target).name} to {Path(query).name}')

    steps = [
        f'foldseek createdb {target} {p}/targetDB',
        f'foldseek createdb {query} {p}/queryDB',
        f'foldseek search {p}/queryDB {p}/targetDB {p}/aln {p}/tmp -a --cov-mode {mode} --tmscore-threshold {minscore}',
        f'foldseek aln2tmscore {p}/queryDB {p}/targetDB {p}/aln {p}/aln_tmscore',
        f'foldseek createtsv {p}/queryDB {p}/targetDB {p}/aln_tmscore {p}/aln_tmscore.tsv'
    ]

    command = '; '.join(steps)
    log = subprocess.run(command, capture_output=True, shell=True)
    assert log.returncode == 0, log.stderr

    with open(f'{p}/aln_tmscore.tsv', 'r') as file:
        qry, rest = next(file).strip().split('\t')
        ref, score, *rest = rest.split(' ')
        rest = [float(i) for i in rest]
        translation = rest[:3]
        rotation = list(chunks(rest[3:], 3))
        # print(rotation, translation)
    tmp.cleanup()

    # Transpose rotation matrix!
    return round(float(score), 4), np.array(rotation).T, np.array(translation)


def transform_(structure: Structure, translation: np.ndarray, rotation: np.ndarray):
    '''
    DEPRECATED, can use "structure.transform(rotation, translation)"
    '''
    for i in structure.get_atoms():
        i.transform(rotation, translation)
    return None


def entropy(labels, base=None):
    '''
    Computes entropy of label distribution.

    https://stackoverflow.com/questions/15450192/fastest-way-to-compute-entropy-in-python

    nextstrain uses entropy:

    https://docs.nextstrain.org/en/latest/guides/share/download-data.html#diversity-entropy-data
    '''

    n_labels = len(labels)

    if n_labels <= 1:
      return 0

    value,counts = np.unique(labels, return_counts=True)
    probs = counts / n_labels
    n_classes = np.count_nonzero(probs)

    if n_classes <= 1:
      return 0

    ent = 0.

    # Compute entropy
    base = e if base is None else base
    for i in probs:
      ent -= i * log(i, base)

    return ent


def mean_pairwise_similarity(labels):
    '''
    Geneious calculates "Identity" as "Mean pairwise identity over all pairs
    in the column." (hover over "Identity" after loading a MSA).

    See also discussions here:

    - https://www.biostars.org/p/3856/
    - https://www.biostars.org/p/5067/
    - https://www.biostars.org/p/223824/

    For en extensive review:

    - https://onlinelibrary.wiley.com/doi/abs/10.1002/prot.10146
    '''
    l = [1 if i == j else 0 for i, j in combinations(labels, 2)]
    return round(sum(l) / len(l), 4)


def search_domains(fold, hmms, cpus=8):
    # hmmsearch -A aln.stk --cpu 8 --cut_ga --tblout 1AAY.tsv --domtblout 1AAY.dom.tsv ../Pfam-A.hmm rcsb_pdb_1AAY.fasta > result.txt

    assert is_tool('hmmsearch'), '"hmmsearch" appears not to be installed'

    tmp = tempfile.TemporaryDirectory()
    p = tmp.name

    with open(f'{p}/query.faa', 'w+') as out:
        out.write(f'>query\n{fold.sequence}\n')

    steps = [
        f'hmmsearch --cpu {cpus} --cut_ga {hmms} {p}/query.faa > {p}/found.txt'
    ]
    command = '; '.join(steps)
    log = subprocess.run(command, capture_output=True, shell=True)
    assert log.returncode == 0, log.stderr
    found = HMMERStandardOutput(f'{p}/found.txt')
    tmp.cleanup()
    return found.dom_hits


def filter_aa(structure):
    residues = []
    aa = 'ARNDCQEGHILKMFPSTWYV'
    for res in structure.get_residues():
        x = res.get_resname()
        x = x[0] + x[1:].lower()  # ALA > Ala
        if SeqUtils.IUPACData.protein_letters_3to1[x] in aa:
            residues.append(res)
    return residues


def flatten(d, parent_key='', sep='_', expected_track_length=None):
    '''
    - https://stackoverflow.com/questions/6027558/flatten-nested-dictionaries-compressing-keys
    - https://stackoverflow.com/questions/12555323/how-to-add-a-new-column-to-an-existing-dataframe

    df = pd.DataFrame.from_dict(flatten(d, expected_track_length=len(model)))
    df = df.assign(selection=pd.Series(anno).values)
    '''
    items = []
    for k, v in d.items():            
        new_key = parent_key + sep + k if parent_key else k
        if isinstance(v, collections.MutableMapping):
            items.extend(flatten(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    
    if not expected_track_length:
        return dict(items)

    else:
        cleaned = {}
        for k, v in dict(items).items():
            try:
                if len(v) == expected_track_length and type(v) != str:
                    cleaned[k] = v
            except TypeError:
                # TypeError: object of type 'float' has no len()
                continue
        return cleaned
                

def is_tool(name):
    '''
    Check whether <name> is on PATH and marked as executable.

    https://stackoverflow.com/questions/11210104/check-if-a-program-exists-from-a-python-script
    '''

    from shutil import which
    return which(name) is not None
