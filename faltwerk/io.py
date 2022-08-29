from collections import defaultdict
import gzip
from io import StringIO
import json
from pathlib import Path
import tempfile
try:
    from typing import Union
except ImportError:
    from typing_extensions import Union

from Bio import PDB, SeqUtils
from Bio.PDB import PDBIO, Structure
from Bio.PDB.PDBParser import PDBParser
from Bio import SeqIO
import pandas as pd
import requests
import screed

from faltwerk.utils import entropy, mean_pairwise_similarity


def save_pdb(structure: Structure, out: Union[str, Path, StringIO]) -> None:
    file = PDBIO()
    file.set_structure(structure)
    
    # Save to stream
    if type(out) == StringIO:
        file.save(out)
        return out

    # Save to file
    else:
        p = Path(out)
        if not p.parent.exists():
            p.parent.mkdir(parents=True)
        file.save(p.resolve().__str__())
        return None


def load_conserved(fp, ref=None, metric=mean_pairwise_similarity):
    '''
    If no reference sequence name is provided, assume the first sequence is
    the reference. Why do we even need to specify the reference? Bc/ in the MSA
    it can contain gaps, which we'll omit bc/ we want to be able to map the
    conservation values to the protein structure, which does not contain gaps
    and we assume is identical to the reference sequence.

    Available functions:

    - mean_pairwise_similarity
    - entropy
    '''
    variance = defaultdict(list)
    cnt, ix = 0, -1

    with screed.open(fp) as file:
        # If no reference is provided, assume that the first sequence is ref.
        if (not ref) and (cnt == 0):
            ix = cnt

        for i in file:
            if (ref == i.name) and (ix != 0):
                ix = cnt
            
            for pos, aa in enumerate(i.sequence):
                variance[pos].append(aa)
            cnt += 1
        
    l = []
    for pos, residues in variance.items():
        ref_aa = residues[ix]
        if not ref_aa == '-':
            l.append(metric(residues))

    return l


def is_gz_file(filepath):
    '''
    https://stackoverflow.com/questions/3703276/how-to-tell-if-a-file-is-gzip-compressed
    '''
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'


def read_sequence(fp: Union[str, Path]):
    # https://www.biostars.org/p/435629/
    with open(fp, 'r') as file:
        l = []
        for record in SeqIO.parse(file, 'pdb-atom'):
            l.append(record.seq)
        if len(l) > 1:
            print('Warning: More than one sequence found, returning first')
        return l[0].__str__()


def read_pdb(fp: Union[str, Path], name: str='x', strict: bool=True) -> Structure:
    '''
    Return structure AND sequence

    # https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ

    p = PDBParser()
    structure = p.get_structure("X", "pdb1fat.ent")
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    print(atom)
    '''
    fp = Path(fp)
    assert fp.exists()

    # print(is_gz_file(fp))
    if is_gz_file(fp):
        tmp = tempfile.NamedTemporaryFile()

        with gzip.open(fp) as file:
            file_content = file.read()
            tmp.write(file_content)
            fp = Path(tmp.name)

    # https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ
    pdb_parser = PDB.PDBParser(QUIET=True, PERMISSIVE=0)
    structure = pdb_parser.get_structure(name, str(fp))
    sequence = read_sequence(fp)

    counts = {
       'models': 0,
       'chains': 0,
    }

    for model in structure:
        counts['models'] += 1
        for chain in model:
            counts['chains'] += 1

    if strict:
        assert sum(counts.values()) == 2, 'More than one model or chain present'

    try:
        tmp.close()
    except NameError:
        pass

    return structure, sequence


def load_bfactor_column(fp):
    '''
    Load annotation data stored in the bfactor column of a .pdb file
    '''
    parser = PDBParser()
    structure = parser.get_structure('', fp)
    d = {}

    for res in structure.get_residues():
        for atom in res:
            # ALA > Ala
            x = atom.parent.resname[0] + atom.parent.resname[1:].lower()  
            try:
                aa = SeqUtils.IUPACData.protein_letters_3to1[x]  
                # same value for all atoms in a resudue
            except KeyError:
                continue

            num = atom.full_id[3][1]
            d[num] = [atom.bfactor, aa]
    
    return [i for i, j in d.values()]


def parse_hyphy(fp, method='meme', direction='positive', skip=[]):
    '''
    TODO: Just count the number of "MEME" or "FUBAR" in the results file to
    infer the program that was used.

    skip .. use e.g. to mask gaps in an alignment, needs to be some form of
    binary iterator (eg [0, 0, 0, 1, 0, 0, 0, ...])
    '''
    with open(fp, 'r') as file:
        d = json.load(file)

    if method == 'meme':
        # p-value
        assert direction == 'positive', 'The MEME method does not estimate negative selection'
        scores = [i[6] for i in d['MLE']['content']['0']]
        

    elif method == 'fubar':
        # posterior probability
        neg, pos = zip(*[i[3:5] for i in d['MLE']['content']['0']])
        if direction == 'positive':
            scores = pos
        elif direction == 'negative':
            scores = neg
        else:
            raise ValueError('Direction must be positive or negative.')

    else:
        raise ValueError('Method not implemented')

    if skip:
        assert len(skip) == len(scores), 'Mask has wrong dimensions'
        scores = [p for p, sk in zip(scores, skip) if not sk]

    return scores


def get_alphafold2_model_from_ena(ID=None, fp=None):
    # eg AF-Q09428-F1
    url = f'https://alphafold.ebi.ac.uk/files/{ID}-model_v3.pdb'
    r = requests.get(url)
    assert r.status_code == 200, 'Download failed'
    pdb = r.text
    if not fp:
        fp = f'{ID}.pdb'
    with open(fp, 'w+') as out:
        out.write(pdb)
    return None
