# Solvent accessibility


# https://biopython.org/docs/1.75/api/Bio.PDB.DSSP.html
# https://anaconda.org/salilab/dssp

# p = PDBParser()
# structure = p.get_structure('', args.model)
# dssp = DSSP(structure[0], args.model, dssp='mkdssp')
# asa = [i[3] for i in dssp.property_list]

from Bio.PDB.DSSP import DSSP


def solvent_access(fold):
    m = fold.structure[0]
    dssp = DSSP(m, fold.path.__str__(), dssp='mkdssp')
    asa = [i[3] for i in dssp.property_list]
    return asa


