import os
import sys


def pad_line(line):
    '''
    https://github.com/haddocking/pdb-tools/blob/master/pdbtools/pdb_reres.py#L107

    Helper function to pad line to 80 characters in case it is shorter
    '''
    size_of_line = len(line)
    if size_of_line < 80:
        padding = 80 - size_of_line + 1
        line = line.strip('\n') + ' ' * padding + '\n'
    return line[:81]  # 80 + newline character


def reindex_pdb(fhandle, starting_resid):
    '''
    https://github.com/haddocking/pdb-tools/blob/master/pdbtools/pdb_reres.py#L116

    Reset the residue number column to start from a specific number.

    This function is a generator.

    Parameters
    ----------
    fhandle : a line-by-line iterator of the original PDB file.

    starting_resid : int
        The starting residue number.

    Yields
    ------
    str (line-by-line)
        The modified (or not) PDB line.
    '''
    _pad_line = pad_line
    prev_resid = None  # tracks chain and resid
    resid = starting_resid - 1  # account for first residue
    records = ('ATOM', 'HETATM', 'TER', 'ANISOU')
    for line in fhandle:
        line = _pad_line(line)
        if line.startswith('MODEL'):
            resid = starting_resid - 1  # account for first residue
            prev_resid = None  # tracks chain and resid
            yield line

        elif line.startswith(records):
            line_resuid = line[17:27]
            if line_resuid != prev_resid:
                prev_resid = line_resuid
                resid += 1
                if resid > 9999:
                    emsg = 'Cannot set residue number above 9999.\n'
                    sys.stderr.write(emsg)
                    sys.exit(1)

            yield line[:22] + str(resid).rjust(4) + line[26:]

        else:
            yield line

