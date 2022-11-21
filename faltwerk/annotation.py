import hashlib

import screed


def hash_aa_sequence(seq, fn=hashlib.md5):
    return fn(seq.encode('utf-8')).hexdigest()


def annotate_from_fasta(cx, path):
    '''
    Annotate the amino acid sequences in a Complex object using the headers
    from a fasta file. Matching is done using sequence hashes.

    Examples:

    >>> from faltwerk.annotation import annotate_from_fasta
    >>> annotate_from_fasta(cx, 'proteins.faa')
    {'B': 'foo', 'C': 'bar'}
    '''
    anno, result = {}, {}
    with screed.open(path) as file:
        for line in file:
            hsh = hash_aa_sequence(line.sequence)
            anno[hsh] = line.name

    for ch, seq in cx.sequences.items():
        hsh = hash_aa_sequence(seq)
        result[ch] = anno.get(hsh)

    return result
