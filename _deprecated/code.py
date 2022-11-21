
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