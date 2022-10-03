# conda activate txmap
# see journal, 2022-08-26
from collections import defaultdict
import os
import re

import pandas as pd
from tqdm import tqdm

import hgvs
from hgvs.assemblymapper import AssemblyMapper
from cdot.hgvs.dataproviders import JSONDataProvider


'''
cdot transcripts:
https://github.com/SACGF/cdot#q-where-can-i-download-the-jsongz-files

wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

wget https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.0/MANE.GRCh38.v1.0.summary.txt.gz

refseq_prot_to_alphafold.csv is in the unpacked bulk download from
"MANE Select dataset" here:
https://alphafold.ebi.ac.uk/download
'''
os.environ['HGVS_SEQREPO_DIR'] = '.../seqrepo/2021-01-29/'
fp = '.../cdot/cdot-0.2.7.refseq.grch37_grch38.json.gz'
fp_variant_summary = 'variant_summary.txt.gz'
fp_mane_summary = 'MANE.GRCh38.v1.0.summary.txt'
fp_refseq_alphafold = 'refseq_prot_to_alphafold.csv'


def parse_hgvs_code(x):
    target = x.split(' ')[0]
    try:
        # Get HGVS c. syntax from eg
        # NM_014855.3(AP5Z1):c.1413_1426del (p.Leu473fs)
        # > NM_014855.3:c.1413_1426del 
        m = re.match('(NM_.*?)\(.*?\)(.*)', target)
        return m.group(1) + m.group(2)
    except AttributeError:
        # odd format
        return None


hdp = JSONDataProvider([fp])

am = AssemblyMapper(
    hdp,
    assembly_name='GRCh38',
    alt_aln_method='splign',
    replace_reference=True)


df = pd.read_csv(fp_variant_summary, compression='gzip', sep='\t', low_memory=False)


mane = pd.read_csv(fp_mane_summary, sep='\t')
d = {}
for i in mane.itertuples():
    d[i.RefSeq_nuc] = i.RefSeq_prot


significance = set([
    'Pathogenic',
    'Likely pathogenic',
    'Pathogenic/Likely pathogenic',
    ])

types = set([
    'single nucleotide variant',
    ])


no_variant, no_mane, g37, xy = 0, 0, 0, 0
# txmap = defaultdict(list)
map_tx, map_aa = {}, {}
codes = {}

for i in tqdm(df.itertuples()):
    if not i.ClinicalSignificance in significance or i.Type not in types:
        continue
    #  109269 .. constrained (likely pathogenic SNVs)
    # 1450297 .. unconstrained
    try:
        tx = re.match('(.*?)\(.*', i.Name).group(1)
        try:
            if i.Assembly == 'GRCh38':
                if not i.VariationID in map_aa:
                    
                    map_aa[i.VariationID] = d[tx]  # causes key error, do 1st
                    map_tx[i.VariationID] = tx

                    code = parse_hgvs_code(i.Name)
                    if code:
                        codes[i.VariationID] = code
                    
                else:
                    # Some MANE tx are on X and Y chromosome
                    xy += 1

            else:
                g37 += 1

            # if i.VariationID == 517:
            #     print(i)
        except KeyError:
            no_mane += 1
            continue

    except AttributeError:
        no_variant += 1
        continue
# 18342 variants don't pass, 121852 have no mane
# (no_variant + no_mane) / 3042140 < 0.05


af = pd.read_csv(fp_refseq_alphafold)


pos = defaultdict(set)
notfound = 0
map_af = {}
hp = hgvs.parser.Parser()
for ID, hgvs_c in tqdm(codes.items()):
    var_c = hp.parse_hgvs_variant(hgvs_c)
    # var_c = am.g_to_c(var_g, map_tx[ID])
    # l.append(var_c)

    # We could now do var_p = am.c_to_p(var_c, 'protein_ID'), but:
    # https://github.com/SACGF/cdot/issues/17
    # https://github.com/biocommons/hgvs/issues/634#issuecomment-1206763697
    # https://github.com/SACGF/cdot/issues/14
    # So we go manual for now.
    base = var_c.posedit.pos.start.base
    offset = var_c.posedit.pos.start.offset 

    if not offset:
        aa = map_aa[ID]
        try:
            model = af[af['refseq_prot'] == aa]['alphafold'].item()
            map_af[aa] = model
            p = int(base / 3)
            if p >= 0:
                pos[model].add(p)
            # 0-based, some values are negative, TODO: investigate:
            # AF-Q8IXV7-F1 -52
        except ValueError:
            # not found in AF2 models
            notfound += 1
            continue


with open('clinvar.tsv', 'w+') as out:
    out.write('ID\tpositions\n')
    for k, v in pos.items():
        out.write(f'{k}\t{",".join(str(i) for i in sorted(v))}\n')


l, nokey = [], []
for k, v in map_tx.items():
    aa = map_aa[k]
    try:
        model = map_af[aa]
        l.append([k, v, aa, model])
    except KeyError:
        nokey.append(aa)
        continue
result = pd.DataFrame(l)
result.columns = 'clinvar mane mane_aa alphafold'.split(' ')
result.to_csv('mane_af2.csv', index=False)
# [84209 rows x 4 columns] -- about 0.8 of clinvar records mapped -- rest? TODO
'''
       clinvar         mane      mane_aa         alphafold
0            5  NM_017547.4  NP_060017.1      AF-Q96CU9-F1
1            6  NM_017547.4  NP_060017.1      AF-Q96CU9-F1


TODO:

There are 16k+ models but we only find 4k?

In [31]: len(map_af)
Out[31]: 4086

In [32]: len(pos)
Out[32]: 4071

In [34]: notfound
Out[34]: 18859
'''



