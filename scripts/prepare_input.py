#! /usr/bin/env python3


from pathlib import Path


p = Path('genomes')
files = list(p.glob('[!.]*.fna.gz'))  # skip hidden files
assert len(files) > 0


with open('input.csv', 'w+') as out:
    for i in files:
        out.write(f'{i.name.split(".")[0]},{i.absolute().__str__()}\n')
