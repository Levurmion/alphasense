import os
import sys
import pandas as pd

# establish files path
ROOT_DIR = os.environ.get('PYTHONPATH')

if len(sys.argv) != 3:
    print('Usage: python get_HGVS-uniprot.py <uniprot variants filepath> <output filepath>')
    exit(1)

INFILEPATH = os.path.join(ROOT_DIR, sys.argv[1])
OUTFILEPATH = os.path.join(ROOT_DIR, sys.argv[2])

print('loading file into dataframe')
UNIPROT_VAR_FILE = pd.read_csv(INFILEPATH, sep='\t', header=0)

with open(OUTFILEPATH, 'w') as outfile:

    headerLine = 'Chromosome Coordinate\tAC\tClinical Significance\tVariant AA Change\n'
    outfile.write(headerLine)

    lineCounter = 0

    for idx, row in UNIPROT_VAR_FILE.iterrows():
        HGVS = row['Chromosome Coordinate']
        uniprot = row['AC']
        clinsig = row['Clinical Significance']
        aaMutation = row['Variant AA Change']
        outfileRow = f'{HGVS}\t{uniprot}\t{clinsig}\t{aaMutation}\n'
        outfile.write(outfileRow)

        lineCounter += 1

        if lineCounter % 100000 == 0:
            print(f'processing line: {lineCounter}\r', end='')
        


