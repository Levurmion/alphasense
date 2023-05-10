import os
import sys
import pandas as pd

# establish files path
ROOT_DIR = os.environ.get('PYTHONPATH')

if len(sys.argv) != 3:
    print('Usage: python get_short_proteins.py <reviewed uniprot filepath> <output filepath>')
    exit(1)

UNIPROT_FILEPATH = os.path.join(ROOT_DIR, sys.argv[1])
OUTPUT_FILEPATH = os.path.join(ROOT_DIR, sys.argv[2])

uniprot_df = pd.read_csv(UNIPROT_FILEPATH, sep='\t', header=0)

proteins_gt2700 = []

with open(OUTPUT_FILEPATH, 'w') as outfile:

    outfile.write("Entry\n")

    for idx, row in uniprot_df.iterrows():

        if int(row['Length']) <= 2700:
            outfile.write(str(row['Entry']) + '\n')
        # else:
        #     proteins_gt2700.append(row)
    
    # outfile.write('---------------PROTEINS >2700aa---------------\n')
    
    # for protein in proteins_gt2700:
    #     uniprotId = protein['Entry']
    #     length = protein['Length']
    #     outfile.write(f'{uniprotId}\t{length}\n')
