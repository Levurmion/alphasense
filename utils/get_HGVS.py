import os
import sys
import pandas as pd
import re

# establish files path
FILES_PATH = os.path.join(os.environ.get('PYTHONPATH'),'files')

humanVariants_path = os.path.join(FILES_PATH, sys.argv[1])
reviewedProteins_path = os.path.join(FILES_PATH, sys.argv[2])
output_path = os.path.join(FILES_PATH, sys.argv[3])

print('loading human variants file into dataframe...')
humanVariants_df = pd.read_csv(humanVariants_path, sep='\t', header=0)
print(humanVariants_df.columns)
print('loading reviewed proteins list into dataframe...')
reviewedProteins_df = pd.read_csv(reviewedProteins_path, sep='\t', header=0)

REVIEWED_PROTEINS = {}

for idx, protein in reviewedProteins_df.loc[:, 'Entry'].items():
    REVIEWED_PROTEINS[protein] = idx

clinicalSigRegex = r'\b(conflicting|uncertain|Conflicting|Uncertain)\b'
indelRegex = r'\b(delins)\b'

line = 0

with open(output_path, 'w') as outfile:
    
    for idx, row in humanVariants_df.iterrows():
        
        if row['Consequence Type'] == 'missense variant' and bool(re.search(clinicalSigRegex, row['Clinical Significance'])) == False:
            
            try:
                proteinReviewed = REVIEWED_PROTEINS[row['AC']]
                if bool(re.search(indelRegex, row['Chromosome Coordinate'])) == False:
                    outfile.write(row['Chromosome Coordinate'] + '\n')
                else:
                    pass
            except KeyError:
                pass
        
        line += 1
      
        if line%100000 == 0:
            print('processing line: ', line, '\r', end='')
    