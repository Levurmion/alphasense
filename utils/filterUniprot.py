import os
import sys

print('starting filtering procedure...')

# establish files path
FILES_PATH = os.path.join(os.environ.get('PYTHONPATH'),'files')

import pandas as pd

if len(sys.argv) != 4:
    print('Usage: python filterUniprot <humanVariants filename> <reviewedProteins filename> <output filename>')
    exit(1)

humanVariants_path = os.path.join(FILES_PATH, sys.argv[1])
reviewedProteins_path = os.path.join(FILES_PATH, sys.argv[2])
output_path = os.path.join(FILES_PATH, sys.argv[3])

humanVariants_df = pd.read_csv(humanVariants_path, sep='\t')
reviewedProteins_df = pd.read_csv(reviewedProteins_path, sep='\t')

import re

# assemble reviewed proteins into a dictionary
REVIEWED_PROTEINS = {}

for idx, protein in reviewedProteins_df.loc[:, 'Entry'].items():
    REVIEWED_PROTEINS[protein] = idx

REVIEWED_PROTEINS_COLS = reviewedProteins_df.columns

FILTERED_DF_COL = [
    'Variant AA Change',
    'Source DB ID',
    'Clinical Significance',
    'Evidence',
]

HEADERS = '\t'.join(list(REVIEWED_PROTEINS_COLS) + FILTERED_DF_COL)

clinicalSigRegex = r'\b(conflicting|uncertain|Conflicting|Uncertain)\b'

originalFileLength = len(humanVariants_df.index)
quantPassedFilter = 0

with open(output_path, 'w') as outfile:

   outfile.write(HEADERS + '\n')

   currentLine = 0

   for idx, row in humanVariants_df.iterrows():

      currentLine += 1

      if currentLine % 100000 == 0:
         print(f'processing line: {currentLine} \r', end='')

      try:
         reviewedProteinsIdx = REVIEWED_PROTEINS[row['AC']]

         if row['Consequence Type'] == 'missense variant' and bool(re.search(clinicalSigRegex, row['Clinical Significance'])) == False:
               
            lineEntries = []

            for col1 in REVIEWED_PROTEINS_COLS:
               lineEntries.append(str(reviewedProteins_df.loc[reviewedProteinsIdx, col1]))

            for col2 in FILTERED_DF_COL:
               lineEntries.append(str(row[col2]))
         
            lineText = '\t'.join(lineEntries)

            outfile.write(lineText + '\n')

            quantPassedFilter += 1
               
      except KeyError:
         pass

print('original file length: ', originalFileLength)
print('filtered dataset length: ', quantPassedFilter)
