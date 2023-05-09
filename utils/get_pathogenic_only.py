import os
import sys
import argparse
import pandas as pd

ROOT_DIR = os.path.join(os.environ.get('PYTHONPATH'))

PATHOGENIC_PATH = os.path.join(ROOT_DIR,'synced_files/uniprot_missense_variants/pathogenic_mutations_all.csv')
OUTFILE_PATH = os.path.join(ROOT_DIR,'synced_files/uniprot_missense_variants/pathogenic_mutations_annotated.csv')

pathogenic_df = pd.read_csv(PATHOGENIC_PATH, sep='\t', header=0)
pathogenic_df_headers = pathogenic_df.columns
outfile_rows = []
         
for idx, row in pathogenic_df.iterrows():
   
   if row['CLIN_SIG'] != '-':
      outfile_rows.append(row)

outfile_df = pd.DataFrame(outfile_rows,columns=pathogenic_df_headers)

outfile_df.to_csv(OUTFILE_PATH, sep='\t', index=False)

print(outfile_df)