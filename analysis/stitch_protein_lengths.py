import os
import sys
import pandas as pd
from betterpicklejar import *

PickleShelf('pickle_shelf')
jar = PickleJar()

ROOT_DIR = os.path.join(os.environ.get('PYTHONPATH'))

UNIPROT_REF_PATH = os.path.join(ROOT_DIR,'files/uniprot_human_reviewed.csv')
INFILE_PATH = os.path.join(ROOT_DIR,'synced_files/uniprot_missense_variants/pathogenic_mutations_all.csv')
OUTFILE_PATH = os.path.join(ROOT_DIR,'synced_files/uniprot_missense_variants/pathogenic_mutations_with_lengths.csv')

uniprot_ref_df = pd.read_csv(UNIPROT_REF_PATH, sep='\t', header=0)
infile_df = pd.read_csv(INFILE_PATH, sep='\t', header=0)
outfile_rows = []

headers = list(infile_df.columns)
headers.append('length')


ENTRY_LENGTHS = {}

for idx, row in uniprot_ref_df.iterrows():
   
   entry = row['Entry']
   length = int(row['Length'])
   
   if length <= 2700:
      ENTRY_LENGTHS[entry] = length
      

for idx, row in infile_df.iterrows():

   try:
      length = pd.Series({'length': ENTRY_LENGTHS[row['uniprot']]})
      new_row = pd.concat([row,length], axis=0)
      outfile_rows.append(new_row)
   except KeyError:
      pass

outfile_df = pd.DataFrame(outfile_rows, columns=headers)
outfile_df.to_csv(OUTFILE_PATH ,index=False, sep='\t')
