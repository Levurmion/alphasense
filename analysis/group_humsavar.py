import os
import sys
import pandas as pd
from pickleconfig import PickleJar
from utils import AA_DICT_LTS
import re

jar = PickleJar()

ROOT_DIR = os.path.join(os.environ.get('PYTHONPATH'))

HUMSAVAR_TSV = os.path.join(ROOT_DIR, 'files/humsavar.tsv')
HUMSAVAR_BENIGN = os.path.join(ROOT_DIR, 'files/humsavar_benign.tsv')
HUMSAVAR_PATHOGENIC = os.path.join(ROOT_DIR, 'files/humsavar_pathogenic.tsv')

humsavar_tsv_df = pd.read_csv(HUMSAVAR_TSV, sep='\t')
humsavar_benign = []
humsavar_pathogenic = []

for idx, row in humsavar_tsv_df.iterrows():
   
   category = row['category']
   
   if category == 'LB/B':
      humsavar_benign.append(row.drop(columns=['Index']))
   elif category == 'LP/P':
      humsavar_pathogenic.append(row.drop(columns=['Index']))

humsavar_benign_df = pd.DataFrame(humsavar_benign)
humsavar_pathogenic_df = pd.DataFrame(humsavar_pathogenic)
   
humsavar_benign_df.to_csv(HUMSAVAR_BENIGN, sep='\t', index=False)
humsavar_pathogenic_df.to_csv(HUMSAVAR_PATHOGENIC, sep='\t', index=False)
   