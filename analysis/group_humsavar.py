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
humsavar_benign_df = humsavar_tsv_df[humsavar_tsv_df['category'] == 'LB/B']
humsavar_pathogenic_df = humsavar_tsv_df[humsavar_tsv_df['category'] == 'LP/P']
   
humsavar_benign_df.to_csv(HUMSAVAR_BENIGN, sep='\t', index=False)
humsavar_pathogenic_df.to_csv(HUMSAVAR_PATHOGENIC, sep='\t', index=False)
   