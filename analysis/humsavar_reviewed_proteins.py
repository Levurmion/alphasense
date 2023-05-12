import os
import sys
import pandas as pd
from pickleconfig import PickleJar
from utils import readFile_as_generator

jar = PickleJar()

ROOT_DIR = os.path.join(os.environ.get('PYTHONPATH'))

UNIPROT_REV_PATH = os.path.join(ROOT_DIR, 'files/uniprot_human_reviewed_lt_2700aa.csv')
HUMSAVAR_TSV_PATH = os.path.join(ROOT_DIR, 'synced_files/uniprot_missense_variants/humsavar_pathogenic.tsv')
HUMSAVAR_REVIEWED = os.path.join(ROOT_DIR, 'synced_files/uniprot_missense_variants/humsavar_pathogenic_reviewed.tsv')

def generate_prot_dict():
   
   dict = {}
   
   for protein in readFile_as_generator(UNIPROT_REV_PATH):
      dict[protein.replace('\n','')] = True
   
   return dict

UNIPROT_DICT = jar.pickle(generate_prot_dict, 'uniprot_dict')

humsavar_tsv_df = pd.read_csv(HUMSAVAR_TSV_PATH, sep='\t')
humsavar_reviewed = []

for row in humsavar_tsv_df.itertuples():
   
   try:
      exists = UNIPROT_DICT[getattr(row, 'AC')]
      humsavar_reviewed.append(row)
   except KeyError:
      pass

humsavar_reviewed_df = pd.DataFrame(humsavar_reviewed).drop(columns=['Index'])

humsavar_reviewed_df.to_csv(HUMSAVAR_REVIEWED, sep='\t', index=False)