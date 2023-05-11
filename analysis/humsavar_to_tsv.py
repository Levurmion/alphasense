import os
import sys
import pandas as pd
from pickleconfig import PickleJar
from utils import readFile_as_generator
import re

jar = PickleJar()

ROOT_DIR = os.path.join(os.environ.get('PYTHONPATH'))

HUMSAVAR_PATH = os.path.join(ROOT_DIR, 'files/humsavar.txt')
HUMSAVAR_TSV = os.path.join(ROOT_DIR, 'files/humsavar.tsv')

HEADERS = ''
ROWS = []

for idx, line in enumerate(readFile_as_generator(HUMSAVAR_PATH)):
   
   if idx == 0:
      HEADERS = '\t'.join([head.strip() for head in re.split(r'\s{2,}',line)])
   else:
      gene_name = line[0:9].strip()
      AC = line[10:21].strip()
      FTId = line[21:33].strip()
      change = line[33:48].strip()
      category = line[48:57].strip()
      dbSNP = line[57:72].strip()
      disease_name = line[72::].strip()
      print(gene_name, AC, FTId, change, category, dbSNP, disease_name)
      ROWS.append('\t'.join([gene_name, AC, FTId, change, category, dbSNP, disease_name]) + '\n')
      

with open(HUMSAVAR_TSV, 'w') as outfile:
   
   outfile.write(HEADERS + '\n')
   
   for row in ROWS:
      outfile.write(row)
      