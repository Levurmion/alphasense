import os
import sys
import argparse
import re
import math
import pandas as pd
from utilities import readFile_as_generator

# establish files path

parser = argparse.ArgumentParser(
    prog='filter_vep_output.py',
    description='filters Ensembl VEP output',
)

parser.add_argument('-f', '--file', help='path to the file relative to the PYTHONPATH dir', required=True)
parser.add_argument('-o', '--output', help='path and name of the output file relative to the PYTHONPATH dir', required=True)

try:
    args = parser.parse_args()
except:
    parser.print_help()
    sys.exit(1)

INPUT_FILENAME = args.file.split('/')[-1]
FILE_PATH = args.file
OUTFILE_PATH = args.output

FILE_CODE = re.search(r'uniprot_human_variants_HGVS_(.+)\.csv',INPUT_FILENAME).group(1)
print(FILE_CODE)

infile_df = pd.read_csv(FILE_PATH, sep='\t', header=0)

outfile_columns = [
   'uniprot',
   'protein_position',
   'WT',
   'Mut',
   'CLIN_SIG',
   'gnomADg_AF'
]

outfile_benign = []
outfile_damaging = []

print(infile_df)

clinicalSigRegex = r'(uncertain|conflicting)'
benignRegex = r'(benign)'
pathogenicRegex = r'(pathogenic)'
consequenceRegex = r'(missense_variant)'

SEEN_MUTATIONS = {}

for idx, row in infile_df.iterrows():
   
   consequence = row['Consequence']
   AF = row['gnomADg_AF'] if row['gnomADg_AF'] == '-' else float(row['gnomADg_AF'])
   clinSig = row['CLIN_SIG']
   uniprot = re.sub(r'\.\s*([^\s.]+)','',row['SWISSPROT'])
   
   mutation = row['Amino_acids'].split('/')
   
   SEEN_KEY = uniprot + row['Amino_acids']
   
   row_dict = {
      'uniprot': uniprot,
      'protein_position': row['Protein_position'],
      'WT': mutation[0] if len(mutation) == 2 else '-',
      'Mut': mutation[1] if len(mutation) == 2 else '-',
      'CLIN_SIG': clinSig,
      'gnomADg_AF': AF
   }
   
   try:
      seen = SEEN_MUTATIONS[SEEN_KEY]
      
   except KeyError:
      SEEN_MUTATIONS[SEEN_KEY] = True
      
      if bool(re.search(consequenceRegex, consequence)) == True and uniprot != '-':
         
         if bool(re.search(clinicalSigRegex, clinSig)) == False:
            
            if bool(re.search(benignRegex, clinSig)) == True:
               outfile_benign.append(row_dict)
               
            elif bool(re.search(pathogenicRegex, clinSig)) == True:
               outfile_damaging.append(row_dict)
               
            else:
               if AF == '-':
                  pass
               
               elif AF >= 0.01:
                  outfile_benign.append(row_dict)
               
               elif AF < 0.01:
                  outfile_damaging.append(row_dict)
            
         else:
            pass
         
      else:
         pass
   
outfile_benign_df = pd.DataFrame(outfile_benign, columns=outfile_columns)
outfile_damaging_df = pd.DataFrame(outfile_damaging, columns=outfile_columns)

benign_name = 'benign_mutations_' + FILE_CODE + '.csv'
damaging_name = 'pathogenic_mutation_' + FILE_CODE + '.csv'

benign_path = os.path.join(OUTFILE_PATH, benign_name)
damaging_path = os.path.join(OUTFILE_PATH, damaging_name)

outfile_benign_df.to_csv(benign_path,sep='\t',index=False)
outfile_damaging_df.to_csv(damaging_path,sep='\t',index=False)