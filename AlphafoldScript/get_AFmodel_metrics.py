from AlphafoldModel import AlphafoldModel
import os
import sys
import pandas as pd
import argparse
import numpy as np
from betterpicklejar import *

PickleShelf('../pickle_shelf')

jar = PickleJar()

# get command line arguments
parser = argparse.ArgumentParser(
   prog='get_AFmodel_metrics.py',
   description='Returns a CSV of the average PAE and plDDT for a given list of variants on Alphafold models. By default, measures the average PAE between all neighbouring residue pairs.'
)

# arguments list
parser.add_argument('--pdb_dir', help='Path to your Alphafold model directory (PDB format).')
parser.add_argument('--pae_dir', help='Path to your PAE matrix directory (JSON format).')
parser.add_argument('-i', '--input', help='Path to your list of missense variants. The file should be a TSV with columns: uniprot, cluster, WT, Mut, position.')
parser.add_argument('-o', '--output', help='Path and output filename.')
parser.add_argument('-r', '--radius', help='A comma-separated argument string of the list of radii to measure the PAE and plDDT scores. ("5.0,5.5,6.0")')
parser.add_argument('--pae_query_only', help='Flag to calculate the average PAE using residue pairs involving the query only. Default=False.')
parser.add_argument('--plddt_window', help='A comma-separated argument string of the window size (must be odd numbers) to calculate plDDT scores. ("3,5,7")')
parser.add_argument('--scrap_ori', help='Flag on whether to keep the original columns. Default= False.')
parser.add_argument('--pickle', help='Use pickle to save alphafold models for future runs.')

args = parser.parse_args()

PDB_PATH: str = args.pdb_dir
PAE_PATH: str = args.pae_dir
INPUT_PATH: str = args.input
OUTPUT_PATH: str = args.output
RADIUS: list = [5] if args.radius == None else [float(radius) for radius in args.radius.split(',')]
PAE_QUERY_ONLY: bool = True if args.pae_query_only != None else False
PLDDT_WINDOW = [float(radius) for radius in args.plddt_window.split(',')] if isinstance(args.plddt_window, str) else False
SCRAP_COL: bool = True if args.scrap_ori != None else False
pickle_name = False if args.pickle == None else args.pickle


# input check
if PLDDT_WINDOW != False:
   
   for window in PLDDT_WINDOW:
      if int(window) % 2 == 0:
         raise Exception('Window needs to be a positive integer odd number.')
         exit(1)
      else:
         pass

# script start


VARIANTS_DF: pd.DataFrame = pd.read_csv(INPUT_PATH, sep='\t')

UNIPROT_IDS = set(VARIANTS_DF['uniprot'])

print(f'Variants file has {len(UNIPROT_IDS)} unique UniProt IDs.')

PROBLEM_MODELS = []
PAE_NO_MATCH = []

# load all required Alphafold models

def load_AF_models():
   
   AF_MODELS: dict = {}
   
   for uniprot in UNIPROT_IDS:
      
      MODEL_PDB = f'AF-{uniprot}-F1-model_v4.pdb'
      MODEL_PAE = f'AF-{uniprot}-F1-predicted_aligned_error_v4.json'
      MODEL_PDB_PATH = os.path.join(PDB_PATH, MODEL_PDB)
      MODEL_PAE_PATH = os.path.join(PAE_PATH, MODEL_PAE)
      
      print(MODEL_PDB_PATH)
      print(MODEL_PAE_PATH)
      
      try:
         AF_MODELS[uniprot] = AlphafoldModel(MODEL_PDB_PATH, MODEL_PAE_PATH)
         print(f'Loading {uniprot} into memory... \r', end='')
      except FileNotFoundError:
         print(f'{uniprot} not found...')
         PROBLEM_MODELS.append(uniprot)
      except ValueError:
         print(f'{MODEL_PAE_PATH} not found...')
         PAE_NO_MATCH.append(uniprot)
   
   return AF_MODELS

AF_MODELS = None

if pickle_name == False:
   AF_MODELS = load_AF_models()
else:
   AF_MODELS = jar.pickle(load_AF_models, pickle_name=pickle_name)

# initialize output columns as dict
OUTPUT_COLS = {}

def gen_pae_colname(rad):
   return f'pae_{str(rad)}A'

def gen_plddt_colname(rad):
   return f'plddt_{str(rad)}A'

def gen_plddt_win_colname(win):
   return f'plddt_{str(int(win))}res'

# queries radii
for rad in RADIUS:
   PAEcolname = gen_pae_colname(rad)
   OUTPUT_COLS[PAEcolname] = []

   PLDDTcolname = gen_plddt_colname(rad)
   OUTPUT_COLS[PLDDTcolname] = []

# if we are going by plDDT window
if isinstance(PLDDT_WINDOW, list):
   
   for win in PLDDT_WINDOW:
      PLDDTcolname = gen_plddt_win_colname(win)
      OUTPUT_COLS[PLDDTcolname] = []

# keep track of the indices of problematic queries
PROBLEM_VARIANTS = []


# get information for each variant
for variant in VARIANTS_DF.itertuples():
   
   currentRow = int(getattr(variant, 'Index'))
   if currentRow % 1000 == 0:
      print(f'Currently processing row: {currentRow}\r', end='')
   
   # get corresponding model in memory
   MODEL: AlphafoldModel = None
   modelNotFound = False
   
   try:
      MODEL = AF_MODELS[getattr(variant, 'uniprot')]
   except KeyError:
      modelNotFound = True
   
   # check if residues are identical
   queryPos: int = None
   
   try:
      queryPos = int(getattr(variant, 'position'))
   except ValueError:
      pass
   
   WTaa = getattr(variant, 'WT')
   
   try:
      if modelNotFound == True:
         PROBLEM_VARIANTS.append(int(getattr(variant, 'Index')))
      elif WTaa != MODEL.get_residue(queryPos)[0] or queryPos == None:
         PROBLEM_VARIANTS.append(int(getattr(variant, 'Index')))
      else:
         for rad in RADIUS:
            if modelNotFound == False:
               PAEcolname = gen_pae_colname(rad)
               PAE_in_rad = np.round(MODEL.get_local_PAE(queryPos, rad, with_query_only=PAE_QUERY_ONLY),3)
               OUTPUT_COLS[PAEcolname].append(PAE_in_rad)

               PLDDTcolname = gen_plddt_colname(rad)
               PLDDT_in_rad = np.round(MODEL.get_local_plddt(queryPos, rad),3)
               OUTPUT_COLS[PLDDTcolname].append(PLDDT_in_rad)
            else:
               pass
         
         if isinstance(PLDDT_WINDOW, list):
            if modelNotFound == False:
               for win in PLDDT_WINDOW:
                  PLDDTcolname = gen_plddt_win_colname(win)
                  PLDDT_for_win = np.round(MODEL.get_plddt_window(queryPos, win)[0], 3)
                  OUTPUT_COLS[PLDDTcolname].append(PLDDT_for_win)
            else:
               pass
   
   except ValueError:
      PROBLEM_VARIANTS.append(int(getattr(variant, 'Index')))

VARIANTS_DF_GOOD = VARIANTS_DF.drop(PROBLEM_VARIANTS)
VARIANTS_DF_GOOD = VARIANTS_DF_GOOD.assign(**OUTPUT_COLS)

if SCRAP_COL == True:
   COLS_ADDED = list(OUTPUT_COLS.keys())
   COLS_TO_KEEP = ['uniprot', 'position', 'WT'] + COLS_ADDED
   VARIANTS_DF_GOOD = VARIANTS_DF_GOOD[COLS_TO_KEEP]

VARIANTS_DF_GOOD.to_csv(OUTPUT_PATH, sep='\t', index=False)

OUTPUT_DIR = os.path.dirname(OUTPUT_PATH)

if len(PROBLEM_VARIANTS) > 0:
   PROBLEM_DF: pd.DataFrame = VARIANTS_DF.iloc[PROBLEM_VARIANTS]
   OUTPUT_FILENAME = OUTPUT_PATH.split('/')[-1]
   PROBLEM_FILENAME = OUTPUT_FILENAME.replace('.csv','') + '_problems.csv'
   PROBLEM_PATH = os.path.join(OUTPUT_DIR, PROBLEM_FILENAME)
   PROBLEM_DF.to_csv(PROBLEM_PATH, sep='\t', index=False)

print(f'{len(PROBLEM_MODELS)} models could not be found.')
print(f'{len(PAE_NO_MATCH)} PAE matrices did not match/not found')
