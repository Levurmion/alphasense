from ModelPDB import ModelPDB
from ModelPAE import ModelPAE
from itertools import combinations
import numpy as np

class AlphafoldModel(ModelPDB,ModelPAE):
   
   def __init__(self, alphafoldPDB: str, alphafoldPAE: str):
      ModelPDB.__init__(self, alphafoldPDB)
      ModelPAE.__init__(self, alphafoldPAE)
   
   
   # evaluate the plddt score of a residue
   def evaluate_plddt(self, residue: int , threshold: float=70):
      
      if residue > self.length or residue <= 0:
         raise ValueError(f'Residue out of range. {self.name} ({self.uniprot}) is {self.length}aa long.')
      else:
         plDDT = self.get_plddt(residue)
         
         return (plDDT, plDDT >= threshold)


   # evaluate the plddt score as an average sliding window around a residue
   def evaluate_plddt_window(self, residue: int, window: int=5, threshold: float=70, return_val: bool=False):
      
      if window % 2 == 0 or window <= 0:
         raise ValueError('Window needs to be an integer positive odd number.')
      if residue > self.length or residue <= 0:
         raise ValueError(f'Residue out of range. {self.name} ({self.uniprot}) is {self.length}aa long.')
      else:
         windowLeft = residue-(window//2) if residue-(window//2) > 1 else 1
         windowRight = residue+(window//2)+1 if residue+(window//2) <= self.length else self.length + 1
         residuePlddts = [self.get_plddt(res) for res in range(windowLeft,windowRight,1)]
         averagePlddt = np.mean(residuePlddts)
         
         if return_val == False:
            return averagePlddt >= threshold
         elif return_val == True:
            return averagePlddt
   
   
   # getter method for all unique combinations of residue pairs in a list
   def get_PAEs(self, residues: list, with_query_only: bool=False):
      
      RESI_PAIRS: list = []
      
      if with_query_only == False:
         RESI_PAIRS = list(combinations(residues,2))
      elif with_query_only == True:
         query = residues[0]
         neighbours = residues
         del neighbours[0]
         for residue in neighbours:
            RESI_PAIRS.append((query, residue))
      
      RESI_PAIR_PAE = [] 
      
      for pair in RESI_PAIRS:
         resName1 = self.get_resID(pair[0])
         resName2 = self.get_resID(pair[1])
         pairName = f'{resName1}-{resName2}'
         # returns a list of tuples [(resName1, PAE1), (resName2, PAE2), ...]
         RESI_PAIR_PAE.append((pairName, self.get_pairwise_PAE(pair[0], pair[1])))
      
      return RESI_PAIR_PAE
   
   
   # getter method to get the average PAE of all residue pairs in a list
   def get_avg_PAE(self, residues: list, with_query_only: bool=False):
      
      pairwisePAEs = self.get_PAEs(residues, with_query_only=with_query_only)
      
      return np.mean(pairwisePAEs)
   
   
   
   
if __name__ == '__main__':
   alphafoldModel = AlphafoldModel('./test_files/AF-P0CG48-F1-model_v4.pdb','./test_files/AF-P0CG48-F1-predicted_aligned_error_v4.json')
   print(alphafoldModel.get_sequence())