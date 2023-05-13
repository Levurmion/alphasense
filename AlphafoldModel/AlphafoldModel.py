from .ModelPDB import ModelPDB
from .ModelPAE import ModelPAE
from itertools import combinations
import numpy as np
import warnings

class AlphafoldModel(ModelPDB,ModelPAE):
   '''
   The `AlphafoldModel` class inherits from the `ModelPDB` and `ModelPAE` classes. Instantiate `AlphafoldModel` with filepaths to the Alphafold PDB structure and its corresponding JSON PAE matrix. The class exposes methods for querying a range of Alphafold model-specific properties such as per-residue plDDT and PAE scores.
   
   Methods:
   - `self.get_plddt()`: returns the plDDT score for a specific residue position.
   - `self.get_plddt_window()`: returns the plDDT score for a window around the queried residue.
   - `self.get_PAE()`: returns the PAE scores for all possible residue pairs in a list.
   - `self.get_avg_PAE()`: returns the average PAE score for all possible residue pairs in a list.
   - `self.get_local_plddt()`: returns the local plDDT score (average/list) within a radius of the queried residue.
   - `self.get_local_PAE()`: returns the local PAE score (average/list) of all possible residue pairs within a radius of the queried residue.
   '''
   
   
   def __init__(self, alphafoldPDB: str, alphafoldPAE: str):
      
      ModelPDB.__init__(self, alphafoldPDB)
      ModelPAE.__init__(self, alphafoldPAE)
      
      if len(self.chains) > 1:
         modelWarning = f'{alphafoldPDB} appears to be made of {len(self.chains)} chains. Are you sure it is an Alphafold model?'
         warnings.warn(modelWarning, category=Warning)
      
      if self.PAEshape[0] != self.chains[0].length:
         error = f'The length of the PDB structure and shape of the PAE matrix do not match. The PDB structure {self.chains[0].dbName} ({self.chains[0].dbAC}) is {self.chains[0].length} residues long and the PAE matrix has the dimensions {self.PAEshape}.'
         raise ValueError(error)
   
   
   def __str__(self):
      return f'{self.title}\n{self.chains[0].dbName} ({self.chains[0].dbSrc}: {self.chains[0].dbAC}) is a <class AlphafoldModel>. {self.chains[0].length} residues long.\n{self.chains[0].sequence}'
   
   # evaluate the plddt score of a residue
   def get_plddt(self, residue: int , threshold: float=70) -> tuple:
      
      if residue > self.chains[0].length or residue <= 0:
         raise ValueError(f'Residue out of range. {self.chains[0].dbName} ({self.chains[0].dbAC}) is {self.chains[0].length} residues long.')
      else:
         queryResidue = self.get_residue(residue, get_instance=True)
         plDDT = queryResidue.atoms[0].temp
         return (plDDT, plDDT >= threshold)


   # evaluate the plddt score as an average sliding window around a residue
   def get_plddt_window(self, residue: int, window: int=5, threshold: float=70) -> tuple:
      
      if window % 2 == 0 or window <= 0:
         raise ValueError('Window needs to be an integer positive odd number.')
      if residue > self.chains[0].length or residue <= 0:
         raise ValueError(f'Residue out of range. {self.chains[0].dbName} ({self.chains[0].dbAC}) is {self.chains[0].length}aa long.')
      else:
         windowLeft = residue-(window//2) if residue-(window//2) > 1 else 1
         windowRight = residue+(window//2)+1 if residue+(window//2) <= self.chains[0].length else self.chains[0].length + 1
         residuePlddts = [self.get_plddt(res) for res in range(int(windowLeft),int(windowRight),1)]
         averagePlddt = np.mean([resRecord[0] for resRecord in residuePlddts])
         
         return (averagePlddt, averagePlddt >= threshold)

   
   # getter method for all unique combinations of residue pairs in a list
   def get_PAE(self, residues: list, with_query_only: bool=False) -> tuple:
      
      if len(residues) <= 1:
         print('Cannot calculate PAE score for only 1 residue!')
         return -1
      
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
         resName1 = self.get_residue(pair[0])
         resName2 = self.get_residue(pair[1])
         pairName = f'{resName1}-{resName2}'
         # returns a list of tuples [(resName1, PAE1), (resName2, PAE2), ...]
         RESI_PAIR_PAE.append((self.get_pairwise_PAE(pair[0], pair[1]), pairName))
      
      return tuple(RESI_PAIR_PAE)
   
   
   # getter method to get the average PAE of all residue pairs in a list
   def get_avg_PAE(self, residues: list, with_query_only: bool=False) -> float:
      
      pairwisePAEs = self.get_PAE(residues, with_query_only=with_query_only)
      
      return np.around(np.mean([PAETuple[0] for PAETuple in pairwisePAEs]), 3) if pairwisePAEs != -1 else pairwisePAEs
   
   
   # get local PAE around a residue based on NN-search
   def get_local_PAE(self, residue: int, radius: float=5, average: bool=True, from_center: bool=False, with_query_only: bool=True):
      '''
      Calculate the PAE score between all residue pairs within a given radius of the query.
      
      - `residue: int` = residue number along the structure
      - `radius: float` = sphere radius (in Å) around the query within which to search for neighbouring residues
      - `average: bool` = Defaults to `True`, returns the calculated average of the PAEs between all residue pairs within the given `radius`. Toggle to `False` to return a list of tuples detailing every `(residue pair, PAE score)` discovered within the sphere.
      - `from_center: bool` = Defaults to `False`, searches for neighbours atom-to-atom from all atoms in a residue. Toggle to 'True' to search only from the coordinate center of the residue.
      - `with_query_only: bool` = Defaults to `True`, only considers residue pairs that involve the query. Toggle to `False` to to consider all possible residue pairs discovered within the sphere.
      '''
      
      residuesInBubble: list = self.get_residues_within(residue, radius, from_center=from_center, get_instance=True)
      residuePositions: list = [residue.position for residue in residuesInBubble]
      
      # reposition query residue to the start of the array
      residuePositions.remove(residue)
      residuePositions.insert(0, residue)
      
      if len(residuePositions) <= 1:
         print(f'No neighbouring residues found within {radius}Å of {self.get_residue(residue)}. `self.get_local_PAE()` returning -1.')
         return -1
      elif average == True:
         return self.get_avg_PAE(residuePositions, with_query_only=with_query_only)
      elif average == False:
         return self.get_PAE(residuePositions, with_query_only=with_query_only)
      
   
   # get local plddt around a residue based on NN-search
   def get_local_plddt(self, residue: int, radius: float=5, average: bool=True, from_center: bool=False):
      '''
      Calculate the average plDDT score for residues within a given radius of the query.
      
      - `residue: int` = residue number along the structure
      - `radius: float` = sphere radius (in Å) around the query within which to search for neighbouring residues
      - `average: bool` = Defaults to `True`, returns the calculated average of the plDDTs of all residues within the given `radius`. Toggle to `False` to return a list of tuples detailing every `(residue, plDDT score)` discovered within the sphere.
      - `from_center: bool` = Defaults to `False`, searches for neighbours atom-to-atom from all atoms in a residue. Toggle to 'True' to search only from the coordinate center of the residue.
      '''
   
      residuesInBubble = self.get_residues_within(residue, radius, from_center=from_center, get_instance=True)
      residuePositions = [residue.position for residue in residuesInBubble]
      
      if len(residuePositions) <= 1:
         print(f'No neighbouring residues found within {radius}Å of {self.get_residue(residue)}. Returning plDDT score for the query.')
         return self.get_plddt(residuePositions[0])[0]
      elif average == True:   # return float of the average
         return np.around(np.mean([self.get_plddt(residue)[0] for residue in residuePositions]),3)
      elif average == False:  # return nested tuples of ((residuePlddt, residue), ...)
         residueNames = [residue.resId for residue in residuesInBubble]
         residuePlddts = [self.get_plddt(residue)[0] for residue in residuePositions]
         return tuple(zip(residuePlddts, residueNames))
         
   
   
if __name__ == '__main__':

   from utils.utilities import function_timer
   
   alphafoldModel = AlphafoldModel('./test_files/AF-P0CG48-F1-model_v4.pdb','./test_files/AF-P0CG48-F1-predicted_aligned_error_v4.json')
   
   print(alphafoldModel.get_plddt_window(20,5))