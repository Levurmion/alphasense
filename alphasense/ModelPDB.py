import math
import numpy as np
from LoadedKDTree import LoadedKDTree
from collections.abc import Iterable
from typing import Union

# create an Atom class to steamline access to residue attributes
class Atom:
    def __init__(self, residue: str, residuePos: int, atomName: str, temp: float, occupancy: float, chain: str, x: float, y: float, z: float):
        self.residuePos = residuePos
        self.atomName = atomName
        self.temp = temp
        self.resId = residue + str(residuePos)
        self.occupancy = occupancy
        self.chain = chain
        self.x = x
        self.y = y
        self.z = z
    
    # static method to calculate the euclidian distance between two Atom instances
    @staticmethod
    def euclid_dist(atom1, atom2) -> float:
        
        if not isinstance(atom1, Atom) or not isinstance(atom2, Atom):
            error = 'Static method Atom.euclid_dist(atom1, atom2) only takes instances of the Atom class'
            raise ValueError(error)
        else:
            return math.sqrt((atom1.x - atom2.x)**2 + (atom1.y - atom2.y)**2 + (atom1.z - atom2.z)**2)


# create a Residue class binding Atom-Residue relationships
class Residue:
    def __init__(self, residue: str, residuePos: int, chain: str):
        self.resId = residue + str(residuePos)
        self.aa = residue
        self.position = residuePos
        self.chain = chain
        self.atoms = []
        self.center = None 
    
    def central_coordinate(self):
        x, y, z = [
            np.around(np.mean([atom.x for atom in self.atoms]), decimals=3),
            np.around(np.mean([atom.y for atom in self.atoms]), decimals=3),
            np.around(np.mean([atom.z for atom in self.atoms]), decimals=3)
        ]
        self.center = (x, y, z)


# create a Chain class to handle multi-chain proteins
class Chain:

    def __init__(self, residues: list[Residue], chainId: str):
        self.residues = residues
        self.chainId = chainId
        self.length = len(residues)
        self.sequence = ''.join([residue.aa for residue in residues])
    
    # if users just want to get the sequence
    def __str__(self):
        printedString = f'chain {self.chainId} ({self.length} residues):\n' + self.sequence
        return printedString
    
    # iterators to loop through the Residue class instances
    def __iter__(self):
        self._currentPos = 0
        return self

    def __next__(self):
        if self._currentPos < self.length:
            currentResidue = self.residues[self._currentPos]
            self._currentPos += 1
            return currentResidue
        elif self._currentPos >= self.length:
            raise StopIteration


# PDB object class
class ModelPDB():

    def __init__(self, model: str):
        self.uniprot, self.name, self.length, self.sequence, self.atoms, self.residues = self.__parse_model(model)
        self.atomicKdtree = LoadedKDTree(self.atoms, self.__retrieve_atomic_coord)


    # dunder methods
    def __str__(self):
        return f'Model {self.name} ({self.uniprot}) is {self.length} residues long.'


    # private methods (mangled)
    def __readFile_as_generator(self, filePath: str) -> str:
        for line in open(filePath):
            yield line


    def __parse_atom_line(self, pdbLine: str) -> list:
        '''
        Note - subtract 1 from these column residuePoss to get index:
        COLUMNS        DATA  TYPE    FIELD        DEFINITION
        -------------------------------------------------------------------------------------
         1 -  4        Record name   "ATOM  "
         7 - 11        Integer       serial       Atom  serial number.
        13 - 16        Atom          name         Atom name.
        17             Character     altLoc       Alternate location indicator.
        18 - 20        Residue name  resName      Residue name.
        22             Character     chainID      Chain identifier.
        23 - 26        Integer       resSeq       Residue sequence number.
        27             AChar         iCode        Code for insertion of residues.
        31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
        39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
        47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
        55 - 60        Real(6.2)     occupancy    Occupancy.
        61 - 66        Real(6.2)     tempFactor   Temperature  factor. <-- B-factor/plDDT
        77 - 78        LString(2)    element      Element symbol, right-justified.
        79 - 80        LString(2)    charge       Charge  on the atom.
        '''

        atomSer: int = int(pdbLine[6:12].strip())
        atomName: str = pdbLine[12:16].strip()
        altLoc: str = pdbLine[16].strip()
        resName: str = pdbLine[17:20].strip()
        chain: str = pdbLine[21]
        resPos: int = int(pdbLine[22:26].strip())
        x: float = float(pdbLine[30:38].strip())
        y: float = float(pdbLine[38:46].strip())
        z: float = float(pdbLine[46:54].strip())
        occupancy: float = float(pdbLine[54:60].strip())
        temp: float = float(pdbLine[60:66].strip())

        return (atomSer, atomName, altLoc, resName, chain, resPos, x, y, z, occupancy, temp)


    def __parse_seqres_line(self, pdbLine: str) -> list:
        residueFields: list = pdbLine[19:70].split()
        return residueFields


    def __parse_model(self, modelFilepath: str) -> list:

        AA_DICT = {'VAL': 'V', 'ILE': 'I', 'LEU': 'L', 'GLU': 'E', 'GLN': 'Q',
                       'ASP': 'D', 'ASN': 'N', 'HIS': 'H', 'TRP': 'W', 'PHE': 'F', 'TYR': 'Y',
                       'ARG': 'R', 'LYS': 'K', 'SER': 'S', 'THR': 'T', 'MET': 'M', 'ALA': 'A',
                       'GLY': 'G', 'PRO': 'P', 'CYS': 'C'}

        uniprot: str = ''
        name: str = ''
        length: int = 0
        sequence: str = ''
        atoms: list = []
        residues: list = []
        
        currentPos = 0
        currentResidue = None

        for line in self.__readFile_as_generator(modelFilepath):

            # split by whitespace characters
            LINE_DTYPE = line[0:6].strip()

            if LINE_DTYPE == 'DBREF':
                uniprot = line[33:41].strip()
                name = line[42:54].strip()

            if LINE_DTYPE == 'SEQRES':
                for res in self.__parse_seqres_line(line):
                    sequence += AA_DICT[res]

            elif LINE_DTYPE == 'ATOM':
                atomSer, atomName, altLoc, resName, chain, resPos, x, y, z, occupancy, temp = self.__parse_atom_line(line)
                
                atomInstance = Atom(AA_DICT[resName],resPos,atomName,temp,occupancy,chain,x,y,z)
                atoms.append(atomInstance)
                
                if resPos > currentPos:
                    if currentResidue != None:
                        currentResidue.central_coordinate()
                        residues.append(currentResidue)
                        currentResidue = Residue(AA_DICT[resName], resPos, chain)
                        currentResidue.atoms.append(atomInstance)
                    else:
                        currentResidue = Residue(AA_DICT[resName], resPos, chain)
                        currentResidue.atoms.append(atomInstance)
                elif resPos == currentPos:
                    currentResidue.atoms.append(atomInstance)
                    
                currentPos = resPos
            
            elif LINE_DTYPE == 'TER':
                currentResidue.central_coordinate()
                residues.append(currentResidue)
                    
        length = len(sequence)

        return [uniprot, name, length, sequence, atoms, residues]


    def __retrieve_atomic_coord(self, atom) -> tuple:
        return (atom.x, atom.y, atom.z)


    # public methods
    def get_residues(self, *residueArgs, get_instance: bool=False) -> list:

        residues = residueArgs[0] if len(residueArgs) == 1 and type(
            residueArgs[0]) == list else residueArgs
        residueList: list = []
        
        firstRes = min(residues)
        lastRes = max(residues)
        
        if firstRes < 1 or lastRes > self.length:
            error = f'Residues provided out of range. {self.name} ({self.uniprot}) is {self.length} residues long.'
            raise ValueError(error)

        if get_instance == False:
            for residue in residues:
                residueList.append(self.residues[residue-1].resId)
        elif get_instance == True:
            for residue in residues:
                residueList.append(self.residues[residue-1])

        return residueList


    def get_residues_within(self, query: Union[type[Atom], Iterable, int], radius: float=5.0, from_center: bool=False, get_instance: bool=False) -> list:
        '''
        `self.get_residues_within(query, radius: float)` finds all residues within a given `radius` (in Amstrongs) of the `query`. The `query` parameter accepts several types as arguments:
        
        - `coordinates: Iterable[float]` = A 3D point in space from which to perform the search.
        - `atom: <class Atom>` = An instance of the `Atom` class.
        - `residue: int` = The residue number along the primary sequence of the protein. By default, this will search a space from all atoms registered with the residue within the given `radius`. Optionally, toggle the `from_center` parameter to `True` to only search from the spatial center of the residue.
        '''
        
        stdQuery = None
        
        if isinstance(query, Atom):
            stdQuery = [self.___retrieve_atomic_coord(query)]
        elif isinstance(query, Iterable) and len(query) == 3 and (isinstance(entry, float) for entry in query):
            stdQuery = [query]
        elif isinstance(query, int) and 0 < query <= self.length:
            if from_center == False:
                stdQuery = [atom for atom in self.get_residues(query, get_instance=True)[0].atoms]
            elif from_center == True:
                stdQuery = [self.get_residues(query, get_instance=True)[0].center]
                    
        if stdQuery == None:
            error = f'The `query` parameter only takes the following as arguments: <class Atom> instance, an iterable representing a 3D point in space, or a residue number. {self.name} ({self.uniprot}) is {self.length} residues long.'
            raise ValueError(error)
        
        neighbours = set()
        
        for point in stdQuery:
            neighbourAtoms = self.atomicKdtree.query_within_radius(point, radius)
            neighbourResidues = [atom.residuePos for atom in neighbourAtoms]
            for neighbour in neighbourResidues:
                neighbours.add(neighbour)
                
        return self.get_residues(list(neighbours)) if get_instance == False else self.get_residues(list(neighbours), get_instance=get_instance)


if __name__ == '__main__':
    alphafold = ModelPDB('./test_files/AF-P04637-F1-model_v4.pdb')
    
    chain = Chain(alphafold.residues, 'A')
    
    print(chain)
