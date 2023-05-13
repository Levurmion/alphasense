import math
import numpy as np
from .LoadedKDTree import LoadedKDTree
from collections.abc import Iterable, Collection
from typing import Union
from collections import defaultdict

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
            error = 'Static method Atom.euclid_dist(atom1, atom2) only takes instances of the <Atom> class.'
            raise ValueError(error)
        else:
            return math.sqrt((atom1.x - atom2.x)**2 + (atom1.y - atom2.y)**2 + (atom1.z - atom2.z)**2)


# create a Residue class binding Atom-Residue relationships
class Residue:
    def __init__(self, residue: str, residuePos: int, refPosition: int, chain: str):
        self.resId = residue + str(residuePos)
        self.refId = residue + str(refPosition)
        self.aa = residue
        self.position = residuePos
        self.refPosition = refPosition
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

    def __init__(self, residues: list, chainId: str, dbSrc: str, dbAC: str, dbName: str, refStart: int):
        self.residues = residues
        self.chainId = chainId
        self.refStart = refStart
        self.dbSrc = dbSrc
        self.dbAC = dbAC
        self.dbName = dbName
        self.sequence = ''.join([residue.aa for residue in residues])
        self.length = len(self.sequence)
    
    # if users just want to get the sequence
    def __str__(self):
        return self.sequence
    
    def get_info(self):
        infoText = f'chain {self.chainId}: {self.dbName}\n{self.dbSrc}: {self.dbAC}\n{self.length} residues\nStructure starts with residue {self.refStart} with respect to the reference sequence.'
        print(infoText)
        return {
            'id': self.chainId,
            'dbSrc': self.dbSrc,
            'dbAC': self.dbAC,
            'dbIdCode': self.dbName,
            'length': self.length,
            'refStart': self.refStart,
            'sequence': self.sequence
        }
    
    # iterators to loop through the Residue class instances
    def __iter__(self):
        self._currentStructPos = 0
        return self

    def __next__(self):
        if self._currentStructPos < self.length:
            currentResidue = self.residues[self._currentStructPos]
            self._currentStructPos += 1
            return currentResidue
        elif self._currentStructPos >= self.length:
            raise StopIteration


##############################################################################################################
#################################### !!!!!!!!! PDB OBJECT CLASS !!!!!!!!! ####################################
##############################################################################################################

class ModelPDB():

    def __init__(self, model: str):
        self.pdbId = None
        self.title, self.atoms, self.chains = self.__parse_model(model)
        self.atomicKdtree = LoadedKDTree(self.atoms, self.__retrieve_atomic_coord)


    # dunder methods
    def __str__(self):
        if len(self.chains) > 1:
            totalResidues = sum([chain.length for chain in self.chains])
            infoString = f'{self.title}\nModel {self.pdbId}, {len(self.chains)} chains.\nTotal length: {totalResidues} residues'
            return infoString
        else:
            return f'Model {self.pdbId} is a single chain protein.\n{self.chains[0].length} residues'


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
        refPos: int = int(pdbLine[22:26].strip())
        x: float = float(pdbLine[30:38].strip())
        y: float = float(pdbLine[38:46].strip())
        z: float = float(pdbLine[46:54].strip())
        occupancy: float = float(pdbLine[54:60].strip())
        temp: float = float(pdbLine[60:66].strip())

        return (atomSer, atomName, altLoc, resName, chain, refPos, x, y, z, occupancy, temp)


    def __parse_dbref_line(self, pdbLine: str) -> list:
        '''
        Note - subtract 1 from these column residuePoss to get index:
        COLUMNS       DATA TYPE     FIELD              DEFINITION
        -----------------------------------------------------------------------------------
        1 -  6       Record name   "DBREF "
        8 - 11       IDcode        idCode             ID code of this entry.
        13            Character     chainID            Chain  identifier.
        15 - 18       Integer       seqBegin           Initial sequence number of the PDB sequence segment.
        19            AChar         insertBegin        Initial  insertion code of the PDB  sequence segment.
        21 - 24       Integer       seqEnd             Ending sequence number of the PDB  sequence segment.
        25            AChar         insertEnd          Ending insertion code of the PDB  sequence segment.
        27 - 32       LString       database           Sequence database name.
        34 - 41       LString       dbAccession        Sequence database accession code.
        43 - 54       LString       dbIdCode           Sequence  database identification code.
        56 - 60       Integer       dbseqBegin         Initial sequence number of the database seqment.
        61            AChar         idbnsBeg           Insertion code of initial residue of the segment, if PDB is the reference.
        63 - 67       Integer       dbseqEnd           Ending sequence number of the database segment.
        68            AChar         dbinsEnd           Insertion code of the ending residue of the segment, if PDB is the reference.
        '''
        pdbId = pdbLine[7:11].strip()
        chainId = pdbLine[12].strip()
        dbSrc = pdbLine[26:32].strip()
        dbAC = pdbLine[33:41].strip()
        dbName = pdbLine[42:54].strip()
        return [pdbId, chainId, dbSrc, dbAC, dbName]


    def __parse_model(self, modelFilepath: str) -> list:

        AA_DICT = {'VAL': 'V', 'ILE': 'I', 'LEU': 'L', 'GLU': 'E', 'GLN': 'Q',
                       'ASP': 'D', 'ASN': 'N', 'HIS': 'H', 'TRP': 'W', 'PHE': 'F', 'TYR': 'Y',
                       'ARG': 'R', 'LYS': 'K', 'SER': 'S', 'THR': 'T', 'MET': 'M', 'ALA': 'A',
                       'GLY': 'G', 'PRO': 'P', 'CYS': 'C'}

        title: str = ''
        atoms: list = []
        chains: list = []
        
        # dictionary for information about different Chains
        chainIdentifiers: dict = {}
        
        # tracker variables for parsing ATOM line
        currentStructPos = 0    # residue number on the actual given structure
        currentRefPos = 0       # residue number on the reference sequence
        currentResidue = None
        currentChain = None
        currentResidues: list = []
        resStart: int = 1

        for line in self.__readFile_as_generator(modelFilepath):

            # split by whitespace characters
            LINE_DTYPE = line[0:6].strip()
            
            if LINE_DTYPE == 'TITLE':
                titleText = line[10::].strip('\n')
                title += titleText

            if LINE_DTYPE == 'DBREF':
                pdbId, chainId, dbSrc, dbAC, dbName = self.__parse_dbref_line(line)
                
                if self.pdbId == None:
                    self.pdbId = pdbId
                
                chainIdentifiers[chainId] = {
                    'dbSrc': dbSrc,
                    'dbAC': dbAC,
                    'dbName': dbName
                }

            elif LINE_DTYPE == 'ATOM':
                atomSer, atomName, altLoc, resName, chain, refPos, x, y, z, occupancy, temp = self.__parse_atom_line(line)

                if currentChain != None:
                        
                    if refPos > currentRefPos:
                        
                        currentStructPos += 1
                        
                        atomInstance = Atom(AA_DICT[resName],currentStructPos,atomName,temp,occupancy,chain,x,y,z)
                        atoms.append(atomInstance)
                        
                        if currentResidue != None:
                            currentResidue.central_coordinate()
                            currentResidues.append(currentResidue)
                            currentResidue = Residue(AA_DICT[resName], currentStructPos, refPos, chain)
                            currentResidue.atoms.append(atomInstance)
                            
                        else:
                            currentResidue = Residue(AA_DICT[resName], currentStructPos, refPos, chain)
                            currentResidue.atoms.append(atomInstance)
                            
                        currentRefPos = refPos
                        
                    elif refPos == currentRefPos:
                        atomInstance = Atom(AA_DICT[resName],currentStructPos,atomName,temp,occupancy,chain,x,y,z)
                        atoms.append(atomInstance)
                        currentResidue.atoms.append(atomInstance)
                            
                else:
                    resStart = refPos
                    currentStructPos += 1
                    currentRefPos = refPos
                    currentResidue = Residue(AA_DICT[resName], currentStructPos, refPos, chain)
                    atomInstance = Atom(AA_DICT[resName],currentStructPos,atomName,temp,occupancy,chain,x,y,z)
                    atoms.append(atomInstance)
                    currentResidue.atoms.append(atomInstance)
                    currentChain = chain
            
            elif LINE_DTYPE == 'TER':
                currentResidue.central_coordinate()
                currentResidues.append(currentResidue)
                chainInfo = chainIdentifiers[currentChain]
                chains.append(Chain(currentResidues, currentChain, chainInfo['dbSrc'], chainInfo['dbAC'], chainInfo['dbName'], resStart))
                currentResidues = []
                currentResidue = None
                currentStructPos = 0
                currentRefPos = 0
                currentChain = None
                resStart = 1

        return [title, atoms, chains]


    def __retrieve_atomic_coord(self, atom) -> tuple:
        return (atom.x, atom.y, atom.z)
    

    # public methods
    
    def get_chain(self, queryChain: str) -> Union[Chain, bool]:
        
        for chain in self.chains:
            if chain.chainId == queryChain:
                return chain
            else:
                pass
        
        print(f'Chain {queryChain} does not exist. Returning False.')
        # if chain not found
        return False
    
    
    def get_residue(self, residue: Union[int, tuple], get_instance: bool=False, by_ref: bool=False) -> Union[str, Residue]:
        
        queryChain = None
        queryResidue = None

        # if multi-chain protein, enforce inputs as tuples
        if len(self.chains) > 1:

            if (isinstance(residue, tuple) == False and isinstance(residue, list) == False) or len(residue) != 2:
                inputError = f'{self.pdbId} is a {len(self.chains)}-chain protein. Please specify chain and residue position as a list/tuple: (chain: str, residue: int).'
                raise ValueError(inputError)
            
            else:
                queryChain, queryResidue = [residue[0], residue[1] - 1]

        # if single-chain protein, allow inputs as integer residue positions
        elif len(self.chains) == 1:
            
            if (isinstance(residue, tuple) or isinstance(residue, list)) and len(residue) == 2:
                queryChain, queryResidue = [residue[0], residue[1] - 1]
            
            elif isinstance(residue, int):
                queryResidue = residue - 1
            
            else:
                inputError = f'{self.pdbId} is a single-chain protein. Query directly with a residue position or the correct chain and residue position as a list/tuple: (chain: str, residue: int)'
                raise ValueError(inputError)
        
        # chain was not found - raise error, queryChain == None if query is an integer residue position for single-chain proteins
        proteinChain = self.get_chain(queryChain) if queryChain != None else self.chains[0]
        
        if proteinChain == False:
            chainNotFound = f'Chain {queryChain} does not exist.'
            raise ValueError(chainNotFound)
        
        # queried residue out of range - raise error
        chainLength = proteinChain.length
        
        if queryResidue < 0 or queryResidue > chainLength - 1:
            print(residue)
            residueOutOfRange = f'Residue out of range. Chain {queryChain} is {chainLength} residues long.'
            raise ValueError(residueOutOfRange)
        
        # if all checks passed, get the residue in the chain
        if get_instance == True:
            return proteinChain.residues[queryResidue]
        elif get_instance == False:
            return proteinChain.residues[queryResidue].resId if by_ref == False else proteinChain.residues[queryResidue].refId
        

    def get_residues_within(self, query: Union[Atom, tuple, int], radius: float=5.0, from_center: bool=False, get_instance: bool=False, by_ref: bool=False) -> Union[list, dict]:
        '''
        `self.get_residues_within(query, radius: float)` finds all residues within a given `radius` (in Amstrongs) of the `query`. The `query` parameter accepts several types as arguments:
        
        - `atom: <class Atom>` = An instance of the `Atom` class.
        - `coordinates: tuple[float, float, float]` = A 3D point in space from which to perform the search.
        - `residue: tuple[str, float] | int` = The residue number along the primary structure of the protein within a specified chain. By default, this will search a space from all atoms registered with the residue within the given `radius`. Optionally, toggle the `from_center` parameter to `True` to only search from the spatial center of the residue.
        - `from_center: bool` = Defaults to `False`, searches for neighbours atom-to-atom from all atoms in a residue. Toggle to 'True' to search only from the coordinate center of the residue.
        '''
        
        stdQuery = None
        
        if isinstance(query, Atom):
            stdQuery = [self.__retrieve_atomic_coord(query)]
        elif isinstance(query, Iterable) and len(query) == 3 and (isinstance(entry, float) for entry in query):
            stdQuery = [query]
        # neither <class Atom> instance nor coordinate tuple, treat as direct residue query
        else:
            if from_center == False:
                stdQuery = [atom for atom in self.get_residue(query, get_instance=True).atoms]
            elif from_center == True:
                stdQuery = [self.get_residue(query, get_instance=True).center]
                    
        if stdQuery == None:
            inputError = f'The `query` parameter only takes the following as arguments: <class Atom> instance, an iterable representing a 3D point in space, or a residue in a specific chain as a list/tuple: (chain: str, residue: int) or an int residue position for single-chain proteins.'
            raise ValueError(inputError)
        
        uniqueNeighRes = set()
        
        # point is a <class Atom> instance
        for point in stdQuery:
            neighbourAtoms = self.atomicKdtree.query_within_radius(point, radius)
            neighbourResidues = [(atom.chain, atom.residuePos) for atom in neighbourAtoms]
            for neighbour in neighbourResidues:
                uniqueNeighRes.add(neighbour)
                
        neighboursPerChain = defaultdict(list)
        
        for res in uniqueNeighRes:
            neighboursPerChain[res[0]].append(self.get_residue(res,get_instance=get_instance,by_ref=by_ref))
            
        # sort residues returned from each chain
        neighbourChains = neighboursPerChain.keys()
        
        for chain in neighbourChains:
            neighboursPerChain[chain].sort(key= lambda res: int(res[1:]) if get_instance == False else res.position)
        
        # if single-chain, don't bother returning the entire dictionary
        if len(self.chains) == 1:
            key = list(neighboursPerChain.keys())[0]
            return neighboursPerChain[key]
        # if multi-chain, return as a dictionary of {'chain': [neighbours]}
        else:
            return dict(neighboursPerChain)
        


if __name__ == '__main__':
    pdb = ModelPDB('./test_files/AF-P04637-F1-model_v4.pdb')

    pdb.get_chain('A').get_info()
