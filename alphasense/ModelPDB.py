
class ModelPDB:
    
    def __init__(self, model: str):
        self.uniprot, self.name, self.length, self.sequence = self.parse_model(model)


    # create a Residue inner class to steamline access to residue attributes
    class Residue:

        def __init__(self, residue: str, position: int, plddt: float, x: float, y: float, z: float):
            self.residue = residue
            self.position = position
            self.plddt = plddt
            self.resID = residue + str(position)
            self.x = x
            self.y = y
            self.z = z


    def readFile_as_generator(self, filePath: str):
        for line in open(filePath):
            yield line


    def parse_model(self, modelFilepath: str):

        AA_DICT_LTS = {'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q',
        'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',
        'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',
        'GLY':'G', 'PRO':'P', 'CYS':'C'}

        uniprot: str = ''
        name: str = ''
        length: int = 0
        sequence: list = []

        currentResNum: int = 0
        currentResName: str = ''
        currentPlddt: float = 0.0
        currentX: float = 0.0
        currentY: float = 0.0
        currentZ: float = 0.0
        
        resCaFound = False
        
        for line in self.readFile_as_generator(modelFilepath):

            # split by whitespace characters
            LINE_FIELDS = line.split()
            LINE_DTYPE = LINE_FIELDS[0]
            
            if LINE_DTYPE == 'DBREF':
                uniprot = LINE_FIELDS[6]
                name = LINE_FIELDS[7]

            elif LINE_DTYPE == 'ATOM':
                resNum = int(LINE_FIELDS[5])
                atomType = LINE_FIELDS[2]
                if resNum > currentResNum:
                    currentPlddt = float(LINE_FIELDS[10])
                    residue = LINE_FIELDS[3]
                    currentResName = AA_DICT_LTS[residue.upper()]
                    currentResNum += 1
                if atomType == 'CA':
                    currentX = float(LINE_FIELDS[6])
                    currentY = float(LINE_FIELDS[7])
                    currentZ = float(LINE_FIELDS[8])
                    resCaFound = True
                else:
                    pass
            
            if resCaFound == True:
                sequence.append(self.Residue(
                    currentResName,
                    currentResNum,
                    currentPlddt,
                    currentX,
                    currentY,
                    currentZ)
                )
                resCaFound = False

        length = len(sequence)

        return [uniprot, name, length, sequence]


    # shorthand getter method to get residue properties
    def get_residues(self, *residueArgs, resId: bool=False, coordinates=False, plddt=False):

        residues = residueArgs[0] if len(residueArgs) == 1 and type(residueArgs[0]) == list else residueArgs
        residueProps: list = []
        
        if resId == False:
            for residue in residues:
                residueProps.append(self.sequence[residue-1].residue)
        elif resId == True:
            for residue in residues:
                residueProps.append(self.get_resID(residue))
        
        coordList: list = []
        if coordinates == True:
            for residue in residues:
                coordList.append(self.get_Ca_coord(residue))
        
        plddtList: list = []
        if plddt == True:
            for residue in residues:
                plddtList.append(self.get_plddt(residue))
            
        if coordinates == True and plddt == True:
            residueProps = list(zip(residueProps, coordList, plddtList))
        elif coordinates == True:
            residueProps = list(zip(residueProps, coordList))
        elif plddt == True:
            residueProps = list(zip(residueProps, plddtList))

        return residueProps
    

    # getter method to get plDDT for a residue
    def get_plddt(self, residue: int):
        return self.sequence[residue - 1].plddt
    
    # getter method to get residue Ca coordinate
    # returns a residue Coordinate instance
    def get_Ca_coord(self, residue: int):
        res = self.sequence[residue - 1]
        return (res.x, res.y, res.z)
    
    
    # getter method to get residue ID in the format resName[position]
    def get_resID(self, residue: int):
        return self.sequence[residue - 1].resID



if __name__ == '__main__':
    alphafold = ModelPDB('./test_files/AF-P04637-F1-model_v4.pdb')
    print(alphafold.get_residues([1,2,8,9],resId=True,plddt=True,coordinates=True))