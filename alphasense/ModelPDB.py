
class ModelPDB:
    
    def __init__(self, model: str):
        self.dbref, self.length, self.sequence, self.plddt, self.coordinates = self.parse_model(model)


    # create a Coordinate inner class to steamline access to residue x,y,z CA coordinates
    class Coordinate:

        def __init__(self, x: float, y: float, z: float):
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

        dbref: dict = {}
        length: int = 0
        sequence: list = []
        plddt: list = []
        coordinates: dict = {}

        currentResNum = 0
        
        for line in self.readFile_as_generator(modelFilepath):

            # split by whitespace characters
            LINE_FIELDS = line.split()
            LINE_DTYPE = LINE_FIELDS[0]
            
            if LINE_DTYPE == 'DBREF':
                dbref['uniprot'] = LINE_FIELDS[6]
                dbref['name'] = LINE_FIELDS[7]

            elif LINE_DTYPE == 'SEQRES':
                RES_FIELDS = LINE_FIELDS[4::]
                for residue in RES_FIELDS:
                    sequence.append(AA_DICT_LTS[residue.upper()])
            
            elif LINE_DTYPE == 'ATOM':
                resNum = int(LINE_FIELDS[5])
                atomType = LINE_FIELDS[2]
                if resNum > currentResNum:
                    plddt.append(float(LINE_FIELDS[10]))
                    currentResNum += 1
                if atomType == 'CA':
                    coorX = LINE_FIELDS[6]
                    coorY = LINE_FIELDS[7]
                    coorZ = LINE_FIELDS[8]
                    coordinates[str(resNum)] = self.Coordinate(coorX, coorY, coorZ)
                else:
                    pass

        length = len(sequence)

        if len(plddt) == length:
            pass
        else:
            raise RuntimeError('Cannot Intialize Alphafold object: Parsed sequence length and number of extracted plDDT scores did not match.')

        return [dbref, length, sequence, plddt, coordinates]


    # getter method to streamline shift in indexing
    def get_residue(self, *args: list, coordinates=False, plddt=False):

        resList = args if type(args[0]) == int else args[0]

        residues: list = []

        if coordinates == True and plddt == True:
            residues = residues = [(self.sequence[resNum - 1], self.coordinates[str(resNum)], self.plddt[resNum - 1]) for resNum in resList]
        elif coordinates == True:
            residues = [(resNum, self.sequence[resNum - 1], self.coordinates[str(resNum)]) for resNum in resList]
        elif plddt == True:
            residues = [(resNum, self.sequence[resNum - 1], self.plddt[resNum - 1]) for resNum in resList]
        else:
            residues = [(resNum, self.sequence[resNum - 1]) for resNum in resList]

        return residues
    

    # getter method to get plDDT for a residue
    def get_plddt(self, residue):
        return self.plddt[residue - 1]
    
    # getter method to get residue Ca coordinate
    # returns a residue Coordinate instance
    def get_Ca_coord(self, residue):
        return self.coordinates[str(residue)]



if __name__ == '__main__':
    alphafold = ModelPDB('./test_files/AF-P04637-F1-model_v4.pdb')
    print(alphafold.get_residue([1,2],plddt=True))
    print(alphafold.get_plddt(2))
    print(alphafold.get_Ca_coord(2).x)