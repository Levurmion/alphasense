import json
import numpy as np

class ModelPAE:

    def __init__(self, PAEMatrix: str):
        self.PAE = self.parse_PAE(PAEMatrix)


    def parse_PAE(self, PAEFilepath: str):

        with open(PAEFilepath, 'r') as PAE_JSON:
            PAE_PARSED = json.load(PAE_JSON)
            return np.array(PAE_PARSED[0]['predicted_aligned_error'])


    # getter method to obtain the average PAE score for a residue pair
    def get_pairwise_PAE(self, residue1: int, residue2: int):
        PAE_1 = self.PAE[residue1 - 1][residue2 - 1]
        PAE_2 = self.PAE[residue2 - 1][residue1 - 1]
        return (PAE_1 + PAE_2)/2


if __name__ == '__main__':
    PAE = ModelPAE('./test_files/AF-P04637-F1-predicted_aligned_error_v4.json')
    print(PAE.get_pairwise_PAE(40, 80))