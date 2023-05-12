import os
import sys
import pandas as pd

# establish files path
FILES_PATH = os.path.join(os.environ.get('PYTHONPATH'),'files')
OUTPUT_PATH = os.path.join(os.environ.get('PYTHONPATH'), 'files/outputs')

# python get_rsID.py <input file path> <output file path>
if len(sys.argv) != 3:
    print('Usage: python get_rsID.py <input file> <output file>')
    exit(1)

inputFilename = sys.argv[1]
outputFilename = sys.argv[2]

print('loading input file into dataframe...')
inputPath = os.path.join(FILES_PATH, inputFilename)
outputPath = os.path.join(OUTPUT_PATH, outputFilename)

input_df = pd.read_csv(inputPath, sep='\t', header=0)

rsID_column = list(input_df.loc[:,'Source DB ID'])

print('writing to outfile...')

lineCounter = 1
with open(outputPath, 'w') as outfile:

    for id in rsID_column:
        outfile.write(id + '\n')

        if lineCounter % 100000 == 0:
            print(f'writing line: {lineCounter}')
        
        lineCounter += 1