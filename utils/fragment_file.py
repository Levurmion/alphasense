import os
import sys
import argparse
import subprocess
import re
import math
from utilities import readFile_as_generator

# establish files path
ROOT_DIR = os.path.join(os.environ.get('PYTHONPATH'))

parser = argparse.ArgumentParser(
    prog='fragment_file.py',
    description='breaks up a large file into a specified number of integer fragments',
)

parser.add_argument('-f', '--file', help='path to the file relative to the PYTHONPATH dir', required=True)
parser.add_argument('-o', '--output', help='path and name of the output file relative to the PYTHONPATH dir', required=True)
parser.add_argument('-n', '--num_fragments', help='target number of fragments to break the file into. If there are extra lines that cannot be included within the target file quantity, the overflow is going to be saved in an extra file.', required=True)
parser.add_argument('--headers', help='[option] to include the first line as headers, defaults to False, apply `--headers True` to include headers.', default=False)

try:
    args = parser.parse_args()
except:
    parser.print_help()
    sys.exit(1)

# get command line arguments
OUTPUT_FILENAME = args.output.split('/')[-1]
FILE_PATH = os.path.join(ROOT_DIR, args.file)
OUTFILE_PATH = os.path.dirname(os.path.join(ROOT_DIR, args.output))
FRAGMENTS = int(args.num_fragments)
HEADERS = bool(args.headers)

# use bash to quickly get the number of lines
# pipe the output of 'cat FILE_PATH' --> wc -l
# pipe the output of wc -l --> awk {print $1}
# awk {print $1} only prints the first column/field ($1) of the output from wc -l which is the number of lines
COMMAND = "wc -l " + FILE_PATH + " | awk '{print $1}'"
NUM_LINES = int(str(subprocess.check_output(COMMAND, shell=True)).strip("b'").strip('\\n'))

if HEADERS == True:
    NUM_LINES -= 1

# round off to the nearest integer
FRAGMENT_SIZE = int(NUM_LINES/FRAGMENTS)
CALC_FRAGMENTS = math.ceil(NUM_LINES/FRAGMENT_SIZE)
print(f'generating {CALC_FRAGMENTS} files from {FILE_PATH} with {NUM_LINES} lines.')
print(f'Every fragment is {FRAGMENT_SIZE} lines each.')

# get number of digits of NUM_LINES
COUNTER_DIGITS = len(str(NUM_LINES))

def format_digits(digits: int, number: int):

    numDigits = len(str(number))
    numLeadingZeroes = digits - numDigits
    return numLeadingZeroes*'0' + str(number)

outputFileExtension = '.' + re.search(r'\.(.*)$', OUTPUT_FILENAME).group(1)
outputFilePrefix = re.sub(r'\.(.*)$','',OUTPUT_FILENAME) + '_'

currentLineStart = 1
currentLineEnd = FRAGMENT_SIZE
currentFilename = outputFilePrefix + f'{format_digits(COUNTER_DIGITS,currentLineStart)}-{format_digits(COUNTER_DIGITS,currentLineEnd)}' + outputFileExtension
currentOutfilePath = os.path.join(OUTFILE_PATH, currentFilename)

FILE_HEADER = ''
currentFileLines = 0

# initialize script level variable to instantiate file object
outfile = None

for counter, line in enumerate(readFile_as_generator(FILE_PATH)):

    # if HEADERS are to be included in each fragment
    if counter == 0:
        if HEADERS == True:
            FILE_HEADER = line
        pass
    else:
        # if we need to open a new file
        if currentFileLines == 0:
            outfile = open(currentOutfilePath, 'w')
            print(f'writing to file: {currentFilename}')
            if HEADERS == True:
                outfile.write(FILE_HEADER)
        
        outfile.write(line)
        
        currentFileLines += 1

        # everytime we reach the FRAGMENT_SIZE, redo filenames and move on to next fragment in a separate file
        if currentFileLines == FRAGMENT_SIZE:
            outfile.close()
            currentFileLines = 0
            currentLineStart += FRAGMENT_SIZE
            currentLineEnd = (currentLineEnd + FRAGMENT_SIZE) if (currentLineEnd + FRAGMENT_SIZE < NUM_LINES) else NUM_LINES
            currentFilename = outputFilePrefix + f'{format_digits(COUNTER_DIGITS,currentLineStart)}-{format_digits(COUNTER_DIGITS,currentLineEnd)}' + outputFileExtension
            currentOutfilePath = os.path.join(OUTFILE_PATH, currentFilename)
    
    


