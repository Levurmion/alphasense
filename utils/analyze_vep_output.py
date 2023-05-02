import os
import sys
import pandas as pd
import re
from utilities import function_timer, AA_DICT_LTS

print('starting analysis...')

# establish files path
ROOT_DIR = os.path.join(os.environ.get('PYTHONPATH'))

# get CLI arguments
if len(sys.argv) != 3:
    print('Usage: python analyze_vep_output.py <vep output filepath> <reference filepath>')
    exit(1)

VEP_OUTPUT_FILEPATH = os.path.join(ROOT_DIR, sys.argv[1])
REF_FILEPATH = os.path.join(ROOT_DIR, sys.argv[2])

print('loading vep output into dataframe')
vep_output_df = pd.read_csv(VEP_OUTPUT_FILEPATH, sep='\t', header=0)
print('loading reference file into dataframe')
ref_df = None

try:
    print('reading reference dataframe from cache')
    ref_df = pd.read_pickle('uniprot_human_variants_trimmed.pkl')
except FileNotFoundError:
    ref_df = pd.read_csv(REF_FILEPATH, sep='\t', header=0)
    print('setting index on reference file...')

    @function_timer('index set in: ')
    def set_index(df):
        return df.set_index(['Chromosome Coordinate','AC'])

    ref_df = set_index(ref_df)

    print('caching reference dataframe')
    # cache dataframe
    ref_df.to_pickle('uniprot_human_variants_trimmed.pkl')


def parse_ref_mutation(mutation):
    parsed_mutation = mutation.replace('p.','')
    parsed_mutation = re.sub(r'(?<=\D)\d+(?=\D)','/',parsed_mutation)
    parsed_mutation = parsed_mutation.split('/')
    try:
        parsed_mutation = [AA_DICT_LTS[aa.upper()] for aa in parsed_mutation]
        return parsed_mutation
    except KeyError:
        return False

# uniprot not found in reference
vep_no_uniprot = []
vep_uniprot_no_match =[]

vep_matching_mutation = []
vep_mutation_no_match = []

lineCounter = 0
for idx, row in vep_output_df.iterrows():

    # vep output fields
    HGVS = row['Uploaded_variation']
    uniprot = re.sub(r'\.\s*([^\s.]+)','',row['SWISSPROT'])
    mutation = row['Amino_acids'].split('/')    # [ref,mutant]

    # reference file records
    if uniprot != '-':
        ref_records_HGVS = ref_df.loc[HGVS]

        try:
            ref_records_uniprot = ref_records_HGVS.loc[uniprot]
            
            ref_records_mutation = ref_records_HGVS.loc[uniprot,'Variant AA Change'].iloc[0]
            recorded_mutation = parse_ref_mutation(ref_records_mutation)
            
            try:
                if recorded_mutation == False:
                    vep_mutation_no_match.append(row)
                elif recorded_mutation[0] == mutation[0] and recorded_mutation[1] == mutation[1]:
                    vep_matching_mutation.append(row)
                else:
                    vep_mutation_no_match.append(row)
            except IndexError:
                vep_mutation_no_match.append(row)
        
        # only one record
        except AttributeError:
            ref_records_mutation = ref_records_HGVS.loc[uniprot,'Variant AA Change']
            recorded_mutation = parse_ref_mutation(ref_records_mutation)

            try:
                if recorded_mutation == False:
                    vep_mutation_no_match.append(row)
                elif recorded_mutation[0] == mutation[0] and recorded_mutation[1] == mutation[1]:
                    vep_matching_mutation.append(row)
                else:
                    vep_mutation_no_match.append(row)
            except IndexError:
                vep_mutation_no_match.append(row)

        except KeyError:
            vep_uniprot_no_match.append(row)

    else:
        vep_no_uniprot.append(row)
    
    lineCounter += 1

    if lineCounter % 1000 == 0:
        print('processing line: ',lineCounter, '\r', end='')


with open('analysis_summary.txt','w') as summary:

    summary.write('vep records with no uniprot: ' + str(len(vep_no_uniprot)) + '\n')
    summary.write('vep records without matching uniprot: ' + str(len(vep_uniprot_no_match)) + '\n')
    summary.write('vep records with matching mutation: ' + str(len(vep_matching_mutation)) + '\n')
    summary.write('vep records without matching mutation: ' + str(len(vep_mutation_no_match)) + '\n')

    summary.write('\n')

    summary.write('### RECORDS WITHOUT MATCHING UNIPROT' + '\n')
    for record in vep_uniprot_no_match:
        summary.write(str(record) + '\n')
    
    summary.write('\n')
    summary.write('----------------------------------\n')
    summary.write('### RECORDS WITHOUT MATCHING MUTATION (BUT UNIPROT MATCHES)' + '\n')
    for record in vep_mutation_no_match:
        summary.write(str(record) + '\n')
