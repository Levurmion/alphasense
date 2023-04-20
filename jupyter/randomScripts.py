import os

# establish files path
FILES_PATH = os.path.join(os.environ.get('PYTHONPATH'),'files')

import pandas as pd

humanVariants_path = os.path.join(FILES_PATH, 'test/uniprot_variants_test.csv')

humanVariants_df = pd.read_csv(humanVariants_path, sep='\t')

expressionSet = set()

for idx, row in humanVariants_df.iterrows():
    
    if row['Clinical Significance'] != '-':
        expressionSet.add(row['Clinical Significance'])

print(expressionSet)