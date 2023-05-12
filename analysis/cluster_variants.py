import os
import sys
import pandas as pd
from utils import function_timer
from collections import defaultdict
from matplotlib import pyplot as plt
from pickleconfig import PickleJar, PickleShelf

jar = PickleJar()

ROOT_DIR = os.path.join(os.environ.get('PYTHONPATH'))

CLUSTER_MAP_PATH = os.path.join(ROOT_DIR, 'files/uniclust30_2018_08.tsv')
VARIANTS_PATH = os.path.join(ROOT_DIR, 'synced_files/vep_variants/vep_benign_variants.csv')

VARS_PER_CLUSTER = defaultdict(list)

UNMAPPABLE = []
cluster_seeds = []

variants_df = pd.read_csv(VARIANTS_PATH, sep='\t', header=0)

def build_cluster_dict():
    
    cluster_df = pd.read_csv(CLUSTER_MAP_PATH, sep='\t', header=0)
    
    clusters = {}
    
    for row in cluster_df.itertuples():
        
        if int(getattr(row, 'Index')) % 100000 == 0:
            print('caching line: ', getattr(row, 'Index'))

        clusters[str(getattr(row,'member'))] = str(getattr(row,'cluster'))
        
    return clusters

CLUSTER_DICT = jar.pickle(build_cluster_dict,'CLUSTER_DICT_2018')

for idx, row in variants_df.iterrows():
    
    if idx % 1000 == 0:
        print(f'processing line: {idx}\r', end='')
    
    uniprot = row['uniprot']
    
    try:
        cluster = CLUSTER_DICT[uniprot]
        cluster_seeds.append(cluster)
    except KeyError:
        cluster_seeds.append('-')

variants_df['cluster'] = cluster_seeds

variants_df.to_csv('../synced_files/vep_variants/vep_benign_variants_with_clusters.csv', sep='\t', index=False)

