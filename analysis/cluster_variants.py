import os
import sys
import pandas as pd
from utils import function_timer
from collections import defaultdict
from matplotlib import pyplot as plt
from pickleconfig import PickleJar, PickleShelf

jar = PickleJar()

ROOT_DIR = os.path.join(os.environ.get('PYTHONPATH'))

CLUSTER_MAP_PATH = os.path.join(ROOT_DIR, 'files/uniclust30_2016_03_cluster_mapping.tsv')
VARIANTS_PATH = os.path.join(ROOT_DIR, 'synced_files/uniprot_missense_variants/benign_mutations_with_lengths.csv')

VARS_PER_CLUSTER = defaultdict(list)

UNMAPPABLE = []
cluster_seeds = []

variants_df = pd.read_csv(VARIANTS_PATH, sep='\t', header=0)

def build_cluster_dict():
    
    cluster_df = None

    try:
        print('reading reference dataframe from cache')
        cluster_df = pd.read_pickle('../dataframe_pickles/uniclust30_2016_03_cluster_mapping.pkl')
    except FileNotFoundError:
        print('cache not found')
        
        cluster_df = pd.read_csv(CLUSTER_MAP_PATH, sep='\t', header=0)
        
        print('headers:', cluster_df.columns)

        print('caching reference dataframe')
        # cache dataframe
        cluster_df.to_pickle('../dataframe_pickles/uniclust30_2016_03_cluster_mapping.pkl')
    
    clusters = {}
    
    for row in cluster_df.itertuples():
        
        if int(getattr(row, 'Index')) % 100000 == 0:
            print('caching line: ', getattr(row, 'Index'))

        clusters[str(getattr(row,'member'))] = str(getattr(row,'cluster'))
        
    return clusters

CLUSTER_DICT = jar.pickle(build_cluster_dict,'CLUSTER_DICT')

for idx, row in variants_df.iterrows():
    
    if idx % 100 == 0:
        print(f'processing line: {idx}\r', end='')
    
    uniprot = row['uniprot']
    WT = row['WT']
    Mut = row['Mut']
    position = row['protein_position']
    mutation = WT + position + Mut
    
    try:
        cluster = CLUSTER_DICT[uniprot]
        VARS_PER_CLUSTER[cluster].append((uniprot, mutation))
        cluster_seeds.append(cluster)
    except KeyError:
        UNMAPPABLE.append((uniprot, mutation))
        cluster_seeds.append('-')

variants_df['cluster'] = cluster_seeds

variants_df.to_csv('../synced_files/uniprot_missense_variants/benign_mutations_with_clusters.csv', sep='\t', index=False)

