import pandas as pd
from itertools import combinations

source_df = pd.read_csv('./ramp2_source_unique.txt', sep=' ', header=None)

source_names = source_df.iloc[:,1].unique()

source_dfs_dict = {
    'genes': {},
    'compounds': {}
}

intersects_dict = {
    'genes': {},
    'compounds': {}
}

for source_name in source_names:
    source_dfs_dict['genes'][source_name] = source_df.loc[(source_df.iloc[:,1] == source_name) & (source_df.iloc[:,0].str.contains('G'))]
    source_dfs_dict['genes'][source_name].set_index(0, inplace=True)
    print(source_dfs_dict['genes'][source_name].head())
    source_dfs_dict['compounds'][source_name] = source_df.loc[(source_df.iloc[:,1] == source_name) & (source_df.iloc[:,0].str.contains('C'))]
    source_dfs_dict['compounds'][source_name].set_index(0, inplace=True)
    print(source_dfs_dict['compounds'][source_name].head())
    intersects_dict['genes'][source_name] = len(source_dfs_dict['genes'][source_name].index)
    intersects_dict['compounds'][source_name] = len(source_dfs_dict['compounds'][source_name].index)

key = source_names[0]
genes_intersection = None
compounds_intersection = None

for source_name, index in enumerate(source_names):


    for other_source_name in source_names:
        if other_source_name is not source_name:


    if index == 0:
        genes_intersection = source_dfs_dict['genes'][source_name].index
        compounds_intersection = source_dfs_dict['compounds'][source_name].index
    else:
        genes_intersection = genes_intersection.intersection(source_dfs_dict['genes'][source_name].index)
        compounds_intersection = compounds_intersection.intersection(source_dfs_dict['compounds'][source_name].index)

        for _key in intersects_dict['genes'].keys():
            intersects_dict['genes'][_key] = intersects_dict['genes'][_key] - len(genes_intersection)
            intersects_dict['compounds'][_key] = intersects_dict['compounds'][_key] - len(compounds_intersection)

        key = f'{key} + {source_name}'

        intersects_dict['genes'][key] = len(genes_intersection)
        intersects_dict['compounds'][key] = len(compounds_intersection)
        
# print(pd.Index(source_dfs_dict['hmdb']))
# intersection = source_dfs_dict['hmdb'].index.intersection(source_dfs_dict['reactome'].index)

# print(intersection)

# intersections_dict = {}

