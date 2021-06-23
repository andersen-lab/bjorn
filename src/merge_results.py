#!/usr/bin/env python
import os
import sys
import glob
import argparse
import time
import json
import pandas as pd
from path import Path
import bjorn_support as bs
import numpy as np

# COLLECTING USER PARAMETERS
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputmutations",
                        type=str,
                        required=True,
                        help="Input mutations csv")
parser.add_argument("-m", "--inputmeta",
                        type=str,
                        required=True,
                        help="Input metadata")
parser.add_argument("-o", "--outfp",
                        type=str,
                        required=True,
                        help="Output filepath")
parser.add_argument("-u", "--unknownvalue",
                        type=str,
                        required=True,
                        help="Unknown value")
parser.add_argument("-n", "--mindate",
                        type=str,
                        required=True,
                        help="Minimum date")
parser.add_argument("-g", "--geojson",
                        type=str,
                        required=True,
                        help="GeoJSON prefix")
parser.add_argument("-t", "--currentdate",
                        type=str,
                        required=True,
                        help="Current date")

args = parser.parse_args()
input_mut = args.inputmutations
input_metadata = args.inputmeta
out_fp = args.outfp
unknown_val = args.unknownvalue
min_date = args.mindate
geojson_prefix = args.geojson
current_datetime = args.currentdate
date_modified = '-'.join(current_datetime.split('-')[:3]) + '-' + ':'.join(current_datetime.split('-')[3:])
max_date = '-'.join(current_datetime.split('-')[:3])

api_data_fp = out_fp            # Output JSON

countries_fp = '{}/gadm_countries.json'.format(geojson_prefix)
divisions_fp = '{}/gadm_divisions.json'.format(geojson_prefix)
locations_fp = '{}/gadm_locations.json'.format(geojson_prefix)

with open(countries_fp) as f:
    countries = json.load(f)
with open(divisions_fp) as f:
    divisions = json.load(f)
with open(locations_fp) as f:
    locations = json.load(f)

# concat with pd
muts = pd.read_csv(input_mut, dtype=str)
# ignore mutations found in non-coding regions
muts = muts.loc[~(muts['gene']=='Non-coding region')]
# fuse with metadata
print(f"Fusing with metadata...")
start = time.time()
meta = pd.read_csv(input_metadata, sep='\t', compression='gzip')
muts = pd.merge(muts, meta, left_on='idx', right_on='strain')
fuse_time = time.time() - start
# clean geographic information
start = time.time()
muts['country'] = muts['country'].astype(str)
muts['country_lower'] = muts['country'].str.lower()
muts['country_normed'] = muts['country_normed'].astype(str)
muts['country_normed_lower'] = muts['country_normed'].str.lower()
muts['division'] = muts['division'].astype(str)
muts['division_lower'] = muts['division'].str.lower()
muts['division_normed'] = muts['division_normed'].astype(str)
muts['division_normed_lower'] = muts['division_normed'].str.lower()
muts['location'] = muts['location'].astype(str)
muts['location_lower'] = muts['location'].str.lower()
muts['location_normed'] = muts['location_normed'].astype(str)
muts['location_normed_lower'] = muts['location_normed'].str.lower()
# clean time information
muts['tmp'] = muts['date_collected'].str.split('-')
muts = muts[muts['tmp'].str.len()>=2]
muts.loc[muts['tmp'].str.len()==2, 'date_collected'] += '-15'
muts['date_collected'] = pd.to_datetime(muts['date_collected'], errors='coerce')
muts['date_collected'] = muts['date_collected'].astype(str)
# Date conditions
muts = muts[muts['date_collected'] <= max_date]
muts = muts[muts['date_collected'] > min_date]

# Filter out sequences with collection dates after submission dates
muts = muts.loc[muts['date_collected']<=muts['date_submitted']]

# rename field names
muts.rename(
    columns={
        'country': 'country_original',
        'division': 'division_original',
        'location': 'location_original',
        'country_lower': 'country_original_lower',
        'division_lower': 'division_original_lower',
        'location_lower': 'location_original_lower',
        'country_normed': 'country',
        'division_normed': 'division',
        'location_normed': 'location',
        'country_normed_lower': 'country_lower',
        'division_normed_lower': 'division_lower',
        'location_normed_lower': 'location_lower',
        'del_len': 'change_length_nt'
    }, inplace=True)

# final cleaning (missing values)
muts.loc[muts['location']=='unk', 'location'] = unknown_val
muts.loc[muts['division']==muts['country'], 'division'] = unknown_val
muts.fillna(unknown_val, inplace=True)
# generate json
meta_info = [
        'strain', 'accession_id',
        'date_modified', 'date_collected', 'date_submitted',
        'country_id', 'country', 'country_original', 'country_lower', 'country_original_lower',
        'division_id', 'division', 'division_original', 'division_lower', 'division_original_lower',
        'location_id', 'location', 'location_original', 'location_lower', 'location_original_lower',
        'pangolin_lineage', 'pangolin_version',
        'clade'
]

muts_info = [
    'type', 'mutation', 'gene',
    'ref_codon', 'pos', 'alt_codon',
    'is_synonymous',
    'ref_aa', 'codon_num', 'alt_aa',
    'absolute_coords',
    'change_length_nt', 'is_frameshift',
    'deletion_codon_coords'
]
muts = muts[~(muts['gene'].isin(['5UTR', '3UTR']))]
muts['date_modified'] = date_modified
muts['country_id'] = muts['country'].apply(lambda x: countries.get(x, unknown_val)).astype(str)
muts['tmp_info1'] = muts['country'] + '-' + muts['division']
muts['division_id'] = muts['tmp_info1'].apply(lambda x: divisions.get(x, unknown_val)).astype(str)
muts['tmp_info2'] = muts['country'] + '-' + muts['division'] + '-' + muts['location']
muts['location_id'] = muts['tmp_info2'].apply(lambda x: locations.get(x, unknown_val)).astype(str)
# muts['country_id'] = muts['country'].apply(lambda x: countries.get(x, unknown_val))
# muts['division_id'] = muts['division'].apply(lambda x: divisions.get(x, unknown_val))
# muts['location_id'] = muts['location'].apply(lambda x: locations.get(x, unknown_val))
muts = muts.drop_duplicates(subset=['accession_id', 'mutation'])
preprocess_time = time.time() - start

# If deletions not in chunk add columns
del_columns = ['is_frameshift', 'change_length_nt', 'deletion_codon_coords', 'absolute_coords']
muts_columns = muts.columns.tolist()
for i in del_columns:
    if i not in muts_columns:
        muts[i] = np.nan

# GENERATE JSON DATA MODEL
start = time.time()
(
    muts.groupby(meta_info, as_index=True)
             .apply(lambda x: x[muts_info].to_dict('records'))
             .reset_index()
             .rename(columns={0:'mutations'})
             .to_json(
                 api_data_fp,
                 orient='records',
                 lines = True,
                 compression = "gzip"
             )
 )
