#!/usr/bin/env python
import os
import re
import sys
import glob
import argparse
import time
import json
import pandas as pd
from path import Path
import bjorn_support as bs
from mappings import COUNTY_CORRECTIONS
import numpy as np
from rapidfuzz import fuzz

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

try:
    meta = pd.read_csv(input_metadata, sep='\t')
    if len(meta.index) == 0:
        raise "empty dataframe"
except Exception as e:
    pd.DataFrame().to_json(out_fp, orient='records', lines=True)
    sys.exit(0)

#parse apart our expected geolocs
countries_std_file = open(os.path.join(geojson_prefix, "gadm_countries.json"),'r')
divisions_std_file = open(os.path.join(geojson_prefix, "gadm_divisions.json"),'r')
locations_std_file = open(os.path.join(geojson_prefix, "gadm_locations.json"), 'r')

countries_std = json.load(countries_std_file)
divisions_std = json.load(divisions_std_file)
locations_std = json.load(locations_std_file)
std_locs = [countries_std, divisions_std, locations_std]


# normalize date information
meta['tmp'] = meta['date_collected'].str.split('-')
meta = meta[meta['tmp'].str.len()>=2]
meta.loc[meta['tmp'].str.len()==2, 'date_collected'] += '-15'
meta['date_collected'] = pd.to_datetime(meta['date_collected'], errors='coerce')
meta['date_collected'] = meta['date_collected'].astype(str)
meta = meta[meta['date_collected'] <= max_date]
meta = meta[meta['date_collected'] > min_date]
meta['date_modified'] = date_modified
#make an empy list for the data to be held in
std_loc = [[0]*3]*len(meta) 

#country location variable
country_index_location = 0

#we just use this loop to find out which string is the country, then go from there
for index, row in meta.iterrows():
    locstring = str(row['locstring']) 
    loc_list = locstring.split("/")
    if len(loc_list) == 0:
        std_loc[index][0] = unknown_val
        std_loc[index][1] = unknown_val
        std_loc[index][2] = unknown_val
        continue
    #we find the country first
    compare_dict = std_locs[0]
    high_ratio = 0
    country_string = ''
    
    for i,loc in enumerate(loc_list): 
        for key in compare_dict.keys():
            ratio = fuzz.ratio(key.lower(),loc.lower())
            if ratio > high_ratio:
                high_ratio = ratio
                country_string = key
                country_index_location = i
        for value in compare_dict.values():
            ratio = fuzz.ratio(value.lower(),loc.lower())
            if ratio > high_ratio:
                high_ratio = ratio
                country_string = list(compare_dict.keys())[list(compare_dict.values()).index(value)]
                country_index_location = i   
 
    high_ratio = 0 
    high_string = ''
    div_or_loc = 0
    #next we figure out if we have a division or location in the next slot
    try:
        next_string = loc_list[country_index_location+1]
    except:
        std_loc[index][0] = country_string
        std_loc[index][1] = unknown_val
        std_loc[index][2] = unknown_val
        continue

    for i,compare_dict in enumerate(std_locs[1:]): 
        for key in compare_dict.keys():
            #make sure the country prefix matches what we already have
            country_key = re.split("-", key)[0]
            if country_key.lower() not in country_string.lower():
                continue
            #get rid of any country level information
            key = re.split("-",key)[1:]
            key = " ".join(key)
            ratio = fuzz.ratio(key.lower(),next_string.lower())
            if ratio > high_ratio:
                high_ratio = ratio
                high_string = key
                div_or_loc = i
        for value in compare_dict.values():
            ratio = fuzz.ratio(value.lower(),next_string.lower())
            if ratio > high_ratio:
                high_ratio = ratio
                high_string = list(compare_dict.keys())[list(compare_dict.values()).index(value)]
                div_or_loc = i

    #we don't have a division only a location
    if div_or_loc == 1:
        std_loc[index][0] = country_string
        std_loc[index][1] = unknown_val
        std_loc[index][2] = high_string
        continue

    #we have a location, we suspect we have a division
    else:
        high_ratio = 0
        location_string = ''
        compare_dict = std_locs[2]
        
        #we have a division but no location
        try:
            next_string = loc_list[country_index_location+2]
        except:
            std_loc[index][0] = country_string
            std_loc[index][1] = high_string
            std_loc[index][2] = unknown_val
            continue
            
        for key in compare_dict.keys():
            #make sure the country prefix matches what we already have
            country_key = re.split("-", key)[0]
            if country_key.lower() not in country_string.lower():
                continue

            #get rid of any country level information
            key = re.split("-",key)[2:]
            key = " ".join(key)

            ratio = fuzz.ratio(key.lower(),next_string.lower())
            if ratio > high_ratio:
                high_ratio = ratio
                location_string = key
        
        for value in compare_dict.values():
            ratio = fuzz.ratio(value.lower(),next_string.lower())
            if ratio > high_ratio:
                high_ratio = ratio
                location_string = list(compare_dict.keys())[list(compare_dict.values()).index(value)]
        try: 
            #we found both a division and location
            std_loc[index][0] = country_string
            std_loc[index][1] = high_string
            std_loc[index][2] = location_string
        except:
            print(index, print(len(std_loc)), print(len(meta)))
            sys.exit(0)

loc_df = pd.DataFrame(std_loc, columns=['country','division','location'])
meta = pd.concat([meta,loc_df], axis=1)

#keep the county corrections mapping
for key, val in COUNTY_CORRECTIONS.items():
    meta.loc[:, 'location'] = meta['location'].str.replace(key, val)
meta['country_id'] = meta['country']
meta['division_id'] = meta['division']
meta['location_id'] = meta['location']
meta.replace({"country_id": std_locs[0]}, inplace=True)
meta.replace({"division_id": std_locs[1]}, inplace=True)
meta.replace({"location_id": std_locs[2]}, inplace=True)
meta['country_lower'] = meta['country'].str.lower()
meta['division_lower'] = meta['division'].str.lower()
meta['location_lower'] = meta['location'].str.lower()
meta.fillna(unknown_val, inplace=True)

meta_info = ['strain', 'accession_id', 'pangolin_lineage', 'date_collected', 'date_submitted', 'date_modified', 'country_id', 'division_id', 'location_id', 'country', 'division', 'location', 'country_lower', 'division_lower', 'location_lower']

muts = pd.read_csv(input_mut, dtype=str)
muts = muts[~(muts['gene'].isin(['5UTR', '3UTR']))]
# ignore mutations found in non-coding regions
muts = muts.loc[~(muts['gene']=='Non-coding region')]
# fuse with metadata
print(f"Fusing muts with metadata...")
# concat with pd
muts_info = [
    'type', 'mutation', 'gene',
    'ref_codon', 'pos', 'alt_codon',
    'is_synonymous',
    'ref_aa', 'codon_num', 'alt_aa',
    'absolute_coords',
    'change_length_nt', 'is_frameshift',
    'deletion_codon_coords'
]
# If deletions not in chunk add columns
del_columns = ['is_frameshift', 'change_length_nt', 'deletion_codon_coords', 'absolute_coords']
muts_columns = muts.columns.tolist()
for i in del_columns:
    if i not in muts_columns:
        muts[i] = np.nan
muts = muts.groupby('idx').apply(lambda x: x[muts_info].to_dict('records')).reset_index().rename(columns={0:'mutations'})
muts['mutations'] = muts['mutations'].map(lambda x: [{k:v for k,v in y.items() if pd.notnull(v)} for y in x])

muts = muts.rename(columns={'idx': 'strain'})
muts['strain'] = muts['strain'].str.strip()
meta['strain'] = meta['strain'].str.strip()
pd.merge(meta[meta_info], muts, on='strain', how='left').to_json(out_fp, orient='records', lines=True)
