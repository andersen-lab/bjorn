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
import uuid
#define global static variables
meta_info = ['strain', 'accession_id', 'pangolin_lineage', 'date_collected', 'date_submitted', \
    'date_modified', 'country_id', 'division_id', 'location_id', 'country', 'division', \
    'location', 'country_lower', 'division_lower', 'location_lower', 'zipcode']
del_columns = ['is_frameshift', 'change_length_nt', 'deletion_codon_coords', 'absolute_coords']
muts_info = ['type', 'mutation', 'gene', 'ref_codon', 'pos', 'alt_codon', 'is_synonymous', \
    'ref_aa', 'codon_num', 'alt_aa', 'absolute_coords', 'change_length_nt', 'is_frameshift','deletion_codon_coords']

#start of functions to support
def parse_jsonl(json_filepath: str):
    """
    Load and parse .jsonl file in read mode.    

    Parameters
    ----------
    json_filepath : str
        Path to the json file to parse.

    Returns
    -------
    json_data : list
        A list of dictionaries.

    """
    json_data = []
    with open(json_filepath, 'r') as jfile:
        for line in jfile:
            json_data.append(json.loads(line))
    return(json_data)

def nextstrain_replacement(nextstrain_filepath: str, meta: pd.DataFrame):
    """
    Opens the nextstrain file for standarization and replaces locstrings
    in the metadata dataframe.

    Parameters
    ---------
    nexstrain_filepath : str
        Filepath to the location of the standard nexstrain locations.
    
    meta : pd.DataFrame
        The metadata dataframe.

    Returns
    -------
    meta : pd.DataFrame
        The metadata dataframe
    """

    #first we open the nextstrain file and get the string conversions
    nextstrain = pd.read_csv(nextstrain_filepath, sep="\t", names=["original","replace"])
    nextstrain['original'] = [item.replace("/*","") for item in nextstrain['original']]
    nextstrain['replace'] = [item.replace("/*","") for item in nextstrain['replace']]
    nextstrain['original'] = [item[:-1] if item[-1] == '/' else item for item in nextstrain['original']]
    next_replace = dict(zip(list(nextstrain['original']), list(nextstrain['replace'])))
    
    locstrings = meta['locstring'].tolist()
    for i,loc in enumerate(locstrings):
        for key, value in next_replace.items():
            if str(key) in str(loc):
                new_str = str(loc).replace(str(key), str(value))
                locstrings[i] = new_str
        
    meta['locstring'] = locstrings
    return(meta)

def normalize_date(current_datetime: str, meta: pd.DataFrame, min_date: str):
    """
    Takes in the current date and normalizes it.
    
    Parameters
    ----------
    current_datetime : str    

    meta : pd.DataFrame
    
    min_date : str

    Returns
    -------
    meta : pd.DataFrame    

    """
    #reformat the string
    date_modified = '-'.join(current_datetime.split('-')[:3]) + '-' + ':'.join(current_datetime.split('-')[3:])
    max_date = '-'.join(current_datetime.split('-')[:3])
    
    # normalize date information
    meta['tmp'] = meta['date_collected'].str.split('-')
    meta = meta[meta['tmp'].str.len()>=2]
    meta.loc[meta['tmp'].str.len()==2, 'date_collected'] += '-15'
    meta['date_collected'] = pd.to_datetime(meta['date_collected'], errors='coerce')
    meta['date_collected'] = meta['date_collected'].astype(str)
    meta = meta[meta['date_collected'] <= max_date]
    meta = meta[meta['date_collected'] > min_date]
    meta['date_modified'] = date_modified
    return(meta)

def off_by_one_location(meta: pd.DataFrame, std_locs: list):
    """
    Takes in all proper countries, divisions, and location and 
    matches each value using fuzzy matching to the locstrings.

    Parameters
    ----------
    meta : pd.DataFrame
        The metadata dataframe.

    std_locs : list
        List of dictionaries for the locations.    

    Returns
    -------
    meta : pd.DataFrame
        The metadata dataframe.
    """
    unknown_val = 'unkn'

    #make an empy list for the data to be held in
    std_loc = []
    std_ids = []
    
    #country location variable
    country_index_location = 0

    #we just use this loop to find out which string is the country, then go from there
    for index, row in meta.iterrows():
        locstring = str(row['locstring']) 
        loc_list = locstring.split("/")
        loc_list = [loc.strip() for loc in loc_list]
        #in the event that we have no normal splits
        if len(loc_list) == 0:
            std_loc.append((unknown_val, unknown_val, unknown_val))
            continue
        
        #we find the country first
        compare_dict = std_locs[0]
        high_ratio = 0
        country_string = ''
        country_id_string = ''

        for i,loc in enumerate(loc_list): 
            #we look for the country in 'name' category
            for country in std_locs:
                country_name = country['name'] 
                ratio = fuzz.ratio(country_name.lower(),loc.lower()) 
                if ratio > high_ratio:
                    high_ratio = ratio
                    country_string = country_name
                    country_index_location = i
                    country_id_string = country['id']
            #we look for the country in 'id' category
            for country in std_locs:
                country_id = country['id']
                ratio = fuzz.ratio(country_id.lower(),loc.lower())
                if ratio > high_ratio:
                    high_ratio = ratio
                    country_string = country['name']
                    country_index_location = i   
                    country_id_string = country['id']
        high_ratio = 0 
        division_string = ''
        division_id_string = ''

        #next we figure out if we have a division or location in the next slot
        try:
            next_string = loc_list[country_index_location+1]
        except:
            std_loc.append((country_string, unknown_val, unknown_val))
            std_ids.append((country_id_string, unknown_val, unknown_val))
            continue
        
        #we go back and get the sub divisions for the country
        for country in std_locs:
            if country['name'] == country_string:
                division_list = country['sub']
        
        #we go through all the sub divisions to see if we can match one highly enough
        for i, division in enumerate(division_list): 
            division_name = division['name']
            division_id = division['id']
            ratio = fuzz.ratio(division_name.lower(),next_string.lower())
            
            if ratio > high_ratio and ratio > 90:
                high_ratio = ratio
                division_string = division_name
                division_id_string = division_id

            ratio = fuzz.ratio(division_id.lower(), next_string.lower())
            if ratio > high_ratio and ratio > 90:
                high_ratio = ratio
                division_string = division_name
                division_id_string = division_id        
        location_list = []
       
        #we didn't find a division, maybe the next string is a location
        if division_string == '':
            division_string = unknown_val
            for country in std_locs:
                if country['name'] == country_string:
                    division_list = country['sub']
                    for div in division_list:
                        location_list.extend(div['sub'])
        #we found a division
        else:
           
            for country in std_locs:
                if country['name'] == country_string:
                    division_list = country['sub']
                    for div in division_list:
                        if div['name'] == division_string:
                            location_list = div['sub']
            
            try:
                next_string = loc_list[country_index_location+2]
            #their is no next string
            except:
                std_loc.append((country_string, division_string, unknown_val))
                std_ids.append((country_id_string, division_id_string, unknown_val))
                continue

        #we looked for the division, now we look for the location
        high_ratio = 0
        location_string = ''
        location_id_string = ''
        #iterate the loications
        for location in location_list:
            location_name = location['name']
            ratio = fuzz.ratio(location_name.lower(),next_string.lower())
            
            if ratio > high_ratio and ratio > 70:
                high_ratio = ratio
                location_string = location_name
                location_id_string = location['id']
         
        if location_string != '':
            #we found both a division and location
            std_loc.append((country_string, division_string, location_string))  
            std_ids.append((country_id_string, division_id_string, location_id_string))
        else:
            std_loc.append((country_string, division_string, unknown_val))
            std_ids.append((country_id_string, division_id_string, unknown_val))
     
    loc_df = pd.DataFrame(std_loc, columns=['country','division','location'])
    ids_df = pd.DataFrame(std_ids, columns=['country_id', 'division_id', 'location_id'])
    if "location" in meta and "country" in meta and "division" in meta:
        meta.drop(["location", "country", "division"], inplace=True, axis=1)
    
    meta = pd.concat([meta, loc_df, ids_df], axis=1)
    return(meta)

def main():
    """
    Main function for use.
    """
    #COLLECTING USER PARAMETERS
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

    #handle the case of an empty dataframe
    try:
        meta = pd.read_csv(input_metadata, sep='\t')
        if len(meta.index) == 0:
            raise "empty dataframe"
    except Exception as e:
        pd.DataFrame().to_json(out_fp, orient='records', lines=True)
        sys.exit(0)

    #normalize the date
    meta = normalize_date(current_datetime, meta, min_date)

    #use nextstrain standard conversions to replace locstrings
    meta = nextstrain_replacement(os.path.join(geojson_prefix, "gisaid_geoLocationRules.tsv"), meta)

    #parse apart our expected geolocs
    std_locs = parse_jsonl(os.path.join(geojson_prefix, "gadm_transformed.jsonl"))

    meta = off_by_one_location(meta, std_locs)
    #keep the county corrections mapping
    #print(meta['location'])
    meta['location'] = [str(item) for item in meta['location'].tolist()]
    meta = meta.replace({"location":COUNTY_CORRECTIONS})

    meta['country_lower'] = meta['country'].str.lower()
    meta['division_lower'] = meta['division'].str.lower()
    meta['location_lower'] = meta['location'].str.lower()
    meta.fillna(unknown_val, inplace=True)

    muts = pd.read_csv(input_mut, dtype=str)
    muts = muts[~(muts['gene'].isin(['5UTR', '3UTR']))]

    #ignore mutations found in non-coding regions
    muts = muts.loc[~(muts['gene']=='Non-coding region')]
    #fuse with metadata
    print(f"Fusing muts with metadata...")
    # concat with pd

    # If deletions not in chunk add columns
    muts_columns = muts.columns.tolist()
    for i in del_columns:
        if i not in muts_columns:
            muts[i] = np.nan
    muts = muts.groupby('idx').apply(lambda x: x[muts_info].to_dict('records')).reset_index().rename(columns={0:'mutations'})
    muts['mutations'] = muts['mutations'].map(lambda x: [{k:v for k,v in y.items() if pd.notnull(v)} for y in x])

    muts = muts.rename(columns={'idx': 'strain'})
    muts['strain'] = muts['strain'].str.strip()
    meta['strain'] = meta['strain'].str.strip()
   
    templ = len(meta)
    
    #lets assign anything without an accession id to have a random string as accesion
    meta['accession_id'] = meta.apply(lambda row: row['accession_id'] if (str(row['accession_id']) != 'None') else (uuid.uuid4().hex[:15].upper()), axis=1)
   
    meta.drop_duplicates(subset='accession_id', keep="last", inplace=True) 
    meta = meta[meta['accession_id'] != 'None']
   
    if "zipcode" in meta.columns:
        td = pd.merge(meta[meta_info], muts, on='strain', how='left')
        if len(td) != len(meta):
            print(td)
            print(meta)
        pd.merge(meta[meta_info], muts, on='strain', how='left').to_json(out_fp, orient='records', lines=True)
    
    else:
        meta_info.remove("zipcode")
        pd.merge(meta[meta_info], muts, on='strain', how='left').to_json(out_fp, orient='records', lines=True)


if __name__ == "__main__":
    main()
