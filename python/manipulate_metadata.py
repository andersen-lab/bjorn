#!/usr/bin/env python
import os
import json
import argparse
import pandas as pd
from fuzzywuzzy import fuzz

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--inputmetadata",
                        type=str,
                        required=True,
                        help="Input filepath metadata csv")
parser.add_argument("-o", "--outfp",
                        type=str,
                        required=True,
                        help="Output filepath")


args = parser.parse_args()
metadata_filepath = args.inputmetadata
output_filepath = args.outfp
#read in the overall metadata file 
meta_df = pd.read_csv(metadata_filepath)
meta_fasta = meta_df['fasta_hdr'].tolist()

#get all .fasta files names
all_fasta = [os.path.splitext(filename)[0] for filename in os.listdir(output_filepath) if os.path.splitext(filename)[1] == '.fasta']

match_indices = [-1] * len(all_fasta)
#first pass is looking for exact match between fileaname and metadata id
for count, fasta in enumerate(all_fasta):
    if fasta in meta_df['ID'].tolist():
        match_indices[count] = meta_df['ID'].tolist().index(fasta)

#if we can't find it we turn to fuzzy matching
for count, fasta in enumerate(all_fasta):
    #we already found an exact match
    if match_indices[count] != -1:
        continue
    high_score = 0
    high_index = 0
    for i,meta in enumerate(meta_fasta):
        if str(meta).lower() == 'nan':
            continue
        #compare and scroe string similarity
        current = fuzz.ratio(fasta,meta)
        if current > high_score:
            high_score = current
            high_index = i
    match_indices[count] = high_index

remove_indices = [x for x in range(len(meta_fasta)) if x not in match_indices]
new_row = [0]*len(meta_df)
for i,f in zip(match_indices, all_fasta):
    new_row[i] = f
meta_df['fasta'] = new_row
#remove what we don't want from metadata
meta_df.drop(meta_df.index[remove_indices], inplace=True)

c = []
d = []
l = []
#manipulate metadata for geo locations accordingly
for index, row in meta_df.iterrows():
    country = ""
    division = ""
    location = ""
    location_list = row['location'].split('/')[1:]
    if len(location_list) > 2:
        country = location_list[0]
        division = location_list[1]
        location = location_list[2]
    elif len(location_list) == 2:
        country = location_list[0]
        location = location_list[1]
    c.append(country)
    d.append(division)
    l.append(location)

#drop and add columns
meta_df.drop(['location', 'zipcode', 'authors', 'originating_lab', 'ID','gb_accession', \
'percent_coverage_cds', 'avg_depth', 'fasta_hdr'], axis=1, inplace=True)
meta_df['country'] = c
meta_df['division'] = d
meta_df['location'] = l
#putting in filler columns to make sure it doesn't crash
meta_df['country_normed'] = c
meta_df['division_normed'] = d
meta_df['location_normed'] = l
meta_df['date_submitted'] = ['']*len(meta_df)
meta_df['pangolin_lineage'] = ['']*len(meta_df)
meta_df['pangolin_version'] = ['']*len(meta_df)
meta_df['clade'] = ['']*len(meta_df)
meta_df.rename(columns={'collection_date':'date_collected', \
    'gisaid_accession':'accession_id'}, inplace=True)

for i,fasta in enumerate(all_fasta):
    with open(os.path.join(output_filepath, fasta+'.fasta'), "r") as file:
        first_line = str(file.readline()).strip().replace('>','')
    temp_df = meta_df[meta_df['fasta'] == fasta]
    temp_df.drop(['fasta'], axis=1, inplace=True)
    temp_df['strain'] = first_line
    temp_df.to_csv(os.path.join(output_filepath, "%s.tsv.gz" %fasta), "\t")
