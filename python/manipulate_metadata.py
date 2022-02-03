#!/usr/bin/env python
import re
import os
import sys
import gzip
import json
import argparse
import pandas as pd
from rapidfuzz import fuzz

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--inputmetadata",
                        type=str,
                        required=True,
                        help="Input filepath metadata csv")

parser.add_argument("-l", "--lineagereport",
                        type=str, 
                        required=False,
                        help="Input filepath lineage report csv")

parser.add_argument("-o", "--outfp",
                        type=str,
                        required=True,
                        help="Output filepath")
parser.add_argument("-f", "--filenames",
                        required=True,
                        nargs='+',
                        help="Filenames of interest")

parser.add_argument("-j", "--jobid",
                        required=True,
                        help="Filenames of interest")



#print("Opening metadata")
args = parser.parse_args()
metadata_filepath = args.inputmetadata
output_filepath = args.outfp
lineage_filepath = args.lineagereport

#read in the overall metadata file
meta_df = pd.read_csv(metadata_filepath)  

all_filenames = args.filenames

#get all .fasta files names
all_fasta = [os.path.basename(os.path.splitext(filename)[0]).replace('.fasta','') for filename in all_filenames]

#get all header strains
job_num =args.jobid

#our ideal match is on strain, else we'll match on filename
strains=[]
for i,fasta in enumerate(all_fasta):
    with gzip.open(os.path.join(output_filepath,fasta +'.fasta.gz'), "r") as file:
        first_line = file.readline().decode("utf-8").strip().replace('>','')
    strains.append(str(first_line).strip())

fasta_hdr = meta_df['fasta_hdr'].tolist()

#these are the fasta files who's hdrs aren't in the metadata file
missing_fasta = [fasta for (strain,fasta) in zip(strains,all_fasta) if strain not in fasta_hdr]
keep_indices = []
scores = []
for fasta in missing_fasta:
    high_score = 0
    high_index = 0
    fasta = str(fasta)
    for index, row in meta_df.iterrows():
        row_val = str(row['ID'])
        score = fuzz.token_set_ratio(row_val, fasta)
        #we need the score to be the highest and more than 97% match
        if float(score) > high_score and float(score) > 97:
            high_score = float(score)
            high_index = index
            #print(score, high_score, row['ID'], high_index, str(fasta))
    
    #we didn't find a good match
    if float(high_score) < 97:
        continue
    keep_indices.append(high_index)
    scores.append(high_score)

drop_indices = [value for value in list(range(0,len(meta_df))) if value not in keep_indices]
temp_df = meta_df[meta_df['fasta_hdr'].isin(strains)]
temp_df_2 = meta_df.drop(meta_df.index[drop_indices])
meta_df = pd.concat([temp_df, temp_df_2], axis=0)

if lineage_filepath:
    lineage_df = pd.read_csv(lineage_filepath)
    taxon=lineage_df['taxon'].tolist()
    taxon_added= []
    for index, row in meta_df.iterrows():
        taxon_string=""
        high_score = 0
        for t in taxon:
            ratio = fuzz.token_set_ratio(str(row['fasta_hdr']),str(t))
            if ratio > high_score:
                high_score = ratio
                taxon_string = t
         
        taxon_added.append(taxon_string)
    meta_df['taxon'] = taxon_added
    meta_df = meta_df.merge(lineage_df, how='left', on='taxon')
c = []
d = []
l = []

#drop and add columns
meta_df.drop(['ID','gb_accession', \
'percent_coverage_cds', 'avg_depth'], axis=1, inplace=True)
meta_df['date_submitted'] = ['']*len(meta_df)


meta_df.rename(columns={'collection_date':'date_collected', \
    'gisaid_accession':'accession_id', 'fasta_hdr':'strain', 'location':'locstring',\
    "lineage":"pangolin_lineage"}, inplace=True)


meta_df.to_csv(os.path.join(output_filepath, "%s.tsv.gz" %job_num), "\t")


#concatenate all these things together
all_fasta_files = [os.path.join(output_filepath, fasta+'.fasta.gz') for fasta in all_fasta]
with gzip.open(os.path.join(output_filepath,str(job_num) + '.fasta.gz'), 'w') as outfile:
    for fname in all_fasta_files:
        with gzip.open(fname, "r") as infile:
            for line in infile:
                outfile.write(line)
        outfile.write(b'\n')
        os.system("rm %s" %fname)
