import os
import sys
import argparse
import time
import json
import pandas as pd
from path import Path


import bjorn_support as bs
import mutations as bm
import data as bd


# COLLECTING USER PARAMETERS
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input",
                        type=str,
                        required=True,
                        help="FASTA filepath containing aligned sequences")
parser.add_argument("-m", "--meta",
                        type=str,
                        required=True,
                        help="Gzipped TSV filepath containing sequence metadata")
parser.add_argument("-r", "--patient-zero",
                        type=str,
                        default="NC_045512.2",
                        help="Sample name of reference sequence")
parser.add_argument("-o", "--outfp",
                        type=str,
                        required=True,
                        help="Output directory")
args = parser.parse_args()
alignment_filepath = args.input
gisaid_meta = args.meta
patient_zero = args.patient_zero
out_fp = Path(args.outfp)

with open('config.json', 'r') as f:
    config = json.load(f)
date = config['date']
patient_zero = config['patient_zero']
data_src = config['data_source']
min_date = config['min_date']
unknown_val = config['unknown_value']
countries_fp = config['countries_fp']
divisions_fp = config['divisions_fp']
locations_fp = config['locations_fp']


print(f"Loading alignment file at {alignment_filepath}")
t0 = time.time()
msa_data = bs.load_fasta(alignment_filepath, is_aligned=True, is_gzip=False)
msa_load_time = time.time() - t0
print(f"Identifying substitution-based mutations...")
t0 = time.time()
subs, _ = bm.identify_replacements_per_sample(msa_data, 
                                              gisaid_meta,  
                                              gene2pos=bd.GENE2POS, 
                                              data_src=data_src,
                                              min_seq_len=20000,
                                              patient_zero=patient_zero
                                            #   test=is_test
                                              )
subs_time = time.time() - t0
print(f"Identifying deletion-based mutations...")
t0 = time.time()
dels, _ = bm.identify_deletions_per_sample(msa_data, 
                                           gisaid_meta,  
                                           gene2pos=bd.GENE2POS, 
                                           data_src=data_src, 
                                           min_del_len=1,
                                           max_del_len=500,
                                           min_seq_len=20000,
                                           patient_zero=patient_zero
                                        #    test=is_test
                                           )
dels_time = time.time() - t0
# QC FILTER: remove seqs with >500 nt deletions
# dels = dels.loc[dels['del_positions'].str.len()<500]
print(subs.shape)
print(dels.shape)
muts = pd.concat([subs, dels])
muts['is_synonymous'] = False
muts.loc[muts['ref_aa']==muts['alt_aa'], 'is_synonymous'] = True
print(muts.shape)
# clean geographic information
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
muts = muts[muts['date_collected']<date]
muts = muts[muts['date_collected']>min_date]
# rename field names
muts.rename(columns={
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
# muts = muts.astype(str) TAKES FOREVER
# muts_filename = alignment_filepath.replace('.aligned.fasta', f'_{date}.mutations.csv')
muts.to_csv(out_fp, index=False)
print(f"Mutations extracted from {alignment_filepath} and saved in {out_fp}\n")