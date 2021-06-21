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


# COLLECTING USER PARAMETERS
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputdir",
                        type=str,
                        required=True,
                        help="Input directory containing chunked mutation results in csv files")
parser.add_argument("-m", "--inputmeta",
                        type=str,
                        required=True,
                        help="Input filepath containing chunked metadata for each sample")
parser.add_argument("-o", "--outfp",
                        type=str,
                        required=True,
                        help="Output filepath")
parser.add_argument("-t", "--time",
                        type=str,
                        required=True,
                        help="Current datetime")

args = parser.parse_args()
input_dir = args.inputdir
meta_fp = args.inputmeta
out_fp = args.outfp
current_datetime = args.time

with open('config.json', 'r') as f:
    config = json.load(f)
is_test = config['feed_test']
out_dir = config['out_dir']
log_fp = out_dir + '/' + 'logs' + '/' + config['log_file'] + '_' + current_datetime + '.txt'
unknown_val = config['unknown_value']
date = '-'.join(current_datetime.split('-')[:3])
min_date = config['min_date']
date_modified = '-'.join(current_datetime.split('-')[:3]) + '-' + ':'.join(current_datetime.split('-')[3:])
api_data_fp = config['outbreak_fp']#/valhalla/gisaid/new_api_data.json.gz'
meta_data_fp = config['meta_outbreak_fp']#'/valhalla/gisaid/new_genomics_metadata.json'
countries_fp = config['countries_fp']#'/home/al/data/geojsons/gadm_countries.json'
divisions_fp = config['divisions_fp']#'/home/al/data/geojsons/gadm_divisions.json'
locations_fp = config['locations_fp']#'/home/al/data/geojsons/gadm_locations.json'
with open(countries_fp) as f:
    countries = json.load(f)
with open(divisions_fp) as f:
    divisions = json.load(f)
with open(locations_fp) as f:
    locations = json.load(f)
# write metadata info to json file
metadata = {'date_modified': date_modified}
if not is_test:
    with open(meta_data_fp, 'w') as fp:
        json.dump(metadata, fp)
# fetch list of mutation csv filepaths
mut_fps = glob.glob(f"{input_dir}/*.mutations.csv")
# concat with pd
muts = pd.concat((pd.read_csv(fp, dtype=str) for fp in mut_fps))
# ignore mutations found in non-coding regions
muts = muts.loc[~(muts['gene']=='Non-coding region')]
# fuse with metadata
print(f"Fusing with metadata...")
start = time.time()
meta = pd.read_csv(meta_fp, sep='\t', compression='gzip')
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
muts = muts[muts['date_collected']<date]
muts = muts[muts['date_collected']>min_date]
# filter out sequences with collection dates after submission dates
muts = muts.loc[muts['date_collected']<=muts['date_submitted']]
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
# generate json
meta_info = [
        'strain', 'accession_id',
        'date_modified', 'date_collected', 'date_submitted',
        'country_id', 'country', 'country_original', 'country_lower', 'country_original_lower',
        'division_id', 'division', 'division_original', 'division_lower', 'division_original_lower',
        'location_id', 'location', 'location_original', 'location_lower', 'location_original_lower',
#         'submitting_lab', 'originating_lab',
#         'authors', 
        'pangolin_lineage', 'pangolin_version',
        'clade', 
#         'nextstrain_clade',
#         'gisaid_epi_isl', 'genbank_accession',
#         'purpose_of_sequencing',
            ]

muts_info = ['type', 'mutation', 'gene', 
             'ref_codon', 'pos', 'alt_codon', 
             'is_synonymous', 
             'ref_aa', 'codon_num', 'alt_aa', 
             'absolute_coords', 
             'change_length_nt', 'is_frameshift',
             'deletion_codon_coords']
muts = muts[~(muts['gene'].isin(['5UTR', '3UTR']))]
muts['date_modified'] = date_modified
# attach geolocation IDs based on GADM
muts['country_id'] = muts['country'].apply(lambda x: countries.get(x, unknown_val)).astype(str)
muts['tmp_info1'] = muts['country'] + '-' + muts['division']
muts['division_id'] = muts['tmp_info1'].apply(lambda x: divisions.get(x, unknown_val)).astype(str)
muts['tmp_info2'] = muts['country'] + '-' + muts['division'] + '-' + muts['location']
muts['location_id'] = muts['tmp_info2'].apply(lambda x: locations.get(x, unknown_val)).astype(str)
# drop record(s) with same sequence ID and mutation (redundancy cleaning)
muts = muts.drop_duplicates(subset=['accession_id', 'mutation'])
preprocess_time = time.time() - start
# GENERATE JSON DATA MODEL
start = time.time()
if not is_test:
    (muts.groupby(meta_info, as_index=True)
             .apply(lambda x: x[muts_info].to_dict('records'))
             .reset_index()
             .rename(columns={0:'mutations'})
             .to_json(api_data_fp,
                      orient='records'))
    io_time = time.time() - start
    start = time.time()
    gzip_cmd = f"gzip -f {api_data_fp}"
    # bs.run_command(gzip_cmd)
    api_data_fp += '.gz'
    gzip_time = time.time() - start
    # print(f'Execution time: {end - start} seconds')
    # upload to gcloud
    upload_cmd = f"/home/al/code/google-cloud-sdk/bin/gsutil -m cp {api_data_fp} gs://andersen-lab_temp/outbreak_genomics/"
    # bs.run_command(upload_cmd)
    # update gcloud access rights
    upload_cmd = f"/home/al/code/google-cloud-sdk/bin/gsutil -m cp {meta_data_fp} gs://andersen-lab_temp/outbreak_genomics/"
    # bs.run_command(upload_cmd)
    # send auto-slack message about it? (nah, too much)
    access_cmd = f"/home/al/code/google-cloud-sdk/bin/gsutil acl ch -R -u AllUsers:R gs://andersen-lab_temp/outbreak_genomics/*"
    # bs.run_command(access_cmd)
    api_data_fn = api_data_fp.split('/')[-1]
    stat_cmd = f"/home/al/code/google-cloud-sdk/bin/gsutil stat gs://andersen-lab_temp/outbreak_genomics/{api_data_fn}"
else:
    io_time = "N/A"
    gzip_time = "N/A"
# upload_stats = bs.run_command_log(stat_cmd)
# GENERATE CSV DATA MODEL
num_records = muts.drop_duplicates(subset=meta_info).shape[0]
num_ids = muts['accession_id'].unique().shape[0]
muts.to_csv(out_fp, index=False)
# Data logging
with open(log_fp, 'w') as f:
    f.write(f"Number of records: {num_records}\n")
    f.write(f"Number of unique accession IDs: {num_ids}\n")
    f.write(f"Metadata Fusion Execution time: {fuse_time} seconds\n")
    f.write(f"Data Preprocessing Execution time: {preprocess_time} seconds\n")
    f.write(f"IO Execution time: {io_time} seconds\n")
    f.write(f"Gzip Execution time: {gzip_time} seconds\n")
    # f.write(f"GCloud upload statistics: \n {upload_stats}")
print(f"Transfer Complete. All results saved in {out_dir}")