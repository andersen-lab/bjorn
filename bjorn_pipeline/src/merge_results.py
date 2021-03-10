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
parser.add_argument("-o", "--outfp",
                        type=str,
                        required=True,
                        help="Output filepath")

args = parser.parse_args()
input_dir = args.inputdir
out_fp = args.outfp

with open('config.json', 'r') as f:
    config = json.load(f)

unknown_val = config['unknown_value']
date = config['date']
date_modified = config['date_modified']
api_data_fp = config['outbreak_fp']#/valhalla/gisaid/new_api_data.json.gz'
meta_data_fp = config['meta_outbreak_fp']#'/valhalla/gisaid/new_genomics_metadata.json'
countries_fp = config['countries_fp']#'/home/al/data/geojsons/gadm_countries.json'
divisions_fp = config['divisions_fp']#'/home/al/data/geojsons/gadm_divisions.json'
locations_fp = config['locations_fp']#'/home/al/data/geojsons/gadm_locations.json'
# with open(countries_fp) as f:
#     countries = json.load(f)
# with open(divisions_fp) as f:
#     divisions = json.load(f)
# with open(locations_fp) as f:
#     locations = json.load(f)
with open(countries_fp) as f:
    countries = json.load(f)
with open(divisions_fp) as f:
    divisions = json.load(f)
with open(locations_fp) as f:
    locations = json.load(f)
# write metadata info to json file
metadata = {'date_modified': date_modified}
with open(meta_data_fp, 'w') as fp:
    json.dump(metadata, fp)
# fetch list of mutation csv filepaths
mut_fps = glob.glob(f"{input_dir}/*.mutations.csv")
# concat with pd
muts = pd.concat((pd.read_csv(fp, dtype=str) for fp in mut_fps))
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
muts['country_id'] = muts['country'].apply(lambda x: countries.get(x, unknown_val)).astype(str)
muts['tmp_info1'] = muts['country'] + '-' + muts['division']
muts['division_id'] = muts['tmp_info1'].apply(lambda x: divisions.get(x, unknown_val)).astype(str)
muts['tmp_info2'] = muts['country'] + '-' + muts['division'] + '-' + muts['location']
muts['location_id'] = muts['tmp_info2'].apply(lambda x: locations.get(x, unknown_val)).astype(str)
# muts['country_id'] = muts['country'].apply(lambda x: countries.get(x, unknown_val))
# muts['division_id'] = muts['division'].apply(lambda x: divisions.get(x, unknown_val))
# muts['location_id'] = muts['location'].apply(lambda x: locations.get(x, unknown_val))
muts = muts.drop_duplicates(subset=['accession_id', 'mutation'])
# GENERATE JSON DATA MODEL
start = time.time()
(muts.groupby(meta_info, as_index=True)
             .apply(lambda x: x[muts_info].to_dict('records'))
             .reset_index()
             .rename(columns={0:'mutations'})
             .to_json(api_data_fp,
                      orient='records',
                      compression='gzip'))
end = time.time()
print(f'Execution time: {end - start} seconds')
# upload to gcloud
upload_cmd = f"/home/al/code/google-cloud-sdk/bin/gsutil -m cp {api_data_fp} gs://andersen-lab_temp/outbreak_genomics/"
bs.run_command(upload_cmd)
# update gcloud access rights
upload_cmd = f"/home/al/code/google-cloud-sdk/bin/gsutil -m cp {meta_data_fp} gs://andersen-lab_temp/outbreak_genomics/"
bs.run_command(upload_cmd)
# send auto-slack message about it? (nah, too much)
access_cmd = f"/home/al/code/google-cloud-sdk/bin/gsutil acl ch -R -u AllUsers:R gs://andersen-lab_temp/outbreak_genomics/*"
bs.run_command(access_cmd)
# GENERATE CSV DATA MODEL
start = time.time()
muts.to_csv(out_fp, index=False)
end = time.time()
print(f'Execution time: {end - start} seconds')