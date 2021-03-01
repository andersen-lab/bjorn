import sys
import time
import json
import numpy as np
import pandas as pd

import bjorn_support as bs
import mutations as bm
import visualize as bv
import reports as br
import data as bd


date = '2021-02-16'
unknown_val = 'None'


countries_fp = '/home/al/data/geojsons/countries.geo.json'
states_fp = '/home/al/data/geojsons/us-states.json'
subs = pd.read_csv(f'/valhalla/gisaid/subs_long_{date}.csv.gz', 
                   compression='gzip')
dels = pd.read_csv(f'/valhalla/gisaid/dels_long_{date}.csv.gz', 
                   compression='gzip')


cols = ['mutation', 'strain', 'country', 'division', 'location', 'date', 'absolute_coords', 'del_len']

dels['pos'] = dels['absolute_coords'].apply(lambda x: int(x.split(':')[0]))
dels['ref_codon'] = dels['del_seq'].copy()

dels['gene_start_pos'] = dels['gene'].apply(lambda x: bd.GENE2POS[x]['start']+2)

dels['pos_in_codon'] = (dels['pos'] - dels['gene_start_pos']) % 3

print('Relative positions in codon counts [internal]')
print(dels['pos_in_codon'].value_counts())

dels['mutation'] = dels[['pos_in_codon', 'gene', 'codon_num', 'del_len']].apply(bm.assign_deletion, axis=1)
dels['deletion_codon_coords'] = dels[['pos_in_codon', 'gene', 'codon_num', 'del_len']].apply(bm.assign_deletion_codon_coords, axis=1)

dels['is_frameshift'] = dels['del_len'].apply(is_frameshift)

print(f'Total substitution count: {subs.shape[0]}')
print(f'Total deletion count: {dels.shape[0]}')
subs['type'] = 'substitution'
muts = pd.concat([subs, dels])
print(f'Total mutation count: {muts.shape}')
# date cleaning
muts['tmp'] = muts['date'].str.split('-')
muts = muts[muts['tmp'].str.len()>=2]
muts.loc[muts['tmp'].str.len()==2, 'date'] += '-15'
muts['date'] = pd.to_datetime(muts['date'], errors='coerce')
muts = muts[muts['date']<date]
muts['date'] = muts['date'].astype(str)
test = muts.groupby('date').agg(num_samples=('strain', 'nunique')).reset_index()
assert test[test.duplicated(subset=['date'])].shape[0]==0
with open(countries_fp) as f:
    countries = json.load(f)
country_map = {x['properties']['name']: x['id'] for x in countries['features']}
# country fips id
muts['country_id'] = muts['country'].apply(lambda x: country_map.get(x, unknown_val))
with open(states_fp) as f:
    states = json.load(f)
state_map = {x['properties']['name']: x['id'] for x in states['features']}
# state fips id
muts['division_id'] = muts['division'].apply(lambda x: state_map.get(x, unknown_val))
# column renaming
muts.rename(columns={
    'date': 'date_collected',
    'GISAID_clade': 'gisaid_clade',
    'Nextstrain_clade': 'nextstrain_clade',
    'del_len': 'change_length_nt'
    }, inplace=True)
muts['date_modified'] = date
muts['is_synonymous'] = False
muts.loc[muts['ref_aa']==muts['alt_aa'], 'is_synonymous'] = True
muts.loc[muts['location']=='unk', 'location'] = unknown_val
muts.loc[muts['purpose_of_sequencing']=='?', 'purpose_of_sequencing'] = unknown_val
muts.loc[muts['genbank_accession']=='?', 'genbank_accession'] = unknown_val
muts.fillna(unknown_val, inplace=True)
muts.loc[muts['division']==muts['country'], 'division'] = unknown_val
sample_ids = muts[['strain']].drop_duplicates().sample(10)['strain'].unique()
test = muts[muts['strain'].isin(sample_ids)]
meta_info = ['strain', 'date_modified',
        'date_collected','date_submitted',
        'country_id', 'country', 
        'division_id', 'division', 'location', 
        'submitting_lab', 'originating_lab',
        'authors', 'pangolin_lineage', 
        'gisaid_clade', 'nextstrain_clade',
        'gisaid_epi_isl', 'genbank_accession',
        'purpose_of_sequencing']

muts_info = ['type', 'mutation', 'gene', 
             'ref_codon', 'pos', 'alt_codon', 
             'is_synonymous', 
             'ref_aa', 'codon_num', 'alt_aa', 
             'absolute_coords', 
             'change_length_nt', 'is_frameshift',
             'deletion_codon_coords']
# TEST 
start = time.time()
test_fp = f'data/TEST_data_model_{date}.json'
(test.groupby(meta_info, as_index=True)
             .apply(lambda x: x[muts_info].to_dict('records'))
             .reset_index()
             .rename(columns={0:'mutations'})
             .to_json(test_fp, 
                      orient='records',
#                       compression='gzip'
                     )
)
end = time.time()
print(f'Test Execution time: {end - start} seconds')
print(f"Test data generated and saved in {test_fp}")
print(f"Latest sample submission: {muts['date_submitted'].max()}")
print(f'Latest sample collection: {muts['date'].max()}')
print(f'Most prominent deletions...')
print(muts[muts['type']=='deletion']['mutation'].value_counts())
print(f'Most prominent substitutions...')
print(muts[muts['type']=='substitution']['mutation'].value_counts())
assert muts.loc[(muts['type']=='deletion')&(muts['mutation']==1)][['pos', 'change_length_nt', 'mutation']].shape[0]==0
# GENERATE JSON DATA MODEL
print(f"Generating REAL data...")
start = time.time()
data_fp = 'data/data.json.gz'
(muts.groupby(meta_info, as_index=True)
             .apply(lambda x: x[muts_info].to_dict('records'))
             .reset_index()
             .rename(columns={0:'mutations'})
             .to_json(data_fp, 
                      orient='records',
                      compression='gzip'))
end = time.time()
print(f'REAL Execution time: {end - start} seconds')
print(f'REAL data generated and saved in {data_fp}')