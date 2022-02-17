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
from mappings import COUNTY_CORRECTIONS
import numpy as np
from rapidfuzz import process, fuzz

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
parser.add_argument("-p", "--inputpango",
                        type=str,
                        required=True,
                        help="Input pango assignments")
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
input_pango = args.inputpango
out_fp = args.outfp
unknown_val = args.unknownvalue
min_date = args.mindate
geojson_prefix = args.geojson
current_datetime = args.currentdate
date_modified = '-'.join(current_datetime.split('-')[:3]) + '-' + ':'.join(current_datetime.split('-')[3:])
max_date = '-'.join(current_datetime.split('-')[:3])

try:
    meta = pd.read_csv(input_metadata, sep='\t')
    pango = pd.read_csv(input_pango)
    muts = pd.read_csv(input_mut, dtype=str)
    if len(meta.index) == 0 or len(muts.index) == 0:
        raise "empty dataframe"
except Exception as e:
    pd.DataFrame().to_json(out_fp, orient='records', lines=True)
    sys.exit(0)

pango = pango[["taxon", "lineage"]].rename(columns={"taxon": "strain", "lineage": "pangolin_lineage"})
meta = meta.drop(columns=["pangolin_lineage"])
meta = pd.merge(meta, pango, on="strain", how="left")
muts = muts[~(muts['gene'].isin(['5UTR', '3UTR']))]
# ignore mutations found in non-coding regions
muts = muts.loc[~(muts['gene']=='Non-coding region')]
# fuse with metadata
print(f"Fusing muts with metadata...")
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
data = pd.merge(meta, muts, on='strain', how='left')

# normalize date information
data['tmp'] = data['date_collected'].str.split('-')
data = data[data['tmp'].str.len()>=2]
data.loc[data['tmp'].str.len()==2, 'date_collected'] += '-15'
data['date_collected'] = pd.to_datetime(data['date_collected'], errors='coerce')
data['date_collected'] = data['date_collected'].astype(str)
data = data[data['date_collected'] <= max_date]
data = data[data['date_collected'] > min_date]
data['date_modified'] = date_modified

# TODO: handle locstring off-by-one errors

data['locstring'] = data['locstring'].str.lower().str.replace('\.', '').str.replace('unknown', '').str.split("/")
data['country'] = data['locstring'].apply(lambda x: x[1] if len(x) >= 2 else '').str.strip()
data['division'] = data['locstring'].apply(lambda x: x[2] if len(x) >= 3 else '').str.strip()
data['location'] = data['locstring'].apply(lambda x: x[3] if len(x) >= 4 else '').str.strip()
data['country'].fillna('', inplace=True)
data.loc[:, 'country'] = data['country'].str.strip()
data.loc[data['country'].str.len() <= 1, 'country'] = ''
data.loc[data['division'].isna(), 'division'] = ''
data['division'] = data['division'].copy()
data.loc[:, 'division'] = data['division'].str.strip()
data.loc[data['division'].str.len() <= 1, 'division'] = ''
data.loc[data['location'].isna(), 'location'] = ''

# TODO: apply NextStrain replacements file

# Patches for known causes of geo errors
data.loc[data['country'] == 'usa', 'country'] = 'united states'
def norm_territories(x):
    if x.country in ['puerto rico', 'guam', 'us virgin islands', 'northern mariana islands', 'american samoa']:
        (x.country, x.division, x.location) = ('united states', x.country, x.division)
    return x
data = data.apply(norm_territories, axis=1)
data.loc[(data['country'].str.contains('congo')) & (data['country'].str.contains(
    'democratic')), 'country'] = 'democratic republic of the congo'
data.loc[(data['country'].str.contains('congo')) & ~(data['country'].str.contains(
    'democratic')), 'country'] = 'republic of congo'
data.loc[data['country'].str.contains(
    'eswatini'), 'country'] = "swaziland"
data.loc[data['country'].str.contains(
    'bonaire'), 'country'] = "bonaire, sint eustatius and saba"
data.loc[data['country'].str.contains(
    'sint eustatius'), 'country'] = "bonaire, sint eustatius and saba"
data.loc[data['country']=="jonavos apskritis", 'country'] = "lithuania"
data.loc[data['division'].str.contains('state of mexico'), 'division'] = 'méxico'
data.loc[data['division'].str.contains('bethlehem'), 'division'] = 'west bank'
data.loc[data['division'].str.contains('hebron'), 'division'] = 'west bank'
data.loc[data['division'].str.contains('jenin'), 'division'] = 'west bank'
data.loc[data['division'].str.contains('jericho'), 'division'] = 'west bank'
data.loc[data['division'].str.contains('ramallah'), 'division'] = 'west bank'
data.loc[data['division'].str.contains('tulkarem'), 'division'] = 'west bank'
data.loc[data['division'].str.contains('nablus'), 'division'] = 'west bank'
data.loc[data['division'].str.contains('copenhagen'), 'division'] = 'hovedstaden'
data.loc[data['division'].str.contains('sjaelland'), 'division'] = 'sjælland'
data.loc[data['division'].str.contains('cape town'), 'division'] = 'western cape'
data.loc[data['division'].str.contains('western cape'), 'division'] = 'western cape'
data.loc[data['division'].str.contains('asturias'), 'division'] = 'principado de asturias'
data.loc[data['division'].str.contains('balear_islands'), 'division'] = 'islas baleamuts'
data.loc[data['division'].str.contains('illes balears'), 'division'] = 'islas baleamuts'
data.loc[data['division'].str.contains('canary islands'), 'division'] = 'canaries'
data.loc[data['division'].str.contains('ceuta'), 'division'] = 'ceuta y melilla'
data.loc[data['division'].str.contains('melilla'), 'division'] = 'ceuta y melilla'
data.loc[data['division'].str.contains('bavaria'), 'division'] = 'bayern'
data.loc[data['division'].str.contains('lower saxony'), 'division'] = 'niedersachsen'
data.loc[data['division'].str.contains('mecklenburg-western pomerania'), 'division'] = 'mecklenburg-vorpommern'
data.loc[data['division'].str.contains('rhineland-palatinate'), 'division'] = 'rheinland-pfalz'
data.loc[(data['division'].str.contains('saxony-anhalt'))
        & (data['country'].str.contains('germany')), 'division'] = 'sachsen-anhalt'
data.loc[(data['division'].str.contains('saxony'))
        & (data['country'].str.contains('germany')), 'division'] = 'sachsen'
data.loc[data['division'].str.contains('north rhine-westphalia'), 'division'] = 'nordrhein-westfalen'
data.loc[data['division'].str.contains('argovie'), 'division'] = 'aargau'
data.loc[data['division'].str.contains('geneva'), 'division'] = 'geneve'
data.loc[data['division'].str.contains('grabunden'), 'division'] = 'graubunden'
data.loc[data['division'].str.contains('luzern'), 'division'] = 'lucerne'
data.loc[data['division'].str.contains('neuenburg'), 'division'] = 'neuchatel'
data.loc[data['division'].str.contains('obwald'), 'division'] = 'obwalden'
data.loc[data['division'].str.contains('saint-gall'), 'division'] = 'sankt gallen'
data.loc[data['division'].str.contains('gallen'), 'division'] = 'sankt gallen'
data.loc[data['division'].str.contains('turgovia'), 'division'] = 'thurgau'
data.loc[data['division'].str.contains('waadt'), 'division'] = 'vaud'
data.loc[data['division'].str.contains('wallis'), 'division'] = 'valais'
data.loc[:, 'location'] = data['location'].str.replace('county', '').str.replace('parish', '').str.replace(',', '')
data.loc[data['location'].str.len() <= 1, 'location'] = ''
data.loc[data['location']=='unk', 'location'] = ''
data.loc[data['division']==data['country'], 'division'] = ''
data.fillna('', inplace=True)

# Match to GADM-derived loc/id database
with open('{}/gadm_transformed.jsonl'.format(geojson_prefix)) as f:
    countries = [json.loads(line) for line in f]
get_geo = lambda b, y: (b[y[2]-1] if 'alias' in b[y[2]] else b[y[2]]) if not y is None else None
s = lambda f, r: lambda a, b, **params: 100 if a == '' and b == '' else (r if a == '' or b == '' else f(a, b, **params))
data.loc[:, 'country_match'] = data['country'].apply(lambda x: get_geo(countries, process.extractOne(x, [country['name'] for country in countries], scorer=s(fuzz.ratio, 80))))
data = data.dropna(axis=0, subset=['country_match'])
data.loc[:, 'division_match'] = data.apply(lambda x: get_geo(x.country_match['sub'], process.extractOne(x.division, [division['name'] for division in x.country_match['sub']], scorer=s(fuzz.ratio, 70))), axis=1)
data.loc[:, 'location_match'] = data.apply(lambda x: get_geo(x.division_match['sub'], process.extractOne(x.location, [location['name'] for location in x.division_match['sub']], scorer=s(fuzz.partial_ratio, 60))), axis=1)
data.loc[:, 'country'] = data['country_match'].apply(lambda x: x['name'])
data.loc[:, 'country_id'] = data['country_match'].apply(lambda x: x['id'])
data.loc[:, 'division'] = data['division_match'].apply(lambda x: x['name'])
data.loc[:, 'division_id'] = data['division_match'].apply(lambda x: x['id'])
data.loc[:, 'location'] = data['location_match'].apply(lambda x: x['name'])
data.loc[:, 'location_id'] = data['location_match'].apply(lambda x: x['id'])
data.loc[data['country_id'] == '', 'country_id'] = unknown_val
data.loc[data['location_id'] == '', 'location_id'] = unknown_val
data.loc[data['division_id'] == '', 'division_id'] = unknown_val
data = data.drop(['country_match', 'division_match', 'location_match'], axis=1)

# TODO: asciify these
data['country_lower'] = data['country'].str.lower()
data['division_lower'] = data['division'].str.lower()
data['location_lower'] = data['location'].str.lower()

data.to_json(out_fp, orient='records', lines=True)
