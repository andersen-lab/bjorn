#!/usr/bin/env python
import argparse
import pandas as pd
import json
import re
from rapidfuzz import process, fuzz

# COLLECTING USER PARAMETERS
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputdata",
                        type=str,
                        required=True,
                        help="Input data")
parser.add_argument("-o", "--outfp",
                        type=str,
                        required=True,
                        help="Output filepath")
parser.add_argument("-u", "--unknownvalue",
                        type=str,
                        required=True,
                        help="Unknown value")
parser.add_argument("-g", "--geojson",
                        type=str,
                        required=True,
                        help="GeoJSON prefix")

args = parser.parse_args()
data = pd.read_csv(args.inputdata, sep='\t', names=['accession_id', 'date_collected', 'date_submitted', 'locstring', 'sequence', 'pangolin_lineage', 'scorpio_constellation', 'full_lineage', 'mutations'])
data = data.drop(['sequence', 'scorpio_constellation', 'full_lineage'], axis=1)
data['accession_id'] = data['accession_id'].astype(str)

# normalize date information
data['tmp'] = data['date_collected'].str.split('-')
data = data[data['tmp'].str.len()>=2]
data.loc[data['tmp'].str.len()==2, 'date_collected'] += '-15'
data['date_collected'] = pd.to_datetime(data['date_collected'], errors='coerce')
data['date_collected'] = data['date_collected'].astype(str)
data = data.drop(['tmp'], axis=1)

# expand gofasta variants mutation string into rich json
# split to minimal mutation descriptions
data['mutations'].fillna('', inplace=True)
data['mutations'] = data['mutations'].str.replace(')','').str.split('|')
def parsemutation(mut):
    mut, *nucs = mut.split('(')
    mut = mut.split(':')
    kind = {'aa':'substitution', 'nuc':'synonymous', 'ins': 'insertion', 'del': 'deletion'}[mut.pop(0)]
    if kind == 'substitution':
        gene = mut.pop(0)
        ref_aa, codon_num, alt_aa = re.split('(?<=\d)(?=\D)|(?<=\D)(?=\d)', mut[0])
    else:
        gene, alt_aa, codon_num, ref_aa = 'None', 'N', 0, 'N'
    if kind == 'synonymous':
        nucs.insert(0, 'nuc:' + mut[0])
    nucs = [nuc for n in nucs for nuc in n.split(';')]
    if len(nucs) == 0: nucs = ['nuc:N0N']
    def parseNuc(nuc):
        nuc = re.split('(?<=\d)(?=\D)|(?<=\D)(?=\d)', nuc.split(':')[1])
        return {'ref_base': nuc[0], 'pos': int(nuc[1]), 'alt_base': nuc[2]}
    return [{ **{
        'is_synonymous': (kind == 'synonymous'),
        'type': kind,
        'mutation': mut[0],
        'gene': gene,
        'codon_num': codon_num,
        'ref_aa': ref_aa,
        'alt_aa': alt_aa }, **parseNuc(nuc)} for nuc in nucs]
data['mutations'] = data['mutations'].map(lambda muts: [m for mut in muts for m in parsemutation(mut)] if isinstance(muts, list) else [])


# data['mutations'] = data['mutations'].str.replace(')','').str.replace('aa:', '').str.split('|')
# def parsemutation(mut):
#     mut = list(reversed([subpart.split(':') for part in mut.split('(') for subpart in part.split(';')]))
#     for part in mut: part[-1] = re.split('(?<=\d)(?=\D)|(?<=\D)(?=\d)', part[-1])
#     out = {'is_synonymous': len(mut) <= 1}
#     if mut[0][0] == 'nuc':
#         out.update({
#             'type': 'substitution',
#             'pos': mut[0][1][1],
#             'ref_base': mut[0][1][0],
#             'alt_base': mut[0][1][2]
#         })
#     return out
# data['mutations'] = data['mutations'].map(lambda muts: [parsemutation(mut) for mut in muts] if isinstance(muts, list) else [])



# below is all geo-handling

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
data = data.drop(['locstring'], axis=1)

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
with open(args.geojson) as f:
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
data.loc[data['country_id'] == '', 'country_id'] = args.unknownvalue
data.loc[data['location_id'] == '', 'location_id'] = args.unknownvalue
data.loc[data['division_id'] == '', 'division_id'] = args.unknownvalue
data = data.drop(['country_match', 'division_match', 'location_match'], axis=1)

# TODO: asciify these
data['country_lower'] = data['country'].str.lower()
data['division_lower'] = data['division'].str.lower()
data['location_lower'] = data['location'].str.lower()

# the export
data.to_json(args.outfp, orient='records', lines=True)
