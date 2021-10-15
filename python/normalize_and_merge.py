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

#make an empy list for the data to be held in
std_loc = [[0]*3]*len(meta) 
for index, row in meta.iterrows():
    locstring = row['locstring'] 
    loc_list = locstring.split("/")

    #we have a continent in the string
    if len(loc_list) == 4:
        loc_list = loc_list[1:]    
     
    #we have the country and either division or location
    if len(loc_list) == 2:    
        sys.exit(0)
    #we have country, division and location
    if len(loc_list) == 3:
        for i,loc in enumerate(loc_list):
            high_ratio = 0
            high_string = ''
            compare_dict = std_locs[i]
            for key in compare_dict.keys():
                #for divisions we may have countries attached
                if i == 1:
                    country = key.split('-')[0]
                    if country.lower() != std_loc[index][0].lower():
                        continue
                #for locations we may have countries attached
                if i == 2:
                    country = key.split('-')[0]
                    if country.lower() != std_loc[index][0].lower():
                        continue
                ratio = fuzz.ratio(key.lower(),loc.lower())
                if ratio > high_ratio:
                    high_ratio = ratio
                    high_string = key
            for value in compare_dict.values():
                ratio = fuzz.ratio(value.lower(),loc.lower())
                if ratio > high_ratio:
                    high_ratio = ratio
                    high_string = list(compare_dict.keys())[list(compare_dict.values()).index(value)]
            std_loc[index][i] = high_string
    print(loc_list, std_loc[index])        
loc_df = pd.DataFrame(std_loc, columns=['country','division','location'])
meta = pd.concat([meta,loc_df], axis=1)
print(meta)

print(f'admin1 standardization...')
meta.loc[meta['division'].isna(), 'division'] = unknown_val+longstring
meta['division'] = meta['division'].copy()
meta.loc[:, 'division'] = meta['division'].str.strip()
meta.loc[meta['division'].str.len() <= 1, 'division'] = unknown_val+longstring
meta.loc[meta['division']=='usa', 'division'] = 'united states'
meta.loc[meta['division'].str.contains('georgia /'), 'division'] = 'georgia'
meta.loc[meta['division'].str.contains('antwerp'), 'division'] = 'vlaanderen'
meta.loc[meta['division'].str.contains('andalu'), 'division'] = 'andalucía'
meta.loc[meta['division'].str.contains('cairo'), 'division'] = 'al qahirah'
meta.loc[meta['division'].str.contains('northern territory'), 'division'] = 'northern territory'
meta.loc[meta['division'].str.contains('fayoum'), 'division'] = 'al fayyum'
meta.loc[meta['division'].str.contains('musca'), 'division'] = 'muscat'
meta.loc[meta['division'].str.contains('kalyoubia'), 'division'] = 'al qalyubiyah'
meta.loc[meta['division'].str.contains('buraymi'), 'division'] = 'al buraymi'
meta.loc[meta['division'].str.contains('buraimi'), 'division'] = 'al buraymi'
meta.loc[meta['division'].str.contains('dakhiliyah'), 'division'] = 'ad dakhliyah'
meta.loc[meta['division'].str.contains('dhahirah'), 'division'] = 'al dhahira'
meta.loc[meta['division'].str.contains('north batinah'), 'division'] = 'al batinah north'
meta.loc[meta['division'].str.contains('south batinah'), 'division'] = 'al batinah south'
meta.loc[meta['division'].str.contains('north sharqiyah'), 'division'] = 'ash sharqiyah north'
meta.loc[meta['division'].str.contains('wuhan'), 'division'] = 'hubei'
meta.loc[meta['division'].str.contains('quebec'), 'division'] = 'québec'
meta.loc[meta['division'].str.contains('toronto'), 'division'] = 'ontario'
meta.loc[meta['division'].str.contains('coahuila de zaragoza'), 'division'] = 'coahuila'
meta.loc[meta['division'].str.contains('mexico city'), 'division'] = 'méxico'
meta.loc[meta['division'].str.contains('michoacan'), 'division'] = 'michoacán'
meta.loc[meta['division'].str.contains('nuevo leon'), 'division'] = 'nuevo león'
meta.loc[meta['division'].str.contains('queretaro'), 'division'] = 'querétaro'
meta.loc[meta['division'].str.contains('sanluispotosi'), 'division'] = 'san luis potosí'
meta.loc[meta['division'].str.contains('san luis potosi'), 'division'] = 'san luis potosí'
meta.loc[meta['division'].str.contains('state of mexico'), 'division'] = 'méxico'
meta.loc[meta['division'].str.contains('yucatan'), 'division'] = 'yucatán'
meta.loc[meta['division'].str.contains('bethlehem'), 'division'] = 'west bank'
meta.loc[meta['division'].str.contains('hebron'), 'division'] = 'west bank'
meta.loc[meta['division'].str.contains('jenin'), 'division'] = 'west bank'
meta.loc[meta['division'].str.contains('jericho'), 'division'] = 'west bank'
meta.loc[meta['division'].str.contains('ramallah'), 'division'] = 'west bank'
meta.loc[meta['division'].str.contains('tulkarem'), 'division'] = 'west bank'
meta.loc[meta['division'].str.contains('nablus'), 'division'] = 'west bank'
meta.loc[meta['division'].str.contains('sharja'), 'division'] = 'sharjah'
meta.loc[meta['division'].str.contains('copenhagen'), 'division'] = 'hovedstaden'
meta.loc[meta['division'].str.contains('sjaelland'), 'division'] = 'sjælland'
meta.loc[meta['division'].str.contains('cape town'), 'division'] = 'western cape'
meta.loc[meta['division'].str.contains('western cape'), 'division'] = 'western cape'
meta.loc[meta['division'].str.contains('amapa'), 'division'] = 'amapá'
meta.loc[meta['division'].str.contains('ceara'), 'division'] = 'ceará'
meta.loc[meta['division'].str.contains('goias'), 'division'] = 'goiás'
meta.loc[meta['division'].str.contains('maranhao'), 'division'] = 'maranhão'
meta.loc[meta['division'].str.contains('paraiba'), 'division'] = 'paraíba'
meta.loc[meta['division'].str.contains('parana'), 'division'] = 'paraná'
meta.loc[meta['division'].str.contains('piaui'), 'division'] = 'piauí'
meta.loc[meta['division'].str.contains('sao paulo'), 'division'] = 'são paulo'
meta.loc[meta['division'].str.contains('aragon'), 'division'] = 'aragón'
meta.loc[meta['division'].str.contains('asturias'), 'division'] = 'principado de asturias'
meta.loc[meta['division'].str.contains('balear islands'), 'division'] = 'islas baleamuts'
meta.loc[meta['division'].str.contains('balear_islands'), 'division'] = 'islas baleamuts'
meta.loc[meta['division'].str.contains('illes balears'), 'division'] = 'islas baleamuts'
meta.loc[meta['division'].str.contains('canary islands'), 'division'] = 'canaries'
meta.loc[meta['division'].str.contains('canary_islands'), 'division'] = 'canaries'
meta.loc[meta['division'].str.contains('castilla la mancha'), 'division'] = 'castilla-la mancha'
meta.loc[meta['division'].str.contains('castilla la mancha'), 'division'] = 'castilla-la mancha'
meta.loc[meta['division'].str.contains('castilla y leon'), 'division'] = 'castilla y león'
meta.loc[meta['division'].str.contains('ceuta'), 'division'] = 'ceuta y melilla'
meta.loc[meta['division'].str.contains('melilla'), 'division'] = 'ceuta y melilla'
meta.loc[meta['division'].str.contains('comunitat valenciana'), 'division'] = 'comunidad valenciana'
meta.loc[meta['division'].str.contains('comunitat_valenciana'), 'division'] = 'comunidad valenciana'
meta.loc[meta['division'].str.contains('la_rioja'), 'division'] = 'la rioja'
meta.loc[meta['division'].str.contains('madrid'), 'division'] = 'comunidad de madrid'
meta.loc[meta['division'].str.contains('murcia'), 'division'] = 'región de murcia'
meta.loc[meta['division'].str.contains('navarra'), 'division'] = 'comunidad foral de navarra'
meta.loc[meta['division'].str.contains('catalunya'), 'division'] = 'cataluña'
meta.loc[meta['division'].str.contains('catalonia'), 'division'] = 'cataluña'
# germany
meta.loc[meta['division'].str.contains('baden-wuerttemberg'), 'division'] = 'baden-württemberg'
meta.loc[meta['division'].str.contains('baden-wurttemberg'), 'division'] = 'baden-württemberg'
meta.loc[meta['division'].str.contains('bavaria'), 'division'] = 'bayern'
meta.loc[meta['division'].str.contains('hesse'), 'division'] = 'hessen'
meta.loc[meta['division'].str.contains('lower saxony'), 'division'] = 'niedersachsen'
meta.loc[meta['division'].str.contains('mecklenburg-western pomerania'), 'division'] = 'mecklenburg-vorpommern'
meta.loc[meta['division'].str.contains('rhineland-palatinate'), 'division'] = 'rheinland-pfalz'
meta.loc[(meta['division'].str.contains('saxony-anhalt'))
        & (meta['country'].str.contains('germany')), 'division'] = 'sachsen-anhalt'
meta.loc[(meta['division'].str.contains('saxony'))
        & (meta['country'].str.contains('germany')), 'division'] = 'sachsen'
meta.loc[meta['division'].str.contains('north rhine-westphalia'), 'division'] = 'nordrhein-westfalen'
meta.loc[meta['division'].str.contains('thuringia'), 'division'] = 'thüringen'
# south africa
meta.loc[(meta['country'].str.contains('south africa'))
        & (meta['division'].str.contains('kwazulu natal')),
        'division'] = 'kwazulu-natal'
meta.loc[(meta['country'].str.contains('south africa'))
        & (meta['division'].str.contains('northern cape province')),
        'division'] = 'northern cape'
# austria
meta.loc[(meta['country']=='austria')
        & (meta['division']=='tyrol'), 'division'] = 'tirol'
# switzerland
meta.loc[meta['division'].str.contains('argovie'), 'division'] = 'aargau'
meta.loc[meta['division'].str.contains('geneva'), 'division'] = 'genève'
meta.loc[meta['division'].str.contains('grabunden'), 'division'] = 'graubünden'
meta.loc[meta['division'].str.contains('luzern'), 'division'] = 'lucerne'
meta.loc[meta['division'].str.contains('neuchatel'), 'division'] = 'neuchâtel'
meta.loc[meta['division'].str.contains('neuenburg'), 'division'] = 'neuchâtel'
meta.loc[meta['division'].str.contains('obwald'), 'division'] = 'obwalden'
meta.loc[meta['division'].str.contains('saint-gall'), 'division'] = 'sankt gallen'
meta.loc[meta['division'].str.contains('st gallen'), 'division'] = 'sankt gallen'
meta.loc[meta['division'].str.contains('st. gallen'), 'division'] = 'sankt gallen'
meta.loc[meta['division'].str.contains('schaffhouse'), 'division'] = 'schaffhausen'
meta.loc[meta['division'].str.contains('turgovia'), 'division'] = 'thurgau'
meta.loc[meta['division'].str.contains('valais'), 'division'] = 'valais'
meta.loc[meta['division'].str.contains('waadt'), 'division'] = 'vaud'
meta.loc[meta['division'].str.contains('wallis'), 'division'] = 'valais'
meta.loc[meta['division'].str.contains('zaerich'), 'division'] = 'zürich'
meta.loc[meta['division'].str.contains('zoerich'), 'division'] = 'zürich'
meta.loc[meta['division'].str.contains('zurich'), 'division'] = 'zürich'
meta.loc[meta['division'].str.contains('unkown'), 'division'] = 'unknown'
meta.loc[meta['division'].str.contains('basel-land'), 'division'] = 'basel-landschaft'
meta.loc[meta['division'].str.contains('grisons'), 'division'] = 'graubünden'
meta.loc[meta['division'].str.contains('st.gallen'), 'division'] = 'sankt gallen'
meta.loc[meta['division'].str.contains('geneve'), 'division'] = 'genève'
meta.loc[meta['division'].str.contains('stankt gallen'), 'division'] = 'sankt gallen'
# sweden
meta.loc[(meta['country']=='sweden')
        &(meta['division']=='sodermanland'), 'division'] = 'södermanland'
meta.loc[(meta['country']=='sweden')
        &(meta['division']=='gavleborg'), 'division'] = 'gävleborg'
meta.loc[(meta['country']=='sweden')
        &(meta['division']=='jamtland'), 'division'] = 'jämtland'
meta.loc[(meta['country']=='sweden')
        &(meta['division']=='jonkoping'), 'division'] = 'jönköping'
meta.loc[(meta['country']=='sweden')
        &(meta['division']=='ostergotland'), 'division'] = 'östergötland'
meta.loc[(meta['country']=='sweden')
        &(meta['division']=='skane'), 'division'] = 'skåne'
meta.loc[(meta['country']=='sweden')
        &(meta['division']=='varmland'), 'division'] = 'värmland'
meta.loc[(meta['country']=='sweden')
        &(meta['division']=='vasterbotten'), 'division'] = 'västerbotten'
meta.loc[(meta['country']=='sweden')
        &(meta['division']=='vasternorrland'), 'division'] = 'västernorrland'
meta.loc[(meta['country']=='sweden')
        &(meta['division']=='vastmanland'), 'division'] = 'västmanland'
meta.loc[(meta['country']=='sweden')
        &(meta['division']=='vastra gotaland'), 'division'] = 'västra götaland'
# puerto rico
meta.loc[(meta['division']=='puerto rico'), 'country'] = 'puerto rico'
meta.loc[(meta['country']=='puerto rico'), 'division'] = meta.loc[(meta['country']=='puerto rico'), 'location'].fillna(unknown_val+longstring)

print(f'admin2 standardization')
meta.loc[meta['location'].isna(), 'location'] = unknown_val+longstring
houston_filter = (meta['division']=='texas') & (meta['location']=='houston')
meta.loc[houston_filter, 'location'] = 'harris county'
for key, val in COUNTY_CORRECTIONS.items():
    meta.loc[:, 'location'] = meta['location'].str.replace(key, val)
meta.loc[:, 'location'] = meta['location'].str.replace('county', '').str.replace('county', '').str.replace(',', '')
meta.loc[:, 'location'] = meta['location'].str.strip().apply(bs.check_state, args=(False,)).str.strip()
meta.loc[meta['location'].str.len() <= 1, 'location'] = unknown_val+longstring
meta.loc[meta['location'].str.contains('anchorage-mat-su'), 'location'] = 'anchorage'
meta.loc[meta['location'].str.contains('anchorage-mat su'), 'location'] = 'anchorage'
meta.loc[meta['location'].str.contains('bra'), 'location'] = 'brazos'
#meta.loc[meta['location'].str.contains('br'), 'location'] = 'brewster'
meta.loc[meta['location'].str.contains('belgrade'), 'location'] = 'gallatin'
meta.loc[meta['location'].str.contains('bozeman'), 'location'] = 'gallatin'
meta.loc[meta['location'].str.contains('big sky'), 'location'] = 'gallatin'
meta.loc[meta['location'].str.contains('belton'), 'location'] = 'bell'
meta.loc[meta['location'].str.contains('brentwood'), 'location'] = 'contra costa'
meta.loc[meta['location'].str.contains('chicago'), 'location'] = 'cook'
meta.loc[meta['location'].str.contains('colombus'), 'location'] = 'franklin'
meta.loc[meta['location'].str.contains('dubois'), 'location'] = 'fremont'
meta.loc[meta['location'].str.contains('dupage'), 'location'] = 'dupage'
meta.loc[meta['location'].str.contains('eau claire'), 'location'] = 'eau claire'
meta.loc[meta['location'].str.contains('ennis'), 'location'] = 'ellis'
meta.loc[meta['location'].str.contains('fond du lac'), 'location'] = 'fond du lac'
meta.loc[meta['location'].str.contains('fond du lac'), 'location'] = 'fond du lac'
meta.loc[meta['location'].str.contains('fonddu lac'), 'location'] = 'fond du lac'
meta.loc[meta['location'].str.contains('frisco'), 'location'] = 'collin'
meta.loc[meta['location'].str.contains('hawai'), 'location'] = 'hawaii'
meta.loc[meta['location'].str.contains('holland'), 'location'] = 'ottawa'
meta.loc[meta['location'].str.contains('honolul'), 'location'] = 'honolulu'
meta.loc[meta['location'].str.contains('indianapolis'), 'location'] = 'marion'
meta.loc[meta['location'].str.contains('interior'), 'location'] = 'fairbanks north star'
meta.loc[meta['location'].str.contains('ithaca'), 'location'] = 'tompkins'
meta.loc[meta['location'].str.contains('kaua'), 'location'] = 'kauai'
meta.loc[meta['location'].str.contains('las vegas'), 'location'] = 'clark'
meta.loc[meta['location'].str.contains('mau'), 'location'] = 'hawaii'
meta.loc[meta['location'].str.contains('mcculloch'), 'location'] = 'mcculloch'
meta.loc[meta['location'].str.contains('mchenry'), 'location'] = 'mchenry'
meta.loc[meta['location'].str.contains('mclennan'), 'location'] = 'mclennan'
meta.loc[meta['location'].str.contains('moris'), 'location'] = 'morris'
meta.loc[meta['location'].str.contains('new york'), 'location'] = 'new york'
meta.loc[meta['location'].str.contains('new york city'), 'location'] = 'new york'
meta.loc[meta['location'].str.contains('new hyde park'), 'location'] = 'nassau'
meta.loc[meta['location'].str.contains('new orleans'), 'location'] = 'orleans'
meta.loc[meta['location'].str.contains('new rochelle'), 'location'] = 'westchester'
meta.loc[meta['location'].str.contains('northern'), 'location'] = 'fairbanks north star'
meta.loc[meta['location'].str.contains('omaha'), 'location'] = 'douglas'
meta.loc[meta['location'].str.contains('ostego'), 'location'] = 'allegan'
meta.loc[meta['location'].str.contains('phoenix'), 'location'] = 'maricopa'
meta.loc[meta['location'].str.contains('san bernadino'), 'location'] = 'san bernardino'
meta.loc[meta['location'].str.contains('seattle'), 'location'] = 'king'
meta.loc[meta['location'].str.contains('st. bernard'), 'location'] = 'saint bernard'
meta.loc[meta['location'].str.contains('st. clair'), 'location'] = 'saint clair'
meta.loc[meta['location'].str.contains('st. lawrence'), 'location'] = 'saint lawrence'
meta.loc[meta['location'].str.contains('st. louis'), 'location'] = 'saint louis'
meta.loc[meta['location'].str.contains('st. tammany'), 'location'] = 'saint tammany'
meta.loc[meta['location'].str.contains('staten island'), 'location'] = 'richmond'
meta.loc[meta['location'].str.contains('thurson'), 'location'] = 'thurston'
meta.loc[meta['location'].str.contains('tucson'), 'location'] = 'pima'
meta.loc[meta['location'].str.contains('west yellowstone'), 'location'] = 'gallatin'
meta.loc[meta['location'].str.contains('adam'), 'location'] = 'adams'
meta.loc[meta['location'].str.contains('alachu'), 'location'] = 'alachua'
meta.loc[meta['location'].str.contains('du bois'), 'location'] = 'dubois'
meta.loc[meta['location'].str.contains('desoto'), 'location'] = 'desoto'
meta.loc[meta['location'].str.contains('pmutsid'), 'location'] = 'pmutsidio'
meta.loc[meta['location'].str.contains('lasalle'), 'location'] = 'la salle'
meta.loc[meta['location'].str.contains('camer'), 'location'] = 'cameron'
meta.loc[meta['location'].str.contains('cast'), 'location'] = 'castro'
meta.loc[meta['location'].str.contains('cros'), 'location'] = 'crosby'
meta.loc[meta['location'].str.contains('ect'), 'location'] = 'ector'
meta.loc[meta['location'].str.contains('galvest'), 'location'] = 'galveston'
meta.loc[meta['location'].str.contains('jeffers'), 'location'] = 'jefferson'
meta.loc[meta['location'].str.contains('kaufm'), 'location'] = 'kaufman'
meta.loc[meta['location'].str.contains('klebe'), 'location'] = 'kleberg'
meta.loc[meta['location'].str.contains('lava'), 'location'] = 'lavaca'
meta.loc[meta['location'].str.contains('mclenn'), 'location'] = 'mclennan'
meta.loc[meta['location'].str.contains('st.clair'), 'location'] = 'saint clair'
meta.loc[meta['location'].str.contains('tarra'), 'location'] = 'tarrant'
meta.loc[meta['location'].str.contains('wall'), 'location'] = 'waller'
meta.loc[meta['location'].str.contains('wichi'), 'location'] = 'wichita'

# final cleaning (missing values)
meta.loc[meta['location']=='unk', 'location'] = unknown_val+longstring
meta.loc[meta['division']==meta['country'], 'division'] = unknown_val+longstring
meta.fillna(unknown_val, inplace=True)
# generate json
with open('{}/gadm_countries.json'.format(geojson_prefix)) as f:
    countries = json.load(f)
with open('{}/gadm_divisions.json'.format(geojson_prefix)) as f:
    divisions = json.load(f)
with open('{}/gadm_locations.json'.format(geojson_prefix)) as f:
    locations = json.load(f)

meta.loc[meta['country'] == 'usa', 'country'] = 'united states'
meta['country_id'] = meta['country'].apply(lambda x: sorted([ (fuzz.ratio(key, x), (key, id)) for key, id in countries.items() ])[-1])
meta.loc[meta['country_id'].apply(lambda x: x[0]) < 80, 'country_id'] = unknown_val
meta.loc[meta['country_id'] != unknown_val, 'country'] = meta.loc[meta['country_id'] != unknown_val, 'country_id'].apply(lambda x: x[1][0]).str.lower()
meta.loc[meta['country_id'] != unknown_val, 'country_id'] = meta.loc[meta['country_id'] != unknown_val, 'country_id'].apply(lambda x: x[1][1]).str.upper()
meta.loc[meta['country_id'] == unknown_val.upper(), 'country_id'] = unknown_val
meta['tmp_info1'] = meta['country'] + '-' + meta['division']
meta['division_id'] = meta['tmp_info1'].apply(lambda x: sorted([ (int(x.split('-')[0] in key) * fuzz.ratio(key, x), (key, id)) for key, id in divisions.items() ])[-1])
meta.loc[meta['division'] == unknown_val, 'division'] = unknown_val+longstring
meta.loc[meta['division_id'].apply(lambda x: x[0]) < 80, 'division_id'] = unknown_val
meta.loc[meta['division_id'] != unknown_val, 'division'] = meta.loc[meta['division_id'] != unknown_val, 'division_id'].apply(lambda x: x[1][0].split("-")[1]).str.lower()
meta.loc[meta['division_id'] != unknown_val, 'division_id'] = meta.loc[meta['division_id'] != unknown_val, 'division_id'].apply(lambda x: x[1][1]).str.upper()
meta.loc[meta['division_id'] == unknown_val.upper(), 'division_id'] = unknown_val
meta['tmp_info2'] = meta['country'] + '-' + meta['division'] + '-' + meta['location']
meta['location_id'] = meta['tmp_info2'].apply(lambda x: locations.get(x, unknown_val)).astype(str)

meta['date_modified'] = date_modified
meta['locstring'] = meta['locstring'].str.split("/")
meta['country'] = meta['locstring'].apply(lambda x: x[1] if len(x) >= 2 else unknown_val).str.strip()
meta['division'] = meta['locstring'].apply(lambda x: x[2] if len(x) >= 3 else unknown_val).str.strip()
meta['location'] = meta['locstring'].apply(lambda x: x[3] if len(x) >= 4 else unknown_val).str.strip()
meta.loc[meta['country'] == 'USA', 'country'] = 'United States'
meta['country_lower'] = meta['country'].str.lower()
meta['division_lower'] = meta['division'].str.lower()
meta['location_lower'] = meta['location'].str.lower()
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
