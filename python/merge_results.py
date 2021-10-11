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


meta = pd.read_csv(input_metadata, sep='\t')
# normalize date information
meta['tmp'] = meta['date_collected'].str.split('-')
meta = meta[meta['tmp'].str.len()>=2]
meta.loc[meta['tmp'].str.len()==2, 'date_collected'] += '-15'
meta['date_collected'] = pd.to_datetime(meta['date_collected'], errors='coerce')
meta['date_collected'] = meta['date_collected'].astype(str)
meta = meta[meta['date_collected'] <= max_date]
meta = meta[meta['date_collected'] > min_date]
# normalize geo information
meta['country'] = meta['country'].copy()
meta['country'].fillna('None', inplace=True)
meta.loc[meta['country'].str.lower() == 'usa',
        'country'] = 'United States'
meta.loc[(meta['country'].str.contains('Congo')) & (meta['country'].str.contains(
    'emocratic')), 'country'] = 'Democratic Republic of the Congo'
meta.loc[(meta['country'].str.contains('Congo')) & ~(meta['country'].str.contains(
    'emocratic')), 'country'] = 'Republic of Congo'
meta.loc[meta['country'].str.contains(
    'Eswatini'), 'country'] = "Swaziland"
meta.loc[meta['country'].str.contains(
    'Bonaire'), 'country'] = "Bonaire, Sint Eustatius and Saba"
meta.loc[meta['country'].str.contains(
    'Sint Eustatius'), 'country'] = "Bonaire, Sint Eustatius and Saba"
meta.loc[meta['country'].str.contains('Cote dIvoire'), 'country'] = "Côte d'Ivoire"
meta.loc[meta['country'].str.contains("Cote d'Ivoire"), 'country'] = "Côte d'Ivoire"
meta.loc[meta['country'].str.contains("México"), 'country'] = "Mexico"
meta.loc[meta['country'].str.contains('North Macedonia'), 'country'] = "Macedonia"
meta.loc[meta['country'].str.contains('Curacao'), 'country'] = "Curaçao"
meta.loc[meta['country'].str.contains('Saint Martin'), 'country'] = "Saint-Martin"
meta.loc[meta['country'].str.contains('Trinidad'), 'country'] = 'Trinidad and Tobago'
meta.loc[meta['country'].str.contains('Czech republic'), 'country'] = 'Czech Republic'
meta.loc[meta['country'].str.contains('St Eustatius'), 'country'] = 'Netherlands'
meta.loc[meta['country'].str.contains('Saint Barthelemy'), 'country'] = 'Saint-Barthélemy'
meta.loc[meta['country'].str.contains('Palestine'), 'country'] = "Palestina"
meta.loc[meta['country'].str.contains("Germany /"), 'country'] = "Germany"
meta.loc[meta['country'].str.contains("France /Nouvelle-Aquitaine"), 'division'] = "Nouvelle-Aquitaine"
meta.loc[meta['country']=="France /Nouvelle-Aquitaine", 'country'] = "France"
meta.loc[meta['country'].str.contains("France /Nouvelle-Aquitaine/ Limoges"), 'division'] = "Nouvelle-Aquitaine"
meta.loc[meta['country'].str.contains("France /Nouvelle-Aquitaine/ Limoges"), 'location'] = "Limoges"
meta.loc[meta['country']=="France /Nouvelle-Aquitaine/ Limoges", 'country'] = "France"
meta.loc[meta['country']=="Kenya /", 'country'] = "Kenya"
meta.loc[meta['country']=="Switzerland/ Schwyz", 'division'] = "Schwyz"
meta.loc[meta['country']=="Switzerland/ Schwyz", 'country'] = "Switzerland"
meta.loc[meta['country']=="USA /Wisconsin", 'division'] = "Wisconsin"
meta.loc[meta['country']=="USA /Wisconsin", 'country'] = "United States"
meta.loc[meta['country']=="Jonavos apskritis", 'country'] = "Lithuania"
meta.loc[meta['country']=="Thailand /Singburi", 'division'] = "Singburi"
meta.loc[meta['country']=="Thailand /Singburi", 'country'] = "Thailand"
meta.loc[meta['country']=="Norway /", 'country'] = "Norway"
meta.loc[meta['country']=="Morocoo", 'country'] = "Morocco"
print(f'Admin1 standardization...')
meta.loc[meta['division'].isna(), 'division'] = 'None'
meta['division'] = meta['division'].copy()
meta.loc[meta['division']=='USA', 'division'] = 'United States'
meta.loc[meta['division'].str.contains('Georgia /'), 'division'] = 'Georgia'
meta.loc[meta['division'].str.contains('Antwerp'), 'division'] = 'Vlaanderen'
meta.loc[meta['division'].str.contains('Andalu'), 'division'] = 'Andalucía'
meta.loc[meta['division'].str.contains('Cairo'), 'division'] = 'Al Qahirah'
meta.loc[meta['division'].str.contains('Northern territory'), 'division'] = 'Northern Territory'
meta.loc[meta['division'].str.contains('Fayoum'), 'division'] = 'Al Fayyum'
meta.loc[meta['division'].str.contains('Musca'), 'division'] = 'Muscat'
meta.loc[meta['division'].str.contains('Kalyoubia'), 'division'] = 'Al Qalyubiyah'
meta.loc[meta['division'].str.contains('Buraymi'), 'division'] = 'Al Buraymi'
meta.loc[meta['division'].str.contains('Buraimi'), 'division'] = 'Al Buraymi'
meta.loc[meta['division'].str.contains('Dakhiliyah'), 'division'] = 'Ad Dakhliyah'
meta.loc[meta['division'].str.contains('Dhahirah'), 'division'] = 'Al Dhahira'
meta.loc[meta['division'].str.contains('North Batinah'), 'division'] = 'Al Batinah North'
meta.loc[meta['division'].str.contains('South Batinah'), 'division'] = 'Al Batinah South'
meta.loc[meta['division'].str.contains('North Sharqiyah'), 'division'] = 'Ash Sharqiyah North'
meta.loc[meta['division'].str.contains('Wuhan'), 'division'] = 'Hubei'
meta.loc[meta['division'].str.contains('Quebec'), 'division'] = 'Québec'
meta.loc[meta['division'].str.contains('Toronto'), 'division'] = 'Ontario'
meta.loc[meta['division'].str.contains('Coahuila de Zaragoza'), 'division'] = 'Coahuila'
meta.loc[meta['division'].str.contains('Mexico City'), 'division'] = 'México'
meta.loc[meta['division'].str.contains('Michoacan'), 'division'] = 'Michoacán'
meta.loc[meta['division'].str.contains('Nuevo Leon'), 'division'] = 'Nuevo León'
meta.loc[meta['division'].str.contains('Queretaro'), 'division'] = 'Querétaro'
meta.loc[meta['division'].str.contains('SanLuisPotosi'), 'division'] = 'San Luis Potosí'
meta.loc[meta['division'].str.contains('San Luis Potosi'), 'division'] = 'San Luis Potosí'
meta.loc[meta['division'].str.contains('State of Mexico'), 'division'] = 'México'
meta.loc[meta['division'].str.contains('Yucatan'), 'division'] = 'Yucatán'
meta.loc[meta['division'].str.contains('Bethlehem'), 'division'] = 'West Bank'
meta.loc[meta['division'].str.contains('Hebron'), 'division'] = 'West Bank'
meta.loc[meta['division'].str.contains('Jenin'), 'division'] = 'West Bank'
meta.loc[meta['division'].str.contains('Jericho'), 'division'] = 'West Bank'
meta.loc[meta['division'].str.contains('Ramallah'), 'division'] = 'West Bank'
meta.loc[meta['division'].str.contains('Tulkarem'), 'division'] = 'West Bank'
meta.loc[meta['division'].str.contains('Nablus'), 'division'] = 'West Bank'
meta.loc[meta['division'].str.contains('Sharja'), 'division'] = 'Sharjah'
meta.loc[meta['division'].str.contains('Copenhagen'), 'division'] = 'Hovedstaden'
meta.loc[meta['division'].str.contains('Sjaelland'), 'division'] = 'Sjælland'
meta.loc[meta['division'].str.contains('Cape Town'), 'division'] = 'Western Cape'
meta.loc[meta['division'].str.contains('Western Cape'), 'division'] = 'Western Cape'
meta.loc[meta['division'].str.contains('Amapa'), 'division'] = 'Amapá'
meta.loc[meta['division'].str.contains('Ceara'), 'division'] = 'Ceará'
meta.loc[meta['division'].str.contains('Goias'), 'division'] = 'Goiás'
meta.loc[meta['division'].str.contains('Maranhao'), 'division'] = 'Maranhão'
meta.loc[meta['division'].str.contains('Paraiba'), 'division'] = 'Paraíba'
meta.loc[meta['division'].str.contains('Parana'), 'division'] = 'Paraná'
meta.loc[meta['division'].str.contains('Piaui'), 'division'] = 'Piauí'
meta.loc[meta['division'].str.contains('Sao Paulo'), 'division'] = 'São Paulo'
meta.loc[meta['division'].str.contains('Aragon'), 'division'] = 'Aragón'
meta.loc[meta['division'].str.contains('Asturias'), 'division'] = 'Principado de Asturias'
meta.loc[meta['division'].str.contains('Balear Islands'), 'division'] = 'Islas Baleamuts'
meta.loc[meta['division'].str.contains('Balear_Islands'), 'division'] = 'Islas Baleamuts'
meta.loc[meta['division'].str.contains('Illes Balears'), 'division'] = 'Islas Baleamuts'
meta.loc[meta['division'].str.contains('Canary Islands'), 'division'] = 'Canaries'
meta.loc[meta['division'].str.contains('Canary_Islands'), 'division'] = 'Canaries'
meta.loc[meta['division'].str.contains('Castilla La Mancha'), 'division'] = 'Castilla-La Mancha'
meta.loc[meta['division'].str.contains('Castilla la Mancha'), 'division'] = 'Castilla-La Mancha'
meta.loc[meta['division'].str.contains('Castilla y Leon'), 'division'] = 'Castilla y León'
meta.loc[meta['division'].str.contains('Ceuta'), 'division'] = 'Ceuta y Melilla'
meta.loc[meta['division'].str.contains('Melilla'), 'division'] = 'Ceuta y Melilla'
meta.loc[meta['division'].str.contains('Comunitat Valenciana'), 'division'] = 'Comunidad Valenciana'
meta.loc[meta['division'].str.contains('Comunitat_Valenciana'), 'division'] = 'Comunidad Valenciana'
meta.loc[meta['division'].str.contains('La_Rioja'), 'division'] = 'La Rioja'
meta.loc[meta['division'].str.contains('Madrid'), 'division'] = 'Comunidad de Madrid'
meta.loc[meta['division'].str.contains('Murcia'), 'division'] = 'Región de Murcia'
meta.loc[meta['division'].str.contains('Navarra'), 'division'] = 'Comunidad Foral de Navarra'
meta.loc[meta['division'].str.contains('Catalunya'), 'division'] = 'Cataluña'
meta.loc[meta['division'].str.contains('Catalonia'), 'division'] = 'Cataluña'
# Germany
meta.loc[meta['division'].str.contains('Baden-Wuerttemberg'), 'division'] = 'Baden-Württemberg'
meta.loc[meta['division'].str.contains('Baden-Wurttemberg'), 'division'] = 'Baden-Württemberg'
meta.loc[meta['division'].str.contains('Bavaria'), 'division'] = 'Bayern'
meta.loc[meta['division'].str.contains('Hesse'), 'division'] = 'Hessen'
meta.loc[meta['division'].str.contains('Lower Saxony'), 'division'] = 'Niedersachsen'
meta.loc[meta['division'].str.contains('Mecklenburg-Western Pomerania'), 'division'] = 'Mecklenburg-Vorpommern'
meta.loc[meta['division'].str.contains('Rhineland-Palatinate'), 'division'] = 'Rheinland-Pfalz'
meta.loc[(meta['division'].str.contains('Saxony-Anhalt'))
        & (meta['country'].str.contains('Germany')), 'division'] = 'Sachsen-Anhalt'
meta.loc[(meta['division'].str.contains('Saxony'))
        & (meta['country'].str.contains('Germany')), 'division'] = 'Sachsen'
meta.loc[meta['division'].str.contains('North Rhine-Westphalia'), 'division'] = 'Nordrhein-Westfalen'
meta.loc[meta['division'].str.contains('Thuringia'), 'division'] = 'Thüringen'
# South Africa
meta.loc[(meta['country'].str.contains('South Africa'))
        & (meta['division'].str.contains('KwaZulu Natal')),
        'division'] = 'KwaZulu-Natal'
meta.loc[(meta['country'].str.contains('South Africa'))
        & (meta['division'].str.contains('Northern Cape Province')),
        'division'] = 'Northern Cape'
# Austria
meta.loc[(meta['country']=='Austria')
        & (meta['division']=='Tyrol'), 'division'] = 'Tirol'
# Switzerland
meta.loc[meta['division'].str.contains('Argovie'), 'division'] = 'Aargau'
meta.loc[meta['division'].str.contains('Geneva'), 'division'] = 'Genève'
meta.loc[meta['division'].str.contains('Grabunden'), 'division'] = 'Graubünden'
meta.loc[meta['division'].str.contains('Luzern'), 'division'] = 'Lucerne'
meta.loc[meta['division'].str.contains('Neuchatel'), 'division'] = 'Neuchâtel'
meta.loc[meta['division'].str.contains('Neuenburg'), 'division'] = 'Neuchâtel'
meta.loc[meta['division'].str.contains('Obwald'), 'division'] = 'Obwalden'
meta.loc[meta['division'].str.contains('Saint-Gall'), 'division'] = 'Sankt Gallen'
meta.loc[meta['division'].str.contains('St Gallen'), 'division'] = 'Sankt Gallen'
meta.loc[meta['division'].str.contains('St. Gallen'), 'division'] = 'Sankt Gallen'
meta.loc[meta['division'].str.contains('Schaffhouse'), 'division'] = 'Schaffhausen'
meta.loc[meta['division'].str.contains('Turgovia'), 'division'] = 'Thurgau'
meta.loc[meta['division'].str.contains('VALAIS'), 'division'] = 'Valais'
meta.loc[meta['division'].str.contains('Waadt'), 'division'] = 'Vaud'
meta.loc[meta['division'].str.contains('Wallis'), 'division'] = 'Valais'
meta.loc[meta['division'].str.contains('Zaerich'), 'division'] = 'Zürich'
meta.loc[meta['division'].str.contains('Zoerich'), 'division'] = 'Zürich'
meta.loc[meta['division'].str.contains('Zurich'), 'division'] = 'Zürich'
meta.loc[meta['division'].str.contains('Unkown'), 'division'] = 'Unknown'
meta.loc[meta['division'].str.contains('Basel-Land'), 'division'] = 'Basel-Landschaft'
meta.loc[meta['division'].str.contains('Grisons'), 'division'] = 'Graubünden'
meta.loc[meta['division'].str.contains('St.Gallen'), 'division'] = 'Sankt Gallen'
meta.loc[meta['division'].str.contains('Geneve'), 'division'] = 'Genève'
meta.loc[meta['division'].str.contains('Stankt Gallen'), 'division'] = 'Sankt Gallen'
# Sweden
meta.loc[(meta['country']=='Sweden')
        &(meta['division']=='Sodermanland'), 'division'] = 'Södermanland'
meta.loc[(meta['country']=='Sweden')
        &(meta['division']=='Gavleborg'), 'division'] = 'Gävleborg'
meta.loc[(meta['country']=='Sweden')
        &(meta['division']=='Jamtland'), 'division'] = 'Jämtland'
meta.loc[(meta['country']=='Sweden')
        &(meta['division']=='Jonkoping'), 'division'] = 'Jönköping'
meta.loc[(meta['country']=='Sweden')
        &(meta['division']=='Ostergotland'), 'division'] = 'Östergötland'
meta.loc[(meta['country']=='Sweden')
        &(meta['division']=='Skane'), 'division'] = 'Skåne'
meta.loc[(meta['country']=='Sweden')
        &(meta['division']=='Varmland'), 'division'] = 'Värmland'
meta.loc[(meta['country']=='Sweden')
        &(meta['division']=='Vasterbotten'), 'division'] = 'Västerbotten'
meta.loc[(meta['country']=='Sweden')
        &(meta['division']=='Vasternorrland'), 'division'] = 'Västernorrland'
meta.loc[(meta['country']=='Sweden')
        &(meta['division']=='Vastmanland'), 'division'] = 'Västmanland'
meta.loc[(meta['country']=='Sweden')
        &(meta['division']=='Vastra Gotaland'), 'division'] = 'Västra Götaland'
# Puerto Rico
meta.loc[(meta['division']=='Puerto Rico'), 'country'] = 'Puerto Rico'
meta.loc[(meta['country']=='Puerto Rico'), 'division'] = meta.loc[(meta['country']=='Puerto Rico'), 'location'].fillna('None')
print(f'Admin2 standardization (U.S. only)')
meta.loc[meta['location'].isna(), 'location'] = 'None'
meta['location'] = meta['location'].copy()
houston_filter = (meta['division']=='Texas') & (meta['location']=='Houston')
meta.loc[houston_filter, 'location'] = 'Harris County'
for key, val in COUNTY_CORRECTIONS.items():
    meta.loc[:, 'location'] = meta['location'].str.replace(key, val)
meta.loc[:, 'location'] = meta['location'].str.replace('County', '').str.replace('county', '').str.replace(',', '')
meta.loc[:, 'location'] = meta['location'].str.strip().apply(bs.check_state, args=(False,)).str.strip()
meta.loc[meta['location'].str.contains('Anchorage-Mat-Su'), 'location'] = 'Anchorage'
meta.loc[meta['location'].str.contains('Anchorage-Mat Su'), 'location'] = 'Anchorage'
meta.loc[meta['location'].str.contains('BRA'), 'location'] = 'Brazos'
meta.loc[meta['location'].str.contains('BR'), 'location'] = 'Brewster'
meta.loc[meta['location'].str.contains('Belgrade'), 'location'] = 'Gallatin'
meta.loc[meta['location'].str.contains('Bozeman'), 'location'] = 'Gallatin'
meta.loc[meta['location'].str.contains('Big Sky'), 'location'] = 'Gallatin'
meta.loc[meta['location'].str.contains('Belton'), 'location'] = 'Bell'
meta.loc[meta['location'].str.contains('Brentwood'), 'location'] = 'Contra Costa'
meta.loc[meta['location'].str.contains('Chicago'), 'location'] = 'Cook'
meta.loc[meta['location'].str.contains('Colombus'), 'location'] = 'Franklin'
meta.loc[meta['location'].str.contains('DuBois'), 'location'] = 'Fremont'
meta.loc[meta['location'].str.contains('DuPage'), 'location'] = 'Dupage'
meta.loc[meta['location'].str.contains('Eau claire'), 'location'] = 'Eau Claire'
meta.loc[meta['location'].str.contains('Ennis'), 'location'] = 'Ellis'
meta.loc[meta['location'].str.contains('Fond Du Lac'), 'location'] = 'Fond du Lac'
meta.loc[meta['location'].str.contains('Fond du lac'), 'location'] = 'Fond du Lac'
meta.loc[meta['location'].str.contains('Fonddu Lac'), 'location'] = 'Fond du Lac'
meta.loc[meta['location'].str.contains('Frisco'), 'location'] = 'Collin'
meta.loc[meta['location'].str.contains('Hawai'), 'location'] = 'Hawaii'
meta.loc[meta['location'].str.contains('Holland'), 'location'] = 'Ottawa'
meta.loc[meta['location'].str.contains('Honolul'), 'location'] = 'Honolulu'
meta.loc[meta['location'].str.contains('Indianapolis'), 'location'] = 'Marion'
meta.loc[meta['location'].str.contains('Interior'), 'location'] = 'Fairbanks North Star'
meta.loc[meta['location'].str.contains('Ithaca'), 'location'] = 'Tompkins'
meta.loc[meta['location'].str.contains('Kaua'), 'location'] = 'Kauai'
meta.loc[meta['location'].str.contains('Las Vegas'), 'location'] = 'Clark'
meta.loc[meta['location'].str.contains('Mau'), 'location'] = 'Hawaii'
meta.loc[meta['location'].str.contains('Mcculloch'), 'location'] = 'McCulloch'
meta.loc[meta['location'].str.contains('Mchenry'), 'location'] = 'McHenry'
meta.loc[meta['location'].str.contains('Mclennan'), 'location'] = 'McLennan'
meta.loc[meta['location'].str.contains('Moris'), 'location'] = 'Morris'
meta.loc[meta['location'].str.contains('New York'), 'location'] = 'New York'
meta.loc[meta['location'].str.contains('New York City'), 'location'] = 'New York'
meta.loc[meta['location'].str.contains('New Hyde Park'), 'location'] = 'Nassau'
meta.loc[meta['location'].str.contains('New Orleans'), 'location'] = 'Orleans'
meta.loc[meta['location'].str.contains('New Rochelle'), 'location'] = 'Westchester'
meta.loc[meta['location'].str.contains('Northern'), 'location'] = 'Fairbanks North Star'
meta.loc[meta['location'].str.contains('Omaha'), 'location'] = 'Douglas'
meta.loc[meta['location'].str.contains('Ostego'), 'location'] = 'Allegan'
meta.loc[meta['location'].str.contains('Phoenix'), 'location'] = 'Maricopa'
meta.loc[meta['location'].str.contains('San Bernadino'), 'location'] = 'San Bernardino'
meta.loc[meta['location'].str.contains('Seattle'), 'location'] = 'King'
meta.loc[meta['location'].str.contains('St. Bernard'), 'location'] = 'Saint Bernard'
meta.loc[meta['location'].str.contains('St. Clair'), 'location'] = 'Saint Clair'
meta.loc[meta['location'].str.contains('St. Lawrence'), 'location'] = 'Saint Lawrence'
meta.loc[meta['location'].str.contains('St. Louis'), 'location'] = 'Saint Louis'
meta.loc[meta['location'].str.contains('St. Tammany'), 'location'] = 'Saint Tammany'
meta.loc[meta['location'].str.contains('Staten Island'), 'location'] = 'Richmond'
meta.loc[meta['location'].str.contains('Thurson'), 'location'] = 'Thurston'
meta.loc[meta['location'].str.contains('Tucson'), 'location'] = 'Pima'
meta.loc[meta['location'].str.contains('West Yellowstone'), 'location'] = 'Gallatin'
meta.loc[meta['location'].str.contains('Adam'), 'location'] = 'Adams'
meta.loc[meta['location'].str.contains('Alachu'), 'location'] = 'Alachua'
meta.loc[meta['location'].str.contains('Du Bois'), 'location'] = 'Dubois'
meta.loc[meta['location'].str.contains('DeSoto'), 'location'] = 'Desoto'
meta.loc[meta['location'].str.contains('PmutsID'), 'location'] = 'Pmutsidio'
meta.loc[meta['location'].str.contains('LaSalle'), 'location'] = 'La Salle'
meta.loc[meta['location'].str.contains('CAMER'), 'location'] = 'Cameron'
meta.loc[meta['location'].str.contains('CAST'), 'location'] = 'Castro'
meta.loc[meta['location'].str.contains('CROS'), 'location'] = 'Crosby'
meta.loc[meta['location'].str.contains('ECT'), 'location'] = 'Ector'
meta.loc[meta['location'].str.contains('GALVEST'), 'location'] = 'Galveston'
meta.loc[meta['location'].str.contains('JEFFERS'), 'location'] = 'Jefferson'
meta.loc[meta['location'].str.contains('KAUFM'), 'location'] = 'Kaufman'
meta.loc[meta['location'].str.contains('KLEBE'), 'location'] = 'Kleberg'
meta.loc[meta['location'].str.contains('LAVA'), 'location'] = 'Lavaca'
meta.loc[meta['location'].str.contains('MCLENN'), 'location'] = 'Mclennan'
meta.loc[meta['location'].str.contains('St.Clair'), 'location'] = 'Saint Clair'
meta.loc[meta['location'].str.contains('TARRA'), 'location'] = 'Tarrant'
meta.loc[meta['location'].str.contains('WALL'), 'location'] = 'Waller'
meta.loc[meta['location'].str.contains('WICHI'), 'location'] = 'Wichita'
# TODO: France, Russia, China, Israel, South Korea
country = 'USA'
if country:
    gisaid_2 = set(meta[meta['country']==country]['location'].unique())
else:
    gisaid_2 = set(meta['location'].unique())

# final cleaning (missing values)
meta.loc[meta['location']=='unk', 'location'] = unknown_val
meta.loc[meta['division']==meta['country'], 'division'] = unknown_val
meta.fillna(unknown_val, inplace=True)
# generate json
with open('{}/gadm_countries.json'.format(geojson_prefix)) as f:
    countries = json.load(f)
with open('{}/gadm_divisions.json'.format(geojson_prefix)) as f:
    divisions = json.load(f)
with open('{}/gadm_locations.json'.format(geojson_prefix)) as f:
    locations = json.load(f)
meta['country_id'] = meta['country'].apply(lambda x: countries.get(x, unknown_val)).astype(str)
meta['tmp_info1'] = meta['country'] + '-' + meta['division']
meta['division_id'] = meta['tmp_info1'].apply(lambda x: divisions.get(x, unknown_val)).astype(str)
meta['tmp_info2'] = meta['country'] + '-' + meta['division'] + '-' + meta['location']
meta['location_id'] = meta['tmp_info2'].apply(lambda x: locations.get(x, unknown_val)).astype(str)
meta['date_modified'] = date_modified
meta_info = ['strain', 'accession_id', 'pangolin_lineage', 'date_collected', 'country_id', 'division_id', 'location_id']

print(meta_info)

muts = pd.read_csv(input_mut, dtype=str)
muts = muts[~(muts['gene'].isin(['5UTR', '3UTR']))]
# ignore mutations found in non-coding regions
muts = muts.loc[~(muts['gene']=='Non-coding region')]
# fuse with metadata
print(f"Fusing with metadata...")
# concat with pd
muts = pd.merge(muts, meta, left_on='idx', right_on='strain')
muts = muts.drop_duplicates(subset=['accession_id', 'mutation'])
# If deletions not in chunk add columns
del_columns = ['is_frameshift', 'change_length_nt', 'deletion_codon_coords', 'absolute_coords']
muts_columns = muts.columns.tolist()
for i in del_columns:
    if i not in muts_columns:
        muts[i] = np.nan
muts_info = [
    'type', 'mutation', 'gene',
    'ref_codon', 'pos', 'alt_codon',
    'is_synonymous',
    'ref_aa', 'codon_num', 'alt_aa',
    'absolute_coords',
    'change_length_nt', 'is_frameshift',
    'deletion_codon_coords'
]

print(muts)

# GENERATE JSON DATA MODEL
(
    muts.groupby(meta_info, as_index=True)
             .apply(lambda x: x[muts_info].to_dict('records'))
             .reset_index()
             .rename(columns={0:'mutations'})
             .to_json(
                 out_fp,
                 orient='records',
                 lines = True,
                 compression = "gzip"
             )
 )
