#!/usr/bin/env python
import gc
import re
import json
import argparse
import pandas as pd
import geopandas as gpd
from path import Path
import bjorn_support as bs
import data as bd




def download_process_data(username, password, chunk_size, current_datetime):
    with open('config.json', 'r') as f:
        config = json.load(f)
    out_dir = Path(config['out_dir'])
    in_fp = out_dir/(config['gisaid_feed'] + '_' + current_datetime + '.json')
    out_fp = out_dir/(config['gisaid_fasta'] + '_' + current_datetime + '.fasta')
    meta_fp = out_dir/(config['gisaid_meta'] + '_' + current_datetime + '.tsv.gz')
    info_fp = out_dir/config['chunk_info']
    gadm_fp = config['gadm']
    is_test = config['feed_test']
    test_size = 100
    metacols = ['covv_virus_name', 'covsurver_prot_mutations', 'covv_location',
                'covv_lineage', 'covv_collection_date', 'covv_accession_id',
                'pangolin_lineages_version', 'covv_clade', 'covv_subm_date']
    corrections = bd.COUNTY_CORRECTIONS
    # Download GISAID API feed
    if not Path.isfile(Path(in_fp)):
        print(f"Downloading GISAID API feed...")
        feed_cmd = f"curl -u {username}:{password} https://www.epicov.org/epi3/3p/scripps/export/provision.json.xz | xz -d -T8 > {in_fp}"
        bs.run_command(feed_cmd)
        print(f"GISAID API feed saved in {in_fp}")
    else:
        print(f"{in_fp} already exists, skipping API feed download step")
    # load sequence data
    print(f"Loading API data...")
    if is_test:
        data = []
        with open(in_fp, 'r') as handle:
            for i, line in enumerate(handle):
                data.append(json.loads(line))
                if i >= test_size:
                    break
    else:
        if Path.isfile(Path(out_fp)) and Path.isfile(Path(meta_fp)):
            data = []
        else:
            data = [json.loads(line) for line in open(in_fp, 'r')]
    print(f"Total number of sequences: {len(data)}")
    # generate fasta file containing all sequences
    if not Path.isfile(Path(out_fp)):
        print(f"Converting to dict...")
        regex = re.compile('[^a-zA-Z]')
        seqs_dict = {sample['covv_virus_name'].replace('hCoV-19/', '').replace(' ', ''): 
                    regex.sub('', sample['sequence'].replace('\n', '')) for sample in data}
        print(f"Converting to FASTA...")
        bs.dict2fasta(seqs_dict, out_fp)
        print(f"FASTA output generated and saved in {out_fp}")
    else:
        print(f"{out_fp} already exists, skipping fasta generation step")
    # generate tsv file containing processed metadata
    if not Path.isfile(Path(meta_fp)):
        # load raw metadata into dataframe
        df = pd.DataFrame(data, columns=metacols)
        # TEST: all accession IDs are unique
        assert df['covv_accession_id'].shape[0]==df['covv_accession_id'].unique().shape[0], f'ERROR: gisaid accession ids not unique'
        num_ids = df['covv_accession_id'].unique().shape[0]
        print(f"Total number of sequences: {num_ids}")
        print(f"Cleaning metadata")
        df.rename(columns={
                        'covv_virus_name': 'strain', 
                        'covv_location': 'location', 
                        'covv_collection_date': 'date_collected',
                        'covv_subm_date': 'date_submitted',
                        'covv_clade': 'clade',
                        'covv_lineage': 'pangolin_lineage',
                        'pangolin_lineages_version': 'pangolin_version',
                        'covv_accession_id': 'accession_id'
                        }, inplace=True)
        print(f"Loading GADM file containing geographic information...")
        gadm = gpd.read_file(gadm_fp)
        gadm_cols = [f'NAME_{i}' for i in range(5)]
        gadm = gadm[gadm_cols]
        print(f"Standardizing location information...")
        res = pd.DataFrame(df['location'].str.split('/').tolist(), 
                    columns=['region',
                            'country', 
                            'division', 
                            'location', 
                            'city', 
                            'town'
                            ])
        df['country'] = res['country'].str.strip()
        df['division'] = res['division'].str.strip()
        df['location'] = res['location'].str.strip()
        print(f"Admin0 standardization...")
        df['country_normed'] = df['country'].copy()
        df['country_normed'].fillna('None', inplace=True)
        df.loc[df['country_normed'].str.lower()=='usa', 'country_normed'] = 'United States'     
        df.loc[(df['country'].str.contains('Congo')) & (df['country'].str.contains('Democratic')), 'country_normed'] = 'Democratic Republic of the Congo'
        df.loc[(df['country'].str.contains('Congo')) & (~df['country'].str.contains('Democratic')), 'country_normed'] = 'Republic of Congo'
        df.loc[df['country_normed'].str.contains('Eswatini'), 'country_normed'] = "Swaziland"
        df.loc[df['country_normed'].str.contains('Bonaire'), 'country_normed'] = "Bonaire, Sint Eustatius and Saba"
        df.loc[df['country_normed'].str.contains('Sint Eustatius'), 'country_normed'] = "Bonaire, Sint Eustatius and Saba"
        df.loc[df['country_normed'].str.contains('Cote dIvoire'), 'country_normed'] = "Côte d'Ivoire"
        df.loc[df['country_normed'].str.contains("Cote d'Ivoire"), 'country_normed'] = "Côte d'Ivoire"
        df.loc[df['country_normed'].str.contains("México"), 'country_normed'] = "Mexico"
        df.loc[df['country_normed'].str.contains('North Macedonia'), 'country_normed'] = "Macedonia"
        df.loc[df['country_normed'].str.contains('Curacao'), 'country_normed'] = "Curaçao"
        df.loc[df['country_normed'].str.contains('Saint Martin'), 'country_normed'] = "Saint-Martin"
        df.loc[df['country_normed'].str.contains('Trinidad'), 'country_normed'] = 'Trinidad and Tobago'
        df.loc[df['country_normed'].str.contains('Czech republic'), 'country_normed'] = 'Czech Republic'
        df.loc[df['country_normed'].str.contains('St Eustatius'), 'country_normed'] = 'Netherlands'
        df.loc[df['country_normed'].str.contains('Saint Barthelemy'), 'country_normed'] = 'Saint-Barthélemy'
        df.loc[df['country_normed'].str.contains('Palestine'), 'country_normed'] = "Palestina"
        df.loc[df['country_normed'].str.contains("Germany /"), 'country_normed'] = "Germany"
        df.loc[df['country_normed'].str.contains("France /Nouvelle-Aquitaine"), 'division'] = "Nouvelle-Aquitaine"
        df.loc[df['country_normed']=="France /Nouvelle-Aquitaine", 'country_normed'] = "France"
        df.loc[df['country_normed'].str.contains("France /Nouvelle-Aquitaine/ Limoges"), 'division'] = "Nouvelle-Aquitaine"
        df.loc[df['country_normed'].str.contains("France /Nouvelle-Aquitaine/ Limoges"), 'location'] = "Limoges"
        df.loc[df['country_normed']=="France /Nouvelle-Aquitaine/ Limoges", 'country_normed'] = "France"
        df.loc[df['country_normed']=="Kenya /", 'country_normed'] = "Kenya"
        df.loc[df['country_normed']=="Switzerland/ Schwyz", 'division'] = "Schwyz"
        df.loc[df['country_normed']=="Switzerland/ Schwyz", 'country_normed'] = "Switzerland"
        df.loc[df['country_normed']=="USA /Wisconsin", 'division'] = "Wisconsin"
        df.loc[df['country_normed']=="USA /Wisconsin", 'country_normed'] = "United States"
        df.loc[df['country_normed']=="Jonavos apskritis", 'country_normed'] = "Lithuania"
        df.loc[df['country_normed']=="Thailand /Singburi", 'division'] = "Singburi"
        df.loc[df['country_normed']=="Thailand /Singburi", 'country_normed'] = "Thailand"
        df.loc[df['country_normed']=="Norway /", 'country_normed'] = "Norway"
        df.loc[df['country_normed']=="Morocoo", 'country_normed'] = "Morocco"
        gisaid_0 = set(df['country_normed'].unique())
        gadm_0 = set(gadm['NAME_0'].unique())
        print(len(gisaid_0))
        print(len(gadm_0))
        print(f'Number of countries captured in GADM: {len(gisaid_0 & gadm_0)}')
        print(f'Countries in GISAID not captured in GADM: {gisaid_0 - gadm_0}')
        missing_countries = sorted(gisaid_0 - gadm_0)
        samples_missing_country = df[df['country'].isin(missing_countries)]
        print(f'Number of samples missing country-level geo-information: {samples_missing_country.shape[0]}')
        print(f'Admin1 standardization...')
        df.loc[df['division'].isna(), 'division'] = 'None'
        df['division_normed'] = df['division'].copy()
        df.loc[df['division_normed']=='USA', 'division_normed'] = 'United States'
        df.loc[df['division_normed'].str.contains('Georgia /'), 'division_normed'] = 'Georgia'
        df.loc[df['division_normed'].str.contains('Antwerp'), 'division_normed'] = 'Vlaanderen'
        df.loc[df['division_normed'].str.contains('Andalu'), 'division_normed'] = 'Andalucía'
        df.loc[df['division_normed'].str.contains('Cairo'), 'division_normed'] = 'Al Qahirah'
        df.loc[df['division_normed'].str.contains('Northern territory'), 'division_normed'] = 'Northern Territory'
        df.loc[df['division_normed'].str.contains('Fayoum'), 'division_normed'] = 'Al Fayyum'
        df.loc[df['division_normed'].str.contains('Musca'), 'division_normed'] = 'Muscat'
        df.loc[df['division_normed'].str.contains('Kalyoubia'), 'division_normed'] = 'Al Qalyubiyah'
        df.loc[df['division_normed'].str.contains('Buraymi'), 'division_normed'] = 'Al Buraymi'
        df.loc[df['division_normed'].str.contains('Buraimi'), 'division_normed'] = 'Al Buraymi'
        df.loc[df['division_normed'].str.contains('Dakhiliyah'), 'division_normed'] = 'Ad Dakhliyah'
        df.loc[df['division_normed'].str.contains('Dhahirah'), 'division_normed'] = 'Al Dhahira'
        df.loc[df['division_normed'].str.contains('North Batinah'), 'division_normed'] = 'Al Batinah North'
        df.loc[df['division_normed'].str.contains('South Batinah'), 'division_normed'] = 'Al Batinah South'
        df.loc[df['division_normed'].str.contains('North Sharqiyah'), 'division_normed'] = 'Ash Sharqiyah North'
        df.loc[df['division_normed'].str.contains('Wuhan'), 'division_normed'] = 'Hubei'
        df.loc[df['division_normed'].str.contains('Quebec'), 'division_normed'] = 'Québec'
        df.loc[df['division_normed'].str.contains('Toronto'), 'division_normed'] = 'Ontario'
        df.loc[df['division_normed'].str.contains('Coahuila de Zaragoza'), 'division_normed'] = 'Coahuila'
        df.loc[df['division_normed'].str.contains('Mexico City'), 'division_normed'] = 'México'
        df.loc[df['division_normed'].str.contains('Michoacan'), 'division_normed'] = 'Michoacán'
        df.loc[df['division_normed'].str.contains('Nuevo Leon'), 'division_normed'] = 'Nuevo León'
        df.loc[df['division_normed'].str.contains('Queretaro'), 'division_normed'] = 'Querétaro'
        df.loc[df['division_normed'].str.contains('SanLuisPotosi'), 'division_normed'] = 'San Luis Potosí'
        df.loc[df['division_normed'].str.contains('San Luis Potosi'), 'division_normed'] = 'San Luis Potosí'
        df.loc[df['division_normed'].str.contains('State of Mexico'), 'division_normed'] = 'México'
        df.loc[df['division_normed'].str.contains('Yucatan'), 'division_normed'] = 'Yucatán'
        df.loc[df['division_normed'].str.contains('Bethlehem'), 'division_normed'] = 'West Bank'
        df.loc[df['division_normed'].str.contains('Hebron'), 'division_normed'] = 'West Bank'
        df.loc[df['division_normed'].str.contains('Jenin'), 'division_normed'] = 'West Bank'
        df.loc[df['division_normed'].str.contains('Jericho'), 'division_normed'] = 'West Bank'
        df.loc[df['division_normed'].str.contains('Ramallah'), 'division_normed'] = 'West Bank'
        df.loc[df['division_normed'].str.contains('Tulkarem'), 'division_normed'] = 'West Bank'
        df.loc[df['division_normed'].str.contains('Nablus'), 'division_normed'] = 'West Bank'
        df.loc[df['division_normed'].str.contains('Sharja'), 'division_normed'] = 'Sharjah'
        df.loc[df['division_normed'].str.contains('Copenhagen'), 'division_normed'] = 'Hovedstaden'
        df.loc[df['division_normed'].str.contains('Sjaelland'), 'division_normed'] = 'Sjælland'
        df.loc[df['division_normed'].str.contains('Cape Town'), 'division_normed'] = 'Western Cape'
        df.loc[df['division_normed'].str.contains('Western Cape'), 'division_normed'] = 'Western Cape'
        df.loc[df['division_normed'].str.contains('Amapa'), 'division_normed'] = 'Amapá'
        df.loc[df['division_normed'].str.contains('Ceara'), 'division_normed'] = 'Ceará'
        df.loc[df['division_normed'].str.contains('Goias'), 'division_normed'] = 'Goiás'
        df.loc[df['division_normed'].str.contains('Maranhao'), 'division_normed'] = 'Maranhão'
        df.loc[df['division_normed'].str.contains('Paraiba'), 'division_normed'] = 'Paraíba'
        df.loc[df['division_normed'].str.contains('Parana'), 'division_normed'] = 'Paraná'
        df.loc[df['division_normed'].str.contains('Piaui'), 'division_normed'] = 'Piauí'
        df.loc[df['division_normed'].str.contains('Sao Paulo'), 'division_normed'] = 'São Paulo'
        df.loc[df['division_normed'].str.contains('Aragon'), 'division_normed'] = 'Aragón'
        df.loc[df['division_normed'].str.contains('Asturias'), 'division_normed'] = 'Principado de Asturias'
        df.loc[df['division_normed'].str.contains('Balear Islands'), 'division_normed'] = 'Islas Baleadf'
        df.loc[df['division_normed'].str.contains('Balear_Islands'), 'division_normed'] = 'Islas Baleadf'
        df.loc[df['division_normed'].str.contains('Illes Balears'), 'division_normed'] = 'Islas Baleadf'
        df.loc[df['division_normed'].str.contains('Canary Islands'), 'division_normed'] = 'Canaries'
        df.loc[df['division_normed'].str.contains('Canary_Islands'), 'division_normed'] = 'Canaries'
        df.loc[df['division_normed'].str.contains('Castilla La Mancha'), 'division_normed'] = 'Castilla-La Mancha'
        df.loc[df['division_normed'].str.contains('Castilla la Mancha'), 'division_normed'] = 'Castilla-La Mancha'
        df.loc[df['division_normed'].str.contains('Castilla y Leon'), 'division_normed'] = 'Castilla y León'
        df.loc[df['division_normed'].str.contains('Ceuta'), 'division_normed'] = 'Ceuta y Melilla'
        df.loc[df['division_normed'].str.contains('Melilla'), 'division_normed'] = 'Ceuta y Melilla'
        df.loc[df['division_normed'].str.contains('Comunitat Valenciana'), 'division_normed'] = 'Comunidad Valenciana'
        df.loc[df['division_normed'].str.contains('Comunitat_Valenciana'), 'division_normed'] = 'Comunidad Valenciana'
        df.loc[df['division_normed'].str.contains('La_Rioja'), 'division_normed'] = 'La Rioja'
        df.loc[df['division_normed'].str.contains('Madrid'), 'division_normed'] = 'Comunidad de Madrid'
        df.loc[df['division_normed'].str.contains('Murcia'), 'division_normed'] = 'Región de Murcia'
        df.loc[df['division_normed'].str.contains('Navarra'), 'division_normed'] = 'Comunidad Foral de Navarra'
        df.loc[df['division_normed'].str.contains('Catalunya'), 'division_normed'] = 'Cataluña'
        df.loc[df['division_normed'].str.contains('Catalonia'), 'division_normed'] = 'Cataluña'
        # Germany
        df.loc[df['division_normed'].str.contains('Baden-Wuerttemberg'), 'division_normed'] = 'Baden-Württemberg'
        df.loc[df['division_normed'].str.contains('Baden-Wurttemberg'), 'division_normed'] = 'Baden-Württemberg'
        df.loc[df['division_normed'].str.contains('Bavaria'), 'division_normed'] = 'Bayern'
        df.loc[df['division_normed'].str.contains('Hesse'), 'division_normed'] = 'Hessen'
        df.loc[df['division_normed'].str.contains('Lower Saxony'), 'division_normed'] = 'Niedersachsen'
        df.loc[df['division_normed'].str.contains('Mecklenburg-Western Pomerania'), 'division_normed'] = 'Mecklenburg-Vorpommern'
        df.loc[df['division_normed'].str.contains('Rhineland-Palatinate'), 'division_normed'] = 'Rheinland-Pfalz'
        df.loc[(df['division_normed'].str.contains('Saxony-Anhalt'))
             & (df['country_normed'].str.contains('Germany')), 'division_normed'] = 'Sachsen-Anhalt'
        df.loc[(df['division_normed'].str.contains('Saxony'))
             & (df['country_normed'].str.contains('Germany')), 'division_normed'] = 'Sachsen'
        df.loc[df['division_normed'].str.contains('North Rhine-Westphalia'), 'division_normed'] = 'Nordrhein-Westfalen'
        df.loc[df['division_normed'].str.contains('Thuringia'), 'division_normed'] = 'Thüringen'
        # South Africa
        df.loc[(df['country_normed'].str.contains('South Africa'))
                & (df['division_normed'].str.contains('KwaZulu Natal')), 
               'division_normed'] = 'KwaZulu-Natal'
        df.loc[(df['country_normed'].str.contains('South Africa'))
             & (df['division_normed'].str.contains('Northern Cape Province')), 
               'division_normed'] = 'Northern Cape'
        # Austria
        df.loc[(df['country_normed']=='Austria') 
             & (df['division_normed']=='Tyrol'), 'division_normed'] = 'Tirol'
        # Switzerland
        df.loc[df['division_normed'].str.contains('Argovie'), 'division_normed'] = 'Aargau'
        df.loc[df['division_normed'].str.contains('Geneva'), 'division_normed'] = 'Genève'
        df.loc[df['division_normed'].str.contains('Grabunden'), 'division_normed'] = 'Graubünden'
        df.loc[df['division_normed'].str.contains('Luzern'), 'division_normed'] = 'Lucerne' 
        df.loc[df['division_normed'].str.contains('Neuchatel'), 'division_normed'] = 'Neuchâtel' 
        df.loc[df['division_normed'].str.contains('Neuenburg'), 'division_normed'] = 'Neuchâtel' 
        df.loc[df['division_normed'].str.contains('Obwald'), 'division_normed'] = 'Obwalden' 
        df.loc[df['division_normed'].str.contains('Saint-Gall'), 'division_normed'] = 'Sankt Gallen' 
        df.loc[df['division_normed'].str.contains('St Gallen'), 'division_normed'] = 'Sankt Gallen' 
        df.loc[df['division_normed'].str.contains('St. Gallen'), 'division_normed'] = 'Sankt Gallen' 
        df.loc[df['division_normed'].str.contains('Schaffhouse'), 'division_normed'] = 'Schaffhausen' 
        df.loc[df['division_normed'].str.contains('Turgovia'), 'division_normed'] = 'Thurgau' 
        df.loc[df['division_normed'].str.contains('VALAIS'), 'division_normed'] = 'Valais' 
        df.loc[df['division_normed'].str.contains('Waadt'), 'division_normed'] = 'Vaud' 
        df.loc[df['division_normed'].str.contains('Wallis'), 'division_normed'] = 'Valais' 
        df.loc[df['division_normed'].str.contains('Zaerich'), 'division_normed'] = 'Zürich' 
        df.loc[df['division_normed'].str.contains('Zoerich'), 'division_normed'] = 'Zürich' 
        df.loc[df['division_normed'].str.contains('Zurich'), 'division_normed'] = 'Zürich' 
        # Sweden
        df.loc[(df['country_normed']=='Sweden')
              &(df['division_normed']=='Sodermanland'), 'division_normed'] = 'Södermanland'
        df.loc[(df['country_normed']=='Sweden')
            &(df['division_normed']=='Gavleborg'), 'division_normed'] = 'Gävleborg'
        df.loc[(df['country_normed']=='Sweden')
            &(df['division_normed']=='Jamtland'), 'division_normed'] = 'Jämtland'
        df.loc[(df['country_normed']=='Sweden')
            &(df['division_normed']=='Jonkoping'), 'division_normed'] = 'Jönköping'
        df.loc[(df['country_normed']=='Sweden')
            &(df['division_normed']=='Ostergotland'), 'division_normed'] = 'Östergötland'
        df.loc[(df['country_normed']=='Sweden')
            &(df['division_normed']=='Skane'), 'division_normed'] = 'Skåne'
        df.loc[(df['country_normed']=='Sweden')
            &(df['division_normed']=='Varmland'), 'division_normed'] = 'Värmland'
        df.loc[(df['country_normed']=='Sweden')
            &(df['division_normed']=='Vasterbotten'), 'division_normed'] = 'Västerbotten'
        df.loc[(df['country_normed']=='Sweden')
            &(df['division_normed']=='Vasternorrland'), 'division_normed'] = 'Västernorrland'
        df.loc[(df['country_normed']=='Sweden')
            &(df['division_normed']=='Vastmanland'), 'division_normed'] = 'Västmanland'
        df.loc[(df['country_normed']=='Sweden')
            &(df['division_normed']=='Vastra Gotaland'), 'division_normed'] = 'Västra Götaland'
        print(f'Admin2 standardization (U.S. only)')
        df.loc[df['location'].isna(), 'location'] = 'None'
        df['location_normed'] = df['location'].copy()
        houston_filter = (df['division']=='Texas') & (df['location']=='Houston')
        df.loc[houston_filter, 'location_normed'] = 'Harris County'
        for key, val in corrections.items():
            df.loc[:, 'location_normed'] = df['location_normed'].str.replace(key, val)
        df.loc[:, 'location_normed'] = df['location_normed'].str.replace('County', '').str.replace('county', '').str.replace(',', '')
        df.loc[:, 'location_normed'] = df['location_normed'].str.strip().apply(bs.check_state, args=(False,)).str.strip()
        df.loc[df['location_normed'].str.contains('Anchorage-Mat-Su'), 'location_normed'] = 'Anchorage'
        df.loc[df['location_normed'].str.contains('Anchorage-Mat Su'), 'location_normed'] = 'Anchorage'
        df.loc[df['location_normed'].str.contains('BRA'), 'location_normed'] = 'Brazos'
        df.loc[df['location_normed'].str.contains('BR'), 'location_normed'] = 'Brewster'
        df.loc[df['location_normed'].str.contains('Belgrade'), 'location_normed'] = 'Gallatin'
        df.loc[df['location_normed'].str.contains('Bozeman'), 'location_normed'] = 'Gallatin'
        df.loc[df['location_normed'].str.contains('Big Sky'), 'location_normed'] = 'Gallatin'
        df.loc[df['location_normed'].str.contains('Belton'), 'location_normed'] = 'Bell'
        df.loc[df['location_normed'].str.contains('Brentwood'), 'location_normed'] = 'Contra Costa'
        df.loc[df['location_normed'].str.contains('Chicago'), 'location_normed'] = 'Cook'
        df.loc[df['location_normed'].str.contains('Colombus'), 'location_normed'] = 'Franklin'
        df.loc[df['location_normed'].str.contains('DuBois'), 'location_normed'] = 'Fremont'
        df.loc[df['location_normed'].str.contains('DuPage'), 'location_normed'] = 'Dupage'
        df.loc[df['location_normed'].str.contains('Eau claire'), 'location_normed'] = 'Eau Claire'
        df.loc[df['location_normed'].str.contains('Ennis'), 'location_normed'] = 'Ellis'
        df.loc[df['location_normed'].str.contains('Fond Du Lac'), 'location_normed'] = 'Fond du Lac'
        df.loc[df['location_normed'].str.contains('Fond du lac'), 'location_normed'] = 'Fond du Lac'
        df.loc[df['location_normed'].str.contains('Fonddu Lac'), 'location_normed'] = 'Fond du Lac'
        df.loc[df['location_normed'].str.contains('Frisco'), 'location_normed'] = 'Collin'
        df.loc[df['location_normed'].str.contains('Hawai'), 'location_normed'] = 'Hawaii'
        df.loc[df['location_normed'].str.contains('Holland'), 'location_normed'] = 'Ottawa'
        df.loc[df['location_normed'].str.contains('Honolul'), 'location_normed'] = 'Honolulu'
        df.loc[df['location_normed'].str.contains('Indianapolis'), 'location_normed'] = 'Marion'
        df.loc[df['location_normed'].str.contains('Interior'), 'location_normed'] = 'Fairbanks North Star'
        df.loc[df['location_normed'].str.contains('Ithaca'), 'location_normed'] = 'Tompkins'
        df.loc[df['location_normed'].str.contains('Kaua'), 'location_normed'] = 'Kauai'
        df.loc[df['location_normed'].str.contains('Las Vegas'), 'location_normed'] = 'Clark'
        df.loc[df['location_normed'].str.contains('Mau'), 'location_normed'] = 'Hawaii'
        df.loc[df['location_normed'].str.contains('Mcculloch'), 'location_normed'] = 'McCulloch'
        df.loc[df['location_normed'].str.contains('Mchenry'), 'location_normed'] = 'McHenry'
        df.loc[df['location_normed'].str.contains('Mclennan'), 'location_normed'] = 'McLennan'
        df.loc[df['location_normed'].str.contains('Moris'), 'location_normed'] = 'Morris'
        df.loc[df['location_normed'].str.contains('New York'), 'location_normed'] = 'New York'
        df.loc[df['location_normed'].str.contains('New York City'), 'location_normed'] = 'New York'
        df.loc[df['location_normed'].str.contains('New Hyde Park'), 'location_normed'] = 'Nassau'
        df.loc[df['location_normed'].str.contains('New Orleans'), 'location_normed'] = 'Orleans'
        df.loc[df['location_normed'].str.contains('New Rochelle'), 'location_normed'] = 'Westchester'
        df.loc[df['location_normed'].str.contains('Northern'), 'location_normed'] = 'Fairbanks North Star'
        df.loc[df['location_normed'].str.contains('Omaha'), 'location_normed'] = 'Douglas'
        df.loc[df['location_normed'].str.contains('Ostego'), 'location_normed'] = 'Allegan'
        df.loc[df['location_normed'].str.contains('Phoenix'), 'location_normed'] = 'Maricopa'
        df.loc[df['location_normed'].str.contains('San Bernadino'), 'location_normed'] = 'San Bernardino'
        df.loc[df['location_normed'].str.contains('Seattle'), 'location_normed'] = 'King'
        df.loc[df['location_normed'].str.contains('St. Bernard'), 'location_normed'] = 'Saint Bernard'
        df.loc[df['location_normed'].str.contains('St. Clair'), 'location_normed'] = 'Saint Clair'
        df.loc[df['location_normed'].str.contains('St. Lawrence'), 'location_normed'] = 'Saint Lawrence'
        df.loc[df['location_normed'].str.contains('St. Louis'), 'location_normed'] = 'Saint Louis'
        df.loc[df['location_normed'].str.contains('St. Tammany'), 'location_normed'] = 'Saint Tammany'
        df.loc[df['location_normed'].str.contains('Staten Island'), 'location_normed'] = 'Richmond'
        df.loc[df['location_normed'].str.contains('Thurson'), 'location_normed'] = 'Thurston'
        df.loc[df['location_normed'].str.contains('Tucson'), 'location_normed'] = 'Pima'
        df.loc[df['location_normed'].str.contains('West Yellowstone'), 'location_normed'] = 'Gallatin'
        df.loc[df['location_normed'].str.contains('Adam'), 'location_normed'] = 'Adams'
        df.loc[df['location_normed'].str.contains('Alachu'), 'location_normed'] = 'Alachua'
        df.loc[df['location_normed'].str.contains('Du Bois'), 'location_normed'] = 'Dubois'
        df.loc[df['location_normed'].str.contains('DeSoto'), 'location_normed'] = 'Desoto'
        df.loc[df['location_normed'].str.contains('PdfID'), 'location_normed'] = 'Pdfidio'
        df.loc[df['location_normed'].str.contains('LaSalle'), 'location_normed'] = 'La Salle'
        df.loc[df['location_normed'].str.contains('CAMER'), 'location_normed'] = 'Cameron'
        df.loc[df['location_normed'].str.contains('CAST'), 'location_normed'] = 'Castro'
        df.loc[df['location_normed'].str.contains('CROS'), 'location_normed'] = 'Crosby'
        df.loc[df['location_normed'].str.contains('ECT'), 'location_normed'] = 'Ector'
        df.loc[df['location_normed'].str.contains('GALVEST'), 'location_normed'] = 'Galveston'
        df.loc[df['location_normed'].str.contains('JEFFERS'), 'location_normed'] = 'Jefferson'
        df.loc[df['location_normed'].str.contains('KAUFM'), 'location_normed'] = 'Kaufman'
        df.loc[df['location_normed'].str.contains('KLEBE'), 'location_normed'] = 'Kleberg'
        df.loc[df['location_normed'].str.contains('LAVA'), 'location_normed'] = 'Lavaca'
        df.loc[df['location_normed'].str.contains('MCLENN'), 'location_normed'] = 'Mclennan'
        df.loc[df['location_normed'].str.contains('St.Clair'), 'location_normed'] = 'Saint Clair'
        df.loc[df['location_normed'].str.contains('TARRA'), 'location_normed'] = 'Tarrant'
        df.loc[df['location_normed'].str.contains('WALL'), 'location_normed'] = 'Waller'
        df.loc[df['location_normed'].str.contains('WICHI'), 'location_normed'] = 'Wichita'
        # TODO: France, Russia, China, Israel, South Korea
        country = 'USA'
        if country:
            gisaid_2 = set(df[df['country']==country]['location_normed'].unique())
        else:
            gisaid_2 = set(df['location_normed'].unique())
        gadm_2 = set(gadm[(~gadm['NAME_2'].isna())&(gadm['NAME_0']=='United States')]['NAME_2'].unique())
        print(len(gisaid_2))
        print(len(gadm_2))
        print(len(gisaid_2&gadm_2))
        locs_missing = sorted(gisaid_2 - gadm_2)
        samples_missing_county = df.loc[(df['location_normed'].isin(locs_missing))&(df['country']=='USA')]
        print(f'Number of samples missing county-level geo-information (U.S. only): {samples_missing_county.shape[0]}')
        print(f'Metadata generated and saving to {meta_fp}')
        df['strain'] = df['strain'].str.replace('hCoV-19/', '').str.replace(' ', '')
        df['country'] = df['country'].astype(str)
        df['country_lower'] = df['country'].str.lower()
        df['country_normed'] = df['country_normed'].astype(str)
        df['country_normed_lower'] = df['country_normed'].str.lower()
        df['division'] = df['division'].astype(str)
        df['division_lower'] = df['division'].str.lower()
        df['division_normed'] = df['division_normed'].astype(str)
        df['division_normed_lower'] = df['division_normed'].str.lower()
        df['location'] = df['location'].astype(str)
        df['location_lower'] = df['location'].str.lower()
        df['location_normed'] = df['location_normed'].astype(str)
        df['location_normed_lower'] = df['location_normed'].str.lower()
        df.to_csv(meta_fp, sep='\t', index=False, compression='gzip')
        del df
        gc.collect();
        print(f'GISAID API feed has been downloaded and processed; ready for variant counting.')
    else:
        print(f"{meta_fp} already exists, skipping metadata generation step")
    info_df = bs.create_chunk_names(meta_fp, chunk_size)
    info_df.to_csv(info_fp, index=False)
    # del data
    # gc.collect();
    return 0


if __name__=="__main__":
    # COLLECTING USER PARAMETERS
    parser = argparse.ArgumentParser()
    parser.add_argument("-u", "--username",
                            type=str,
                            required=True,
                            help="Username to use for downloading API feed from GISAID")
    parser.add_argument("-p", "--password",
                            type=str,
                            required=True,
                            help="Password to use for downloading API feed from GISAID")
    parser.add_argument("-s", "--size",
                            type=int,
                            required=True,
                            help="Chunk size")
    parser.add_argument("-t", "--time",
                            type=str,
                            required=True,
                            help="Current datetime")
    args = parser.parse_args()
    username = args.username
    password = args.password
    chunk_size = args.size
    current_datetime = args.time
    result = download_process_data(username, password, chunk_size, current_datetime)
    assert result==0, "ERROR: downloading GISAID data incomplete. Please inspect log."
    # info_df = bs.create_chunk_names(meta_fp, chunk_size)


# with open('config.json', 'r') as f:
#     config = json.load(f)

# out_dir = Path(config['out_dir'])
# in_fp = out_dir/config['gisaid_feed']
# out_fp = out_dir/config['gisaid_fasta']
# meta_fp = out_dir/config['gisaid_meta']
# gadm_fp = config['gadm']
# is_test = config['feed_test']
# test_size = 100
# metacols = ['covv_virus_name', 'covsurver_prot_mutations', 'covv_location',
#              'covv_lineage', 'covv_collection_date', 'covv_accession_id',
#              'pangolin_lineages_version', 'covv_clade', 'covv_subm_date']
# corrections = bd.COUNTY_CORRECTIONS



# # Download GISAID API feed
# if not Path.isfile(Path(in_fp)):
#     print(f"Downloading GISAID API feed...")
#     feed_cmd = f"curl -u {username}:{password} https://www.epicov.org/epi3/3p/scripps/export/provision.json.xz | xz -d -T8 > {in_fp}"
#     bs.run_command(feed_cmd)
#     print(f"GISAID API feed saved in {in_fp}")
# else:
#     print(f"{in_fp} already exists, skipping API feed download step")
# # load sequence data
# print(f"Loading API data...")
# if is_test:
#     data = []
#     with open(in_fp, 'r') as handle:
#         for i, line in enumerate(handle):
#             data.append(json.loads(line))
#             if i >= test_size:
#                 break
# else:
#     data = [json.loads(line) for line in open(in_fp, 'r')]
# print(f"Total number of sequences: {len(data)}")
# # generate fasta file containing all sequences
# if not Path.isfile(Path(out_fp)):
#     print(f"Converting to dict...")
#     regex = re.compile('[^a-zA-Z]')
#     seqs_dict = {sample['covv_virus_name'].replace('hCoV-19/', '').replace(' ', ''): 
#                 regex.sub('', sample['sequence'].replace('\n', '')) for sample in data}
#     print(f"Converting to FASTA...")
#     bs.dict2fasta(seqs_dict, out_fp)
#     print(f"FASTA output generated and saved in {out_fp}")
# else:
#     print(f"{out_fp} already exists, skipping fasta generation step")
# # generate tsv file containing processed metadata
# if not Path.isfile(Path(meta_fp)):
#     # load raw metadata into dataframe
#     df = pd.DataFrame(data, columns=metacols)
#     # TEST: all accession IDs are unique
#     assert df['covv_accession_id'].shape[0]==df['covv_accession_id'].unique().shape[0], f'ERROR: gisaid accession ids not unique'
#     num_ids = df['covv_accession_id'].unique().shape[0]
#     print(f"Total number of sequences: {num_ids}")
#     print(f"Cleaning metadata")
#     df.rename(columns={
#                     'covv_virus_name': 'strain', 
#                     'covv_location': 'location', 
#                     'covv_collection_date': 'date_collected',
#                     'covv_subm_date': 'date_submitted',
#                     'covv_clade': 'clade',
#                     'covv_lineage': 'pangolin_lineage',
#                     'pangolin_lineages_version': 'pangolin_version',
#                     'covv_accession_id': 'accession_id'
#                     }, inplace=True)
#     print(f"Loading GADM file containing geographic information...")
#     gadm = gpd.read_file(gadm_fp)
#     gadm_cols = [f'NAME_{i}' for i in range(5)]
#     gadm = gadm[gadm_cols]
#     print(f"Standardizing location information...")
#     res = pd.DataFrame(df['location'].str.split('/').tolist(), 
#                 columns=['region',
#                         'country', 
#                         'division', 
#                         'location', 
#                         'city', 
#                         'town'
#                         ])
#     df['country'] = res['country'].str.strip()
#     df['division'] = res['division'].str.strip()
#     df['location'] = res['location'].str.strip()
#     print(f"Admin0 standardization...")
#     df['country_normed'] = df['country'].copy()
#     df['country_normed'].fillna('None', inplace=True)
#     df.loc[df['country_normed']=='USA', 'country_normed'] = 'United States'
#     df.loc[df['country_normed'].str.contains('Congo'), 'country_normed'] = 'Republic of Congo'
#     df.loc[df['country_normed'].str.contains('Cote dIvoire'), 'country_normed'] = "Côte d'Ivoire"
#     df.loc[df['country_normed'].str.contains('North Macedonia'), 'country_normed'] = "Macedonia"
#     df.loc[df['country_normed'].str.contains('Curacao'), 'country_normed'] = "Curaçao"
#     df.loc[df['country_normed'].str.contains('Saint Martin'), 'country_normed'] = "Saint-Martin"
#     df.loc[df['country_normed'].str.contains('Trinidad'), 'country_normed'] = 'Trinidad and Tobago'
#     df.loc[df['country_normed'].str.contains('Czech republic'), 'country_normed'] = 'Czech Republic'
#     df.loc[df['country_normed'].str.contains('St Eustatius'), 'country_normed'] = 'Netherlands'
#     df.loc[df['country_normed'].str.contains('Saint Barthelemy'), 'country_normed'] = 'Saint-Barthélemy'
#     df.loc[df['country_normed'].str.contains('Palestine'), 'country_normed'] = "Palestina"
#     df.loc[df['country_normed'].str.contains("Germany /"), 'country_normed'] = "Germany"
#     df.loc[df['country_normed'].str.contains("France /Nouvelle-Aquitaine"), 'division'] = "Nouvelle-Aquitaine"
#     df.loc[df['country_normed']=="France /Nouvelle-Aquitaine", 'country_normed'] = "France"
#     df.loc[df['country_normed'].str.contains("France /Nouvelle-Aquitaine/ Limoges"), 'division'] = "Nouvelle-Aquitaine"
#     df.loc[df['country_normed'].str.contains("France /Nouvelle-Aquitaine/ Limoges"), 'location'] = "Limoges"
#     df.loc[df['country_normed']=="France /Nouvelle-Aquitaine/ Limoges", 'country_normed'] = "France"
#     df.loc[df['country_normed']=="Kenya /", 'country_normed'] = "Kenya"
#     df.loc[df['country_normed']=="Switzerland/ Schwyz", 'division'] = "Schwyz"
#     df.loc[df['country_normed']=="Switzerland/ Schwyz", 'country_normed'] = "Switzerland"
#     df.loc[df['country_normed']=="USA /Wisconsin", 'division'] = "Wisconsin"
#     df.loc[df['country_normed']=="USA /Wisconsin", 'country_normed'] = "United States"
#     df.loc[df['country_normed']=="Jonavos apskritis", 'country_normed'] = "Lithuania"
#     df.loc[df['country_normed']=="Thailand /Singburi", 'division'] = "Singburi"
#     df.loc[df['country_normed']=="Thailand /Singburi", 'country_normed'] = "Thailand"
#     df.loc[df['country_normed']=="Norway /", 'country_normed'] = "Norway"
#     df.loc[df['country_normed']=="Morocoo", 'country_normed'] = "Morocco"
#     gisaid_0 = set(df['country_normed'].unique())
#     gadm_0 = set(gadm['NAME_0'].unique())
#     print(len(gisaid_0))
#     print(len(gadm_0))
#     print(f'Number of countries captured in GADM: {len(gisaid_0 & gadm_0)}')
#     print(f'Countries in GISAID not captured in GADM: {gisaid_0 - gadm_0}')
#     missing_countries = sorted(gisaid_0 - gadm_0)
#     samples_missing_country = df[df['country'].isin(missing_countries)]
#     print(f'Number of samples missing country-level geo-information: {samples_missing_country.shape[0]}')
#     print(f'Admin1 standardization...')
#     df.loc[df['division'].isna(), 'division'] = 'None'
#     df['division_normed'] = df['division'].copy()
#     df.loc[df['division_normed']=='USA', 'division_normed'] = 'United States'
#     df.loc[df['division_normed'].str.contains('Georgia /'), 'division_normed'] = 'Georgia'
#     df.loc[df['division_normed'].str.contains('Antwerp'), 'division_normed'] = 'Vlaanderen'
#     df.loc[df['division_normed'].str.contains('Andalu'), 'division_normed'] = 'Andalucía'
#     df.loc[df['division_normed'].str.contains('Cairo'), 'division_normed'] = 'Al Qahirah'
#     df.loc[df['division_normed'].str.contains('Northern territory'), 'division_normed'] = 'Northern Territory'
#     df.loc[df['division_normed'].str.contains('Fayoum'), 'division_normed'] = 'Al Fayyum'
#     df.loc[df['division_normed'].str.contains('Musca'), 'division_normed'] = 'Muscat'
#     df.loc[df['division_normed'].str.contains('Kalyoubia'), 'division_normed'] = 'Al Qalyubiyah'
#     df.loc[df['division_normed'].str.contains('Buraymi'), 'division_normed'] = 'Al Buraymi'
#     df.loc[df['division_normed'].str.contains('Buraimi'), 'division_normed'] = 'Al Buraymi'
#     df.loc[df['division_normed'].str.contains('Dakhiliyah'), 'division_normed'] = 'Ad Dakhliyah'
#     df.loc[df['division_normed'].str.contains('Dhahirah'), 'division_normed'] = 'Al Dhahira'
#     df.loc[df['division_normed'].str.contains('North Batinah'), 'division_normed'] = 'Al Batinah North'
#     df.loc[df['division_normed'].str.contains('South Batinah'), 'division_normed'] = 'Al Batinah South'
#     df.loc[df['division_normed'].str.contains('North Sharqiyah'), 'division_normed'] = 'Ash Sharqiyah North'
#     df.loc[df['division_normed'].str.contains('Wuhan'), 'division_normed'] = 'Hubei'
#     df.loc[df['division_normed'].str.contains('Quebec'), 'division_normed'] = 'Québec'
#     df.loc[df['division_normed'].str.contains('Toronto'), 'division_normed'] = 'Ontario'
#     df.loc[df['division_normed'].str.contains('Coahuila de Zaragoza'), 'division_normed'] = 'Coahuila'
#     df.loc[df['division_normed'].str.contains('Mexico City'), 'division_normed'] = 'México'
#     df.loc[df['division_normed'].str.contains('Michoacan'), 'division_normed'] = 'Michoacán'
#     df.loc[df['division_normed'].str.contains('Nuevo Leon'), 'division_normed'] = 'Nuevo León'
#     df.loc[df['division_normed'].str.contains('Queretaro'), 'division_normed'] = 'Querétaro'
#     df.loc[df['division_normed'].str.contains('SanLuisPotosi'), 'division_normed'] = 'San Luis Potosí'
#     df.loc[df['division_normed'].str.contains('San Luis Potosi'), 'division_normed'] = 'San Luis Potosí'
#     df.loc[df['division_normed'].str.contains('State of Mexico'), 'division_normed'] = 'México'
#     df.loc[df['division_normed'].str.contains('Yucatan'), 'division_normed'] = 'Yucatán'
#     df.loc[df['division_normed'].str.contains('Bethlehem'), 'division_normed'] = 'West Bank'
#     df.loc[df['division_normed'].str.contains('Hebron'), 'division_normed'] = 'West Bank'
#     df.loc[df['division_normed'].str.contains('Jenin'), 'division_normed'] = 'West Bank'
#     df.loc[df['division_normed'].str.contains('Jericho'), 'division_normed'] = 'West Bank'
#     df.loc[df['division_normed'].str.contains('Ramallah'), 'division_normed'] = 'West Bank'
#     df.loc[df['division_normed'].str.contains('Tulkarem'), 'division_normed'] = 'West Bank'
#     df.loc[df['division_normed'].str.contains('Nablus'), 'division_normed'] = 'West Bank'
#     df.loc[df['division_normed'].str.contains('Sharja'), 'division_normed'] = 'Sharjah'
#     df.loc[df['division_normed'].str.contains('Copenhagen'), 'division_normed'] = 'Hovedstaden'
#     df.loc[df['division_normed'].str.contains('Sjaelland'), 'division_normed'] = 'Sjælland'
#     df.loc[df['division_normed'].str.contains('Cape Town'), 'division_normed'] = 'Western Cape'
#     df.loc[df['division_normed'].str.contains('Western Cape'), 'division_normed'] = 'Western Cape'
#     df.loc[df['division_normed'].str.contains('Amapa'), 'division_normed'] = 'Amapá'
#     df.loc[df['division_normed'].str.contains('Ceara'), 'division_normed'] = 'Ceará'
#     df.loc[df['division_normed'].str.contains('Goias'), 'division_normed'] = 'Goiás'
#     df.loc[df['division_normed'].str.contains('Maranhao'), 'division_normed'] = 'Maranhão'
#     df.loc[df['division_normed'].str.contains('Paraiba'), 'division_normed'] = 'Paraíba'
#     df.loc[df['division_normed'].str.contains('Parana'), 'division_normed'] = 'Paraná'
#     df.loc[df['division_normed'].str.contains('Piaui'), 'division_normed'] = 'Piauí'
#     df.loc[df['division_normed'].str.contains('Sao Paulo'), 'division_normed'] = 'São Paulo'
#     df.loc[df['division_normed'].str.contains('Aragon'), 'division_normed'] = 'Aragón'
#     df.loc[df['division_normed'].str.contains('Asturias'), 'division_normed'] = 'Principado de Asturias'
#     df.loc[df['division_normed'].str.contains('Balear Islands'), 'division_normed'] = 'Islas Baleadf'
#     df.loc[df['division_normed'].str.contains('Balear_Islands'), 'division_normed'] = 'Islas Baleadf'
#     df.loc[df['division_normed'].str.contains('Illes Balears'), 'division_normed'] = 'Islas Baleadf'
#     df.loc[df['division_normed'].str.contains('Canary Islands'), 'division_normed'] = 'Canaries'
#     df.loc[df['division_normed'].str.contains('Canary_Islands'), 'division_normed'] = 'Canaries'
#     df.loc[df['division_normed'].str.contains('Castilla La Mancha'), 'division_normed'] = 'Castilla-La Mancha'
#     df.loc[df['division_normed'].str.contains('Castilla la Mancha'), 'division_normed'] = 'Castilla-La Mancha'
#     df.loc[df['division_normed'].str.contains('Castilla y Leon'), 'division_normed'] = 'Castilla y León'
#     df.loc[df['division_normed'].str.contains('Ceuta'), 'division_normed'] = 'Ceuta y Melilla'
#     df.loc[df['division_normed'].str.contains('Melilla'), 'division_normed'] = 'Ceuta y Melilla'
#     df.loc[df['division_normed'].str.contains('Comunitat Valenciana'), 'division_normed'] = 'Comunidad Valenciana'
#     df.loc[df['division_normed'].str.contains('Comunitat_Valenciana'), 'division_normed'] = 'Comunidad Valenciana'
#     df.loc[df['division_normed'].str.contains('La_Rioja'), 'division_normed'] = 'La Rioja'
#     df.loc[df['division_normed'].str.contains('Madrid'), 'division_normed'] = 'Comunidad de Madrid'
#     df.loc[df['division_normed'].str.contains('Murcia'), 'division_normed'] = 'Región de Murcia'
#     df.loc[df['division_normed'].str.contains('Navarra'), 'division_normed'] = 'Comunidad Foral de Navarra'
#     df.loc[df['division_normed'].str.contains('Catalunya'), 'division_normed'] = 'Cataluña'
#     df.loc[df['division_normed'].str.contains('Catalonia'), 'division_normed'] = 'Cataluña'
#     df.loc[df['division_normed'].str.contains('Baden-Wuerttemberg'), 'division_normed'] = 'Baden-Württemberg'
#     df.loc[df['division_normed'].str.contains('Baden-Wurttemberg'), 'division_normed'] = 'Baden-Württemberg'
#     df.loc[df['division_normed'].str.contains('Bavaria'), 'division_normed'] = 'Bayern'
#     df.loc[df['division_normed'].str.contains('Hesse'), 'division_normed'] = 'Hessen'
#     df.loc[df['division_normed'].str.contains('Lower Saxony'), 'division_normed'] = 'Niedersachsen'
#     df.loc[df['division_normed'].str.contains('Mecklenburg-Western Pomerania'), 'division_normed'] = 'Mecklenburg-Vorpommern'
#     df.loc[df['division_normed'].str.contains('Rhineland-Palatinate'), 'division_normed'] = 'Rheinland-Pfalz'
#     df.loc[df['division_normed'].str.contains('Saxony'), 'division_normed'] = 'Sachsen'
#     df.loc[df['division_normed'].str.contains('Saxony-Anhalt'), 'division_normed'] = 'Sachsen-Anhalt'
#     df.loc[df['division_normed'].str.contains('North Rhine-Westphalia'), 'division_normed'] = 'Nordrhein-Westfalen'
#     df.loc[df['division_normed'].str.contains('Thuringia'), 'division_normed'] = 'Thüringen'
#     print(f'Admin2 standardization (U.S. only)')
#     df.loc[df['location'].isna(), 'location'] = 'None'
#     df['location_normed'] = df['location'].copy()
#     for key, val in corrections.items():
#         df.loc[:, 'location_normed'] = df['location_normed'].str.replace(key, val)
#     df.loc[:, 'location_normed'] = df['location_normed'].str.replace('County', '').str.replace('county', '').str.replace(',', '')
#     df.loc[:, 'location_normed'] = df['location_normed'].str.strip().apply(bv.check_state, args=(False,)).str.strip()
#     df.loc[df['location_normed'].str.contains('Anchorage-Mat-Su'), 'location_normed'] = 'Anchorage'
#     df.loc[df['location_normed'].str.contains('Anchorage-Mat Su'), 'location_normed'] = 'Anchorage'
#     df.loc[df['location_normed'].str.contains('BRA'), 'location_normed'] = 'Brazos'
#     df.loc[df['location_normed'].str.contains('BR'), 'location_normed'] = 'Brewster'
#     df.loc[df['location_normed'].str.contains('Belgrade'), 'location_normed'] = 'Gallatin'
#     df.loc[df['location_normed'].str.contains('Bozeman'), 'location_normed'] = 'Gallatin'
#     df.loc[df['location_normed'].str.contains('Big Sky'), 'location_normed'] = 'Gallatin'
#     df.loc[df['location_normed'].str.contains('Belton'), 'location_normed'] = 'Bell'
#     df.loc[df['location_normed'].str.contains('Brentwood'), 'location_normed'] = 'Contra Costa'
#     df.loc[df['location_normed'].str.contains('Chicago'), 'location_normed'] = 'Cook'
#     df.loc[df['location_normed'].str.contains('Colombus'), 'location_normed'] = 'Franklin'
#     df.loc[df['location_normed'].str.contains('DuBois'), 'location_normed'] = 'Fremont'
#     df.loc[df['location_normed'].str.contains('DuPage'), 'location_normed'] = 'Dupage'
#     df.loc[df['location_normed'].str.contains('Eau claire'), 'location_normed'] = 'Eau Claire'
#     df.loc[df['location_normed'].str.contains('Ennis'), 'location_normed'] = 'Ellis'
#     df.loc[df['location_normed'].str.contains('Fond Du Lac'), 'location_normed'] = 'Fond du Lac'
#     df.loc[df['location_normed'].str.contains('Fond du lac'), 'location_normed'] = 'Fond du Lac'
#     df.loc[df['location_normed'].str.contains('Fonddu Lac'), 'location_normed'] = 'Fond du Lac'
#     df.loc[df['location_normed'].str.contains('Frisco'), 'location_normed'] = 'Collin'
#     df.loc[df['location_normed'].str.contains('Hawai'), 'location_normed'] = 'Hawaii'
#     df.loc[df['location_normed'].str.contains('Holland'), 'location_normed'] = 'Ottawa'
#     df.loc[df['location_normed'].str.contains('Honolul'), 'location_normed'] = 'Honolulu'
#     df.loc[df['location_normed'].str.contains('Indianapolis'), 'location_normed'] = 'Marion'
#     df.loc[df['location_normed'].str.contains('Interior'), 'location_normed'] = 'Fairbanks North Star'
#     df.loc[df['location_normed'].str.contains('Ithaca'), 'location_normed'] = 'Tompkins'
#     df.loc[df['location_normed'].str.contains('Kaua'), 'location_normed'] = 'Kauai'
#     df.loc[df['location_normed'].str.contains('Las Vegas'), 'location_normed'] = 'Clark'
#     df.loc[df['location_normed'].str.contains('Mau'), 'location_normed'] = 'Hawaii'
#     df.loc[df['location_normed'].str.contains('Mcculloch'), 'location_normed'] = 'McCulloch'
#     df.loc[df['location_normed'].str.contains('Mchenry'), 'location_normed'] = 'McHenry'
#     df.loc[df['location_normed'].str.contains('Mclennan'), 'location_normed'] = 'McLennan'
#     df.loc[df['location_normed'].str.contains('Moris'), 'location_normed'] = 'Morris'
#     df.loc[df['location_normed'].str.contains('New York'), 'location_normed'] = 'New York'
#     df.loc[df['location_normed'].str.contains('New York City'), 'location_normed'] = 'New York'
#     df.loc[df['location_normed'].str.contains('New Hyde Park'), 'location_normed'] = 'Nassau'
#     df.loc[df['location_normed'].str.contains('New Orleans'), 'location_normed'] = 'Orleans'
#     df.loc[df['location_normed'].str.contains('New Rochelle'), 'location_normed'] = 'Westchester'
#     df.loc[df['location_normed'].str.contains('Northern'), 'location_normed'] = 'Fairbanks North Star'
#     df.loc[df['location_normed'].str.contains('Omaha'), 'location_normed'] = 'Douglas'
#     df.loc[df['location_normed'].str.contains('Ostego'), 'location_normed'] = 'Allegan'
#     df.loc[df['location_normed'].str.contains('Phoenix'), 'location_normed'] = 'Maricopa'
#     df.loc[df['location_normed'].str.contains('San Bernadino'), 'location_normed'] = 'San Bernardino'
#     df.loc[df['location_normed'].str.contains('Seattle'), 'location_normed'] = 'King'
#     df.loc[df['location_normed'].str.contains('St. Bernard'), 'location_normed'] = 'Saint Bernard'
#     df.loc[df['location_normed'].str.contains('St. Clair'), 'location_normed'] = 'Saint Clair'
#     df.loc[df['location_normed'].str.contains('St. Lawrence'), 'location_normed'] = 'Saint Lawrence'
#     df.loc[df['location_normed'].str.contains('St. Louis'), 'location_normed'] = 'Saint Louis'
#     df.loc[df['location_normed'].str.contains('St. Tammany'), 'location_normed'] = 'Saint Tammany'
#     df.loc[df['location_normed'].str.contains('Staten Island'), 'location_normed'] = 'Richmond'
#     df.loc[df['location_normed'].str.contains('Thurson'), 'location_normed'] = 'Thurston'
#     df.loc[df['location_normed'].str.contains('Tucson'), 'location_normed'] = 'Pima'
#     df.loc[df['location_normed'].str.contains('West Yellowstone'), 'location_normed'] = 'Gallatin'
#     df.loc[df['location_normed'].str.contains('Adam'), 'location_normed'] = 'Adams'
#     df.loc[df['location_normed'].str.contains('Alachu'), 'location_normed'] = 'Alachua'
#     df.loc[df['location_normed'].str.contains('Du Bois'), 'location_normed'] = 'Dubois'
#     df.loc[df['location_normed'].str.contains('DeSoto'), 'location_normed'] = 'Desoto'
#     df.loc[df['location_normed'].str.contains('PdfID'), 'location_normed'] = 'Pdfidio'
#     df.loc[df['location_normed'].str.contains('LaSalle'), 'location_normed'] = 'La Salle'
#     df.loc[df['location_normed'].str.contains('CAMER'), 'location_normed'] = 'Cameron'
#     df.loc[df['location_normed'].str.contains('CAST'), 'location_normed'] = 'Castro'
#     df.loc[df['location_normed'].str.contains('CROS'), 'location_normed'] = 'Crosby'
#     df.loc[df['location_normed'].str.contains('ECT'), 'location_normed'] = 'Ector'
#     df.loc[df['location_normed'].str.contains('GALVEST'), 'location_normed'] = 'Galveston'
#     df.loc[df['location_normed'].str.contains('JEFFERS'), 'location_normed'] = 'Jefferson'
#     df.loc[df['location_normed'].str.contains('KAUFM'), 'location_normed'] = 'Kaufman'
#     df.loc[df['location_normed'].str.contains('KLEBE'), 'location_normed'] = 'Kleberg'
#     df.loc[df['location_normed'].str.contains('LAVA'), 'location_normed'] = 'Lavaca'
#     df.loc[df['location_normed'].str.contains('MCLENN'), 'location_normed'] = 'Mclennan'
#     df.loc[df['location_normed'].str.contains('St.Clair'), 'location_normed'] = 'Saint Clair'
#     df.loc[df['location_normed'].str.contains('TARRA'), 'location_normed'] = 'Tarrant'
#     df.loc[df['location_normed'].str.contains('WALL'), 'location_normed'] = 'Waller'
#     df.loc[df['location_normed'].str.contains('WICHI'), 'location_normed'] = 'Wichita'
#     # TODO: France, Russia, China, Israel, South Korea
#     country = 'USA'
#     if country:
#         gisaid_2 = set(df[df['country']==country]['location_normed'].unique())
#     else:
#         gisaid_2 = set(df['location_normed'].unique())
#     gadm_2 = set(gadm[(~gadm['NAME_2'].isna())&(gadm['NAME_0']=='United States')]['NAME_2'].unique())
#     print(len(gisaid_2))
#     print(len(gadm_2))
#     print(len(gisaid_2&gadm_2))
#     locs_missing = sorted(gisaid_2 - gadm_2)
#     samples_missing_county = df.loc[(df['location_normed'].isin(locs_missing))&(df['country']=='USA')]
#     print(f'Number of samples missing county-level geo-information (U.S. only): {samples_missing_county.shape[0]}')
#     print(f'Metadata generated and saving to {meta_fp}')
#     df['strain'] = df['strain'].str.replace('hCoV-19/', '').str.replace(' ', '')
#     df['country'] = df['country'].astype(str)
#     df['country_lower'] = df['country'].str.lower()
#     df['country_normed'] = df['country_normed'].astype(str)
#     df['country_normed_lower'] = df['country_normed'].str.lower()
#     df['division'] = df['division'].astype(str)
#     df['division_lower'] = df['division'].str.lower()
#     df['division_normed'] = df['division_normed'].astype(str)
#     df['division_normed_lower'] = df['division_normed'].str.lower()
#     df['location'] = df['location'].astype(str)
#     df['location_lower'] = df['location'].str.lower()
#     df['location_normed'] = df['location_normed'].astype(str)
#     df['location_normed_lower'] = df['location_normed'].str.lower()
#     df.to_csv(meta_fp, sep='\t', index=False, compression='gzip')
#     print(f'GISAID API feed has been downloaded and processed; ready for variant counting.')
# else:
#     print(f"{meta_fp} already exists, skipping metadata generation step")