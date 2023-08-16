from dask.distributed import Client, LocalCluster
cluster = LocalCluster(n_workers=72, threads_per_worker=2, processes=True, memory_limit=0)
# cluster.scale(144)
client = Client(cluster)

import dask
import dask.dataframe as dd

from datetime import date, timedelta
import requests

who_cases = dd.read_csv('https://covid19.who.int/WHO-COVID-19-global-data.csv')
cdc_cases = dd.read_csv('https://data.cdc.gov/api/views/3nnm-4jni/rows.csv?accessType=DOWNLOAD', blocksize=None, assume_missing=True)
cdc_deaths = dd.read_csv('https://data.cdc.gov/api/views/ite7-j2w7/rows.csv?accessType=DOWNLOAD', blocksize=None, assume_missing=True)
jhu_cases = dd.read_csv('time_series_covid19_confirmed_US.csv')
fips = {
       "ALABAMA": 1,
       "ALASKA": 2,
       "ARIZONA": 4,
       "ARKANSAS": 5,
       "CALIFORNIA": 6,
       "COLORADO": 8,
       "CONNECTICUT": 9,
       "DELAWARE": 10,
       "DISTRICT OF COLUMBIA": 11,
       "FLORIDA": 12,
       "GEORGIA": 13,
       "HAWAII": 15,
       "IDAHO": 16,
       "ILLINOIS": 17,
       "INDIANA": 18,
       "IOWA": 19,
       "KANSAS": 20,
       "KENTUCKY": 21,
       "LOUISIANA": 22,
       "MAINE": 23,
       "MARYLAND": 24,
       "MASSACHUSETTS": 25,
       "MICHIGAN": 26,
       "MINNESOTA": 27,
       "MISSISSIPPI": 28,
       "MISSOURI": 29,
       "MONTANA": 30,
       "NEBRASKA": 31,
       "NEVADA": 32,
       "NEW HAMPSHIRE": 33,
       "NEW JERSEY": 34,
       "NEW MEXICO": 35,
       "NEW YORK": 36,
       "NORTH CAROLINA": 37,
       "NORTH DAKOTA": 38,
       "OHIO": 39,
       "OKLAHOMA": 40,
       "OREGON": 41,
       "PENNSYLVANIA": 42,
       "RHODE ISLAND": 44,
       "SOUTH CAROLINA": 45,
       "SOUTH DAKOTA": 46,
       "TENNESSEE": 47,
       "TEXAS": 48,
       "UTAH": 49,
       "VERMONT": 50,
       "VIRGINIA": 51,
       "WASHINGTON": 53,
       "WEST VIRGINIA": 54,
       "WISCONSIN": 55,
       "WYOMING": 56,
       "PUERTO RICO": 72,
       "UNITED STATES VIRGIN ISLANDS": 78,
       "AMERICAN SAMOA": 60,
       "COMMONWEALTH OF THE NORTHERN MARIANA ISLANDS": 69,
       "GUAM": 66
}

def group_dates_who(x):
    b = x.iloc[0].copy().drop(['New_cases', 'Cumulative_cases', 'New_deaths', 'Cumulative_deaths', 'Date_reported'])
    b['dates_list'] = x.loc[:, 'Date_reported'].values
    b['cases'] = x.loc[:, 'New_cases'].fillna(0).values
    b['losses'] = x.loc[:, 'New_deaths'].fillna(0).values
    return b
c=who_cases.groupby('Country_code').apply(group_dates_who).compute()
def group_dates_cdc(x):
    b = x.iloc[0].copy().drop(['date', 'state_fips', 'county_fips', 'covid_cases_per_100k', 'date_updated', 'Week-Ending Date', 'FIPS Code', 'COVID-19 Deaths'])
    b['dates_list'] = x.loc[:, 'date'].values
    b['cases'] = x.loc[:, 'cases'].fillna(0).values
    b['losses'] = x.loc[:, 'COVID-19 Deaths'].fillna(0).values
    return b
d=cdc_cases[['county', 'state', 'county_fips', 'county_population', 'covid_cases_per_100k', 'date_updated']]
d['state_fips'] = d['state'].str.upper().apply(lambda x: fips[x], meta='int')
d['FIPS'] = d['county_fips'].astype('int')
d['cases'] = (d['covid_cases_per_100k'] * d['county_population'] / 1e5).round(0)
d['date'] = dd.to_datetime(d['date_updated'])
#d['FD'] = d['FIPS'].astype(str) + d['date'].astype(str)
e=cdc_deaths[['Week-Ending Date', 'FIPS Code', 'COVID-19 Deaths']]
e['FIPS'] = e['FIPS Code'].astype('int')
e['date'] = dd.to_datetime(e['Week-Ending Date'])
#e['FD'] = e['FIPS'].astype(str) + e['date'].astype(str)
f = dd.merge_asof(d.sort_values('date'), e.sort_values('date'), by='FIPS', on='date', direction='nearest')
# dd.merge(d, e, on=['FD'], how='left')
f=f.groupby('FIPS').apply(group_dates_cdc).compute()

global_epi = c

# global geo matcher
# global_epi["UID"] = global_epi.index.values
global_epi["iso2"] = global_epi["Country_code"]
global_epi = global_epi.set_index('iso2', drop=True)
global_epi.merge(ll[["iso2", "iso3", "Population", "FIPS", "Country_Region", "Lat", "Long_"]], on="iso2", how="inner")

#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import json
import sys

workdir = rundir = '.'

# prep geo data
print("prepping geo data", file=sys.stderr)
loc_lookup = pd.read_csv(workdir+"/UID_ISO_FIPS_LookUp_Table.csv")
countries_geo = pd.read_csv(rundir+"/countries_geo.tsv", sep="\t", header=None, names=["ID", "name", "geometry"])
countries_geo.loc[:, "geometry"] = countries_geo["geometry"].apply(json.loads)
wb_region_table = {"AFG":"South Asia","ALB":"Europe & Central Asia","DZA":"Middle East & North Africa","ASM":"East Asia & Pacific","AND":"Europe & Central Asia","AGO":"Sub-Saharan Africa","ATG":"Latin America & Caribbean","ARG":"Latin America & Caribbean","ARM":"Europe & Central Asia","ABW":"Latin America & Caribbean","AUS":"East Asia & Pacific","AUT":"Europe & Central Asia","AZE":"Europe & Central Asia","BHS":"Latin America & Caribbean","BHR":"Middle East & North Africa","BGD":"South Asia","BRB":"Latin America & Caribbean","BLR":"Europe & Central Asia","BEL":"Europe & Central Asia","BLZ":"Latin America & Caribbean","BEN":"Sub-Saharan Africa","BMU":"North America","BTN":"South Asia","BOL":"Latin America & Caribbean","BIH":"Europe & Central Asia","BWA":"Sub-Saharan Africa","BRA":"Latin America & Caribbean","VGB":"Latin America & Caribbean","BRN":"East Asia & Pacific","BGR":"Europe & Central Asia","BFA":"Sub-Saharan Africa","BDI":"Sub-Saharan Africa","CPV":"Sub-Saharan Africa","KHM":"East Asia & Pacific","CMR":"Sub-Saharan Africa","CAN":"North America","CYM":"Latin America & Caribbean","CAF":"Sub-Saharan Africa","TCD":"Sub-Saharan Africa","CHI":"Europe & Central Asia","CHL":"Latin America & Caribbean","CHN":"East Asia & Pacific","COL":"Latin America & Caribbean","COM":"Sub-Saharan Africa","COD":"Sub-Saharan Africa","COG":"Sub-Saharan Africa","CRI":"Latin America & Caribbean","CIV":"Sub-Saharan Africa","HRV":"Europe & Central Asia","CUB":"Latin America & Caribbean","CUW":"Latin America & Caribbean","CYP":"Europe & Central Asia","CZE":"Europe & Central Asia","DNK":"Europe & Central Asia","DJI":"Middle East & North Africa","DMA":"Latin America & Caribbean","DOM":"Latin America & Caribbean","ECU":"Latin America & Caribbean","EGY":"Middle East & North Africa","SLV":"Latin America & Caribbean","GNQ":"Sub-Saharan Africa","ERI":"Sub-Saharan Africa","EST":"Europe & Central Asia","SWZ":"Sub-Saharan Africa","ETH":"Sub-Saharan Africa","FRO":"Europe & Central Asia","FJI":"East Asia & Pacific","FIN":"Europe & Central Asia","FRA":"Europe & Central Asia","PYF":"East Asia & Pacific","GAB":"Sub-Saharan Africa","GMB":"Sub-Saharan Africa","GEO":"Europe & Central Asia","DEU":"Europe & Central Asia","GHA":"Sub-Saharan Africa","GIB":"Europe & Central Asia","GRC":"Europe & Central Asia","GRL":"Europe & Central Asia","GRD":"Latin America & Caribbean","GUM":"East Asia & Pacific","GTM":"Latin America & Caribbean","GIN":"Sub-Saharan Africa","GNB":"Sub-Saharan Africa","GUY":"Latin America & Caribbean","HTI":"Latin America & Caribbean","HND":"Latin America & Caribbean","HKG":"East Asia & Pacific","HUN":"Europe & Central Asia","ISL":"Europe & Central Asia","IND":"South Asia","IDN":"East Asia & Pacific","IRN":"Middle East & North Africa","IRQ":"Middle East & North Africa","IRL":"Europe & Central Asia","IMN":"Europe & Central Asia","ISR":"Middle East & North Africa","ITA":"Europe & Central Asia","JAM":"Latin America & Caribbean","JPN":"East Asia & Pacific","JOR":"Middle East & North Africa","KAZ":"Europe & Central Asia","KEN":"Sub-Saharan Africa","KIR":"East Asia & Pacific","PRK":"East Asia & Pacific","KOR":"East Asia & Pacific","XKX":"Europe & Central Asia","KWT":"Middle East & North Africa","KGZ":"Europe & Central Asia","LAO":"East Asia & Pacific","LVA":"Europe & Central Asia","LBN":"Middle East & North Africa","LSO":"Sub-Saharan Africa","LBR":"Sub-Saharan Africa","LBY":"Middle East & North Africa","LIE":"Europe & Central Asia","LTU":"Europe & Central Asia","LUX":"Europe & Central Asia","MAC":"East Asia & Pacific","MDG":"Sub-Saharan Africa","MWI":"Sub-Saharan Africa","MYS":"East Asia & Pacific","MDV":"South Asia","MLI":"Sub-Saharan Africa","MLT":"Middle East & North Africa","MHL":"East Asia & Pacific","MRT":"Sub-Saharan Africa","MUS":"Sub-Saharan Africa","MEX":"Latin America & Caribbean","FSM":"East Asia & Pacific","MDA":"Europe & Central Asia","MCO":"Europe & Central Asia","MNG":"East Asia & Pacific","MNE":"Europe & Central Asia","MAR":"Middle East & North Africa","MOZ":"Sub-Saharan Africa","MMR":"East Asia & Pacific","NAM":"Sub-Saharan Africa","NRU":
iso3_name_table = dict(zip(countries_geo["ID"], countries_geo["name"]))

# download US geo data
states_geo = pd.read_csv(rundir + "/states_geo.tsv", sep="\t", header=None, names=["ID", "state", "geometry"])
states_geo.loc[:, "geometry"] = states_geo["geometry"].apply(json.loads)
counties_geo = pd.read_csv(rundir + "/counties_geo.tsv", sep="\t", header=None, names=["ID", "state", "geometry"])
counties_geo.loc[:, "geometry"] = counties_geo["geometry"].apply(json.loads)
metros_geo = pd.read_csv(rundir + "/metros_geo.tsv", sep="\t")
metros_geo.loc[:, "geometry"] = metros_geo["geometry"].apply(json.loads)
us_state_to_abbrev = {"Alabama": "AL", "Alaska": "AK", "Arizona": "AZ", "Arkansas": "AR", "California": "CA", "Colorado": "CO", "Connecticut": "CT", "Delaware": "DE", "Florida": "FL", "Georgia": "GA", "Hawaii": "HI", "Idaho": "ID", "Illinois": "IL", "Indiana": "IN", "Iowa": "IA", "Kansas": "KS", "Kentucky": "KY", "Louisiana": "LA", "Maine": "ME", "Maryland": "MD", "Massachusetts": "MA", "Michigan": "MI", "Minnesota": "MN", "Mississippi": "MS", "Missouri": "MO", "Montana": "MT", "Nebraska": "NE", "Nevada": "NV", "New Hampshire": "NH", "New Jersey": "NJ", "New Mexico": "NM", "New York": "NY", "North Carolina": "NC", "North Dakota": "ND", "Ohio": "OH", "Oklahoma": "OK", "Oregon": "OR", "Pennsylvania": "PA", "Rhode Island": "RI", "South Carolina": "SC", "South Dakota": "SD", "Tennessee": "TN", "Texas": "TX", "Utah": "UT", "Vermont": "VT", "Virginia": "VA", "Washington": "WA", "West Virginia": "WV", "Wisconsin": "WI", "Wyoming": "WY", "District of Columbia": "DC", "American Samoa": "AS", "Guam": "GU", "Northern Mariana Islands": "MP", "Puerto Rico": "PR", "United States Minor Outlying Islands": "UM", "U.S. Virgin Islands": "VI"}
tostate = lambda x: us_state_to_abbrev[x] if x in us_state_to_abbrev else x.replace(" ", "")
us_county_to_metro = pd.read_csv(rundir + "/county_to_metro.csv")

# define shared columns
metakeys = ["geometry", "name", "FIPS", "population", "iso3", "admin1", "admin2", "lat", "long", "cases", "dates_list", "losses", "ID", "country_name"]

# global data read-in & numpy conversion
print("reading in global epi data", file=sys.stderr)
# global_epi = pd.read_csv(workdir+"/time_series_covid19_confirmed_global.csv")
# dates_list_global = pd.to_datetime(list(global_epi.loc[:, "1/22/20":]), format="%m/%d/%y").strftime("%Y-%m-%d")
# global_cases = list(global_epi.loc[:, "1/22/20":].to_numpy())
# global_epi = global_epi.loc[:, :"Long"]
# global_epi["cases"] = global_cases
# global_epi["losses"] = list(pd.read_csv(workdir+"/time_series_covid19_deaths_global.csv").loc[:, "1/22/20":].to_numpy())

global_epi = c

# global geo matcher
# global_epi["UID"] = global_epi.index.values
global_epi["iso2"] = global_epi["Country_code"]
global_data = global_epi.merge(ll[["iso2", "iso3", "Population", "FIPS", "Country_Region", "Lat", "Long_"]], on="iso2", how="inner")
global_data = global_data.rename(columns={"Country_Region":"country_name","Lat":"lat","Long_":"long","Population":"population"})
# global_data.loc[:, "country_name"] = global_data["country_name"].str.replace("*","", regex=False)
global_data["name"] = global_data["country_name"]
# global_data.loc[~pd.isnull(global_data["admin1"]), "name"] = global_data["admin1"]
# global_data = pd.merge(global_data, countries_geo[["name", "ID"]].rename(columns={"name":"country_name"}), on="country_name", how="left")
global_data = global_data.merge(cg[["name", "ID", "geometry"]], left_on="iso3", right_on="ID", how="inner")
global_data["iso3"] = global_data["ID"]
global_data["name"] = global_data["name_x"]
# global_data["ID"] = None
global_data["admin1"] = None
global_data["admin2"] = None
# global_data["dates_list"] = [dates_list_global]*len(global_data)
global_data = global_data[metakeys]
#missing_filter = global_data.drop(columns=["admin1", "population"]).isna().any(axis=1)
#missing_provinces = global_data[missing_filter]
#global_data = global_data[~missing_filter]

# US data read-in & numpy conversion
print("reading in county-level epi data", file=sys.stderr)
# counties_epi = pd.read_csv(workdir+"/time_series_covid19_confirmed_US.csv")
# dates_list_counties = pd.to_datetime(list(counties_epi.loc[:, "1/22/20":]), format="%m/%d/%y").strftime("%Y-%m-%d")
# counties_cases = list(counties_epi.loc[:, "1/22/20":].to_numpy())
# counties_epi = counties_epi.loc[:, :"Long_"]
# counties_epi["cases"] = counties_cases
# counties_epi["losses"] = list(pd.read_csv(workdir+"/time_series_covid19_deaths_US.csv").loc[:, "1/22/20":].to_numpy())

# US geo matcher
counties_epi = f
counties_epi["ID"] = "USA_US-" + counties_epi["state"].apply(tostate) + "_" + \
    counties_epi[["FIPS", "county"]].apply(lambda x: "{:05d}".format(int(x["FIPS"])) if pd.notna(x["FIPS"]) else x["county"], axis=1)
county_data = pd.merge(counties_epi, counties_geo.drop(["state"], axis=1), on="ID", how="left")
county_data = pd.merge(county_data, loc_lookup[["UID", "Population", "FIPS", "Lat", "Long_"]], on="FIPS", how="left")
#print(counties_epi["FIPS"], file=sys.stderr)
county_data = county_data.rename(columns={"county":"admin2","state":"admin1","Lat":"lat","Long_":"long","Population":"population"})
# county_data["dates_list"] = [dates_list_counties]*len(county_data)
county_data["country_name"] = "United States of America"
county_data["name"] = county_data["admin2"]
county_data["iso3"] = "USA"
county_data = county_data[metakeys]
#missing_filter = county_data.drop(columns=["admin1", "population"]).isna().any(axis=1)
#missing_counties = county_data[missing_filter]
#county_data = county_data[~missing_filter]

# gather global and US data
print("gathering and patching", file=sys.stderr)
all_data = global_data.append(county_data)
all_data = all_data.loc[~pd.isnull(all_data["iso3"])]
all_data = all_data.fillna(np.nan).replace([np.nan], [None])
set_loc_id = lambda x: "{}_{}_{}".format(x["iso3"], x["admin1"], x["admin2"])
all_data["location_id"] = all_data.apply(set_loc_id, axis=1)
all_data.loc[all_data["iso3"] == "USA", "country_name"] = "United States of America"
all_data["admin_level"] = 2
all_data["FIPS"] = pd.to_numeric(all_data["FIPS"])
all_data.loc[pd.isna(all_data["admin2"]), "admin_level"] = 1
all_data.loc[pd.isna(all_data["admin1"]), "admin_level"] = 0

# generate aggregate entries
def sumdata(x):
    if len(x) < 1:
        return None
    #print(x, file=sys.stderr)
    y = x.iloc[0].copy()
    y["cases"] = np.sum(x["cases"].apply(lambda y: y[:min(x["cases"].apply(len))]).fillna(0))
    y["losses"] = np.sum(x["losses"].apply(lambda y: y[:min(x["losses"].apply(len))]).fillna(0))
    y["population"] = np.sum(x["population"].fillna(0))
    #print(len(pd.DataFrame([y])), file=sys.stderr)
    return pd.DataFrame([y])

summed_metros = all_data.loc[all_data["admin_level"] == 2].merge(us_county_to_metro, on="FIPS", how="left").groupby(["CBSA"]).apply(sumdata)
summed_metros = summed_metros[["CBSA", "admin1", "iso3", "dates_list", "lat", "long", "country_name", "cases", "losses", "population"]].reset_index(drop=True)
summed_metros = summed_metros.merge(metros_geo, on="CBSA", how="left")
summed_metros.loc[:, "CBSA"] = summed_metros["CBSA"].apply(int)
summed_metros["admin2"] = summed_metros["name"]
summed_metros["admin_level"] = 1.5
summed_metros["ID"] = "USA_US-" + summed_metros["admin1"].apply(tostate) + "-" + summed_metros["CBSA"].apply(str)
summed_metros["location_id"] = summed_metros["ID"]
all_data = all_data.append(summed_metros).reset_index(drop=True)
#summed_metros = summed_metros.drop("metro")
# sum admin2 divisions to provinces and backfill in the USA and as needed
summed_provinces = all_data.loc[all_data["admin_level"] == 2].groupby(["iso3", "admin1"]).apply(sumdata).drop(columns=["geometry"])
summed_provinces["admin2"] = None
summed_provinces["admin_level"] = 1
summed_provinces["name"] = summed_provinces["admin1"]
summed_provinces["ID"] = "USA_US-" + summed_provinces["admin1"].apply(tostate)
summed_provinces = pd.merge(summed_provinces, states_geo[["ID", "geometry"]], on="ID", how="left")
# TODO: handle non-US-state (ie CHN, CAN, AUS province) case ^ (lower priority; need geo data)
summed_provinces["location_id"] = summed_provinces.apply(set_loc_id, axis=1)
summed_provinces = summed_provinces.loc[(~summed_provinces["location_id"].isin(all_data["location_id"])) | (summed_provinces["iso3"] == "USA")]
all_data = all_data.loc[~all_data["location_id"].isin(summed_provinces["location_id"])]
all_data = all_data.append(summed_provinces).reset_index(drop=True)

# aggregate provinces to countries and backfill as needed
summed_countries = all_data.loc[all_data["admin_level"] == 1].groupby(["iso3"]).apply(sumdata).drop(columns=["geometry"])
summed_countries["admin1"] = None
summed_countries["admin2"] = None
summed_countries["admin_level"] = 0
summed_countries["ID"] = summed_countries["iso3"]
summed_countries["name"] = summed_countries["iso3"].apply(lambda x: iso3_name_table[x] if x in iso3_name_table else "")
summed_countries = pd.merge(summed_countries, countries_geo[["ID", "geometry"]], on="ID", how="left")
summed_countries["location_id"] = summed_countries.apply(set_loc_id, axis=1)
summed_countries = summed_countries.loc[(~summed_countries["location_id"].isin(all_data["location_id"])) | (summed_countries["iso3"] == "USA")]
all_data = all_data.loc[~all_data["location_id"].isin(summed_countries["location_id"])]
all_data = all_data.append(summed_countries).reset_index(drop=True)
all_data["country_name"] = all_data["iso3"].apply(lambda x: iso3_name_table[x] if x in iso3_name_table else "")

# aggregate countries to regions
all_data["wb_region"] = all_data["iso3"].apply(lambda x: wb_region_table[x] if x in wb_region_table else "")
summed_regions = all_data.loc[all_data["admin_level"] == 0].groupby(["wb_region"]).apply(sumdata).drop(columns=["geometry"])
summed_regions["iso3"] = None
summed_regions["admin1"] = None
summed_regions["admin2"] = None
summed_regions["admin_level"] = -1
summed_regions["ID"] = summed_regions["wb_region"]
summed_regions["name"] = summed_regions["wb_region"]
summed_regions["location_id"] = summed_regions["wb_region"]
all_data = all_data.append(summed_regions).reset_index(drop=True)

# fill missing values
all_data["iso3"] = all_data["iso3"].fillna("None")
all_data["admin1"] = all_data["admin1"].fillna("None")
all_data["admin2"] = all_data["admin2"].fillna("None")
missing_id = pd.isnull(all_data["ID"]) | (len(all_data["ID"]) < 1)
all_data.loc[(all_data["admin_level"] == -1) & missing_id, "ID"] = all_data["wb_region"]
all_data.loc[(all_data["admin_level"] == 0) & missing_id, "ID"] = all_data["iso3"]
all_data.loc[(all_data["admin_level"] == 1) & missing_id, "ID"] = all_data["iso3"] + "_" + all_data["admin1"]
all_data.loc[(all_data["admin_level"] == 2) & missing_id, "ID"] = all_data["iso3"] + "_" + all_data["admin1"] + "_" + all_data["admin2"]
all_data.loc[:, "ID"] = all_data["ID"].apply(lambda x: x.replace(" ", ""))
all_data.loc[pd.isnull(all_data["name"]), "name"] = all_data["ID"]
bpop = loc_lookup[["iso3", "Province_State", "Admin2", "Population"]].rename(columns={"Province_State": "admin1", "Admin2": "admin2", "Population": "bpop"}).fillna("None")
bpop["location_id"] = bpop.apply(set_loc_id, axis=1)
bpop = pd.to_numeric(all_data.merge(bpop, on="location_id", how="left")["bpop"], errors='coerce')
bpop = pd.DataFrame([pd.to_numeric(all_data["population"], errors="coerce"), bpop]).transpose()
all_data.loc[:, "population"] = bpop.max(axis=1, skipna=True)
all_data.loc[:, "location_id"] = all_data["location_id"].apply(lambda x: x.replace(" ", ""))

# query for num_subnational
all_data["num_subnational"] = 0
def get_sub(level):
    def f(x):
        x.loc[pd.isna(x[level]), "num_subnational"] = len(x) - 1
        return x
    return f
all_data = all_data.groupby(["iso3"]).apply(get_sub("admin1"))
all_data = all_data.groupby(["iso3", "admin1"]).apply(get_sub("admin2"))

# fill blank with None and export to next step
print("exporting to parallel", file=sys.stderr)
all_data = all_data.fillna(np.nan).replace([np.nan], [None])
all_data = all_data.rename(columns={"ID":"location_id", "location_id":"ID"})
# all_data.to_json("/dev/stdout", orient="records", lines=True)


all_data['dates_list'] = all_data['dates_list'].map(lambda dateslist: [pd.to_datetime(d).strftime('%Y-%m-%d') for d in dateslist])

all_data.to_json("o6.jsonl", orient="records", lines=True)
