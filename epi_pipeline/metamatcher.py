#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import json
import sys

workdir = sys.argv[1]
rundir = sys.argv[2]

# prep geo data
print("prepping geo data", file=sys.stderr)
loc_lookup = pd.read_csv(workdir+"/UID_ISO_FIPS_LookUp_Table.csv")
countries_geo = pd.read_csv(rundir+"/countries_geo.tsv", sep="\t", header=None, names=["ID", "name", "geometry"])
countries_geo.loc[:, "geometry"] = countries_geo["geometry"].apply(json.loads)
wb_region_table = {"AFG":"South Asia","ALB":"Europe & Central Asia","DZA":"Middle East & North Africa","ASM":"East Asia & Pacific","AND":"Europe & Central Asia","AGO":"Sub-Saharan Africa","ATG":"Latin America & Caribbean","ARG":"Latin America & Caribbean","ARM":"Europe & Central Asia","ABW":"Latin America & Caribbean","AUS":"East Asia & Pacific","AUT":"Europe & Central Asia","AZE":"Europe & Central Asia","BHS":"Latin America & Caribbean","BHR":"Middle East & North Africa","BGD":"South Asia","BRB":"Latin America & Caribbean","BLR":"Europe & Central Asia","BEL":"Europe & Central Asia","BLZ":"Latin America & Caribbean","BEN":"Sub-Saharan Africa","BMU":"North America","BTN":"South Asia","BOL":"Latin America & Caribbean","BIH":"Europe & Central Asia","BWA":"Sub-Saharan Africa","BRA":"Latin America & Caribbean","VGB":"Latin America & Caribbean","BRN":"East Asia & Pacific","BGR":"Europe & Central Asia","BFA":"Sub-Saharan Africa","BDI":"Sub-Saharan Africa","CPV":"Sub-Saharan Africa","KHM":"East Asia & Pacific","CMR":"Sub-Saharan Africa","CAN":"North America","CYM":"Latin America & Caribbean","CAF":"Sub-Saharan Africa","TCD":"Sub-Saharan Africa","CHI":"Europe & Central Asia","CHL":"Latin America & Caribbean","CHN":"East Asia & Pacific","COL":"Latin America & Caribbean","COM":"Sub-Saharan Africa","COD":"Sub-Saharan Africa","COG":"Sub-Saharan Africa","CRI":"Latin America & Caribbean","CIV":"Sub-Saharan Africa","HRV":"Europe & Central Asia","CUB":"Latin America & Caribbean","CUW":"Latin America & Caribbean","CYP":"Europe & Central Asia","CZE":"Europe & Central Asia","DNK":"Europe & Central Asia","DJI":"Middle East & North Africa","DMA":"Latin America & Caribbean","DOM":"Latin America & Caribbean","ECU":"Latin America & Caribbean","EGY":"Middle East & North Africa","SLV":"Latin America & Caribbean","GNQ":"Sub-Saharan Africa","ERI":"Sub-Saharan Africa","EST":"Europe & Central Asia","SWZ":"Sub-Saharan Africa","ETH":"Sub-Saharan Africa","FRO":"Europe & Central Asia","FJI":"East Asia & Pacific","FIN":"Europe & Central Asia","FRA":"Europe & Central Asia","PYF":"East Asia & Pacific","GAB":"Sub-Saharan Africa","GMB":"Sub-Saharan Africa","GEO":"Europe & Central Asia","DEU":"Europe & Central Asia","GHA":"Sub-Saharan Africa","GIB":"Europe & Central Asia","GRC":"Europe & Central Asia","GRL":"Europe & Central Asia","GRD":"Latin America & Caribbean","GUM":"East Asia & Pacific","GTM":"Latin America & Caribbean","GIN":"Sub-Saharan Africa","GNB":"Sub-Saharan Africa","GUY":"Latin America & Caribbean","HTI":"Latin America & Caribbean","HND":"Latin America & Caribbean","HKG":"East Asia & Pacific","HUN":"Europe & Central Asia","ISL":"Europe & Central Asia","IND":"South Asia","IDN":"East Asia & Pacific","IRN":"Middle East & North Africa","IRQ":"Middle East & North Africa","IRL":"Europe & Central Asia","IMN":"Europe & Central Asia","ISR":"Middle East & North Africa","ITA":"Europe & Central Asia","JAM":"Latin America & Caribbean","JPN":"East Asia & Pacific","JOR":"Middle East & North Africa","KAZ":"Europe & Central Asia","KEN":"Sub-Saharan Africa","KIR":"East Asia & Pacific","PRK":"East Asia & Pacific","KOR":"East Asia & Pacific","XKX":"Europe & Central Asia","KWT":"Middle East & North Africa","KGZ":"Europe & Central Asia","LAO":"East Asia & Pacific","LVA":"Europe & Central Asia","LBN":"Middle East & North Africa","LSO":"Sub-Saharan Africa","LBR":"Sub-Saharan Africa","LBY":"Middle East & North Africa","LIE":"Europe & Central Asia","LTU":"Europe & Central Asia","LUX":"Europe & Central Asia","MAC":"East Asia & Pacific","MDG":"Sub-Saharan Africa","MWI":"Sub-Saharan Africa","MYS":"East Asia & Pacific","MDV":"South Asia","MLI":"Sub-Saharan Africa","MLT":"Middle East & North Africa","MHL":"East Asia & Pacific","MRT":"Sub-Saharan Africa","MUS":"Sub-Saharan Africa","MEX":"Latin America & Caribbean","FSM":"East Asia & Pacific","MDA":"Europe & Central Asia","MCO":"Europe & Central Asia","MNG":"East Asia & Pacific","MNE":"Europe & Central Asia","MAR":"Middle East & North Africa","MOZ":"Sub-Saharan Africa","MMR":"East Asia & Pacific","NAM":"Sub-Saharan Africa","NRU":"East Asia & Pacific","NPL":"South Asia","NLD":"Europe & Central Asia","NCL":"East Asia & Pacific","NZL":"East Asia & Pacific","NIC":"Latin America & Caribbean","NER":"Sub-Saharan Africa","NGA":"Sub-Saharan Africa","MKD":"Europe & Central Asia","MNP":"East Asia & Pacific","NOR":"Europe & Central Asia","OMN":"Middle East & North Africa","PAK":"South Asia","PLW":"East Asia & Pacific","PAN":"Latin America & Caribbean","PNG":"East Asia & Pacific","PRY":"Latin America & Caribbean","PER":"Latin America & Caribbean","PHL":"East Asia & Pacific","POL":"Europe & Central Asia","PRT":"Europe & Central Asia","PRI":"Latin America & Caribbean","QAT":"Middle East & North Africa","ROU":"Europe & Central Asia","RUS":"Europe & Central Asia","RWA":"Sub-Saharan Africa","WSM":"East Asia & Pacific","SMR":"Europe & Central Asia","STP":"Sub-Saharan Africa","SAU":"Middle East & North Africa","SEN":"Sub-Saharan Africa","SRB":"Europe & Central Asia","SYC":"Sub-Saharan Africa","SLE":"Sub-Saharan Africa","SGP":"East Asia & Pacific","SXM":"Latin America & Caribbean","SVK":"Europe & Central Asia","SVN":"Europe & Central Asia","SLB":"East Asia & Pacific","SOM":"Sub-Saharan Africa","ZAF":"Sub-Saharan Africa","SSD":"Sub-Saharan Africa","ESP":"Europe & Central Asia","LKA":"South Asia","KNA":"Latin America & Caribbean","LCA":"Latin America & Caribbean","MAF":"Latin America & Caribbean","VCT":"Latin America & Caribbean","SDN":"Sub-Saharan Africa","SUR":"Latin America & Caribbean","SWE":"Europe & Central Asia","CHE":"Europe & Central Asia","SYR":"Middle East & North Africa","TWN":"East Asia & Pacific","TJK":"Europe & Central Asia","TZA":"Sub-Saharan Africa","THA":"East Asia & Pacific","TLS":"East Asia & Pacific","TGO":"Sub-Saharan Africa","TON":"East Asia & Pacific","TTO":"Latin America & Caribbean","TUN":"Middle East & North Africa","TUR":"Europe & Central Asia","TKM":"Europe & Central Asia","TCA":"Latin America & Caribbean","TUV":"East Asia & Pacific","UGA":"Sub-Saharan Africa","UKR":"Europe & Central Asia","ARE":"Middle East & North Africa","GBR":"Europe & Central Asia","USA":"North America","URY":"Latin America & Caribbean","UZB":"Europe & Central Asia","VUT":"East Asia & Pacific","VEN":"Latin America & Caribbean","VNM":"East Asia & Pacific","VIR":"Latin America & Caribbean","PSE":"Middle East & North Africa","YEM":"Middle East & North Africa","ZMB":"Sub-Saharan Africa","ZWE":"Sub-Saharan Africa"}
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
global_epi = pd.read_csv(workdir+"/time_series_covid19_confirmed_global.csv")
dates_list_global = pd.to_datetime(list(global_epi.loc[:, "1/22/20":]), format="%m/%d/%y").strftime("%Y-%m-%d")
global_cases = list(global_epi.loc[:, "1/22/20":].to_numpy())
global_epi = global_epi.loc[:, :"Long"]
global_epi["cases"] = global_cases
global_epi["losses"] = list(pd.read_csv(workdir+"/time_series_covid19_deaths_global.csv").loc[:, "1/22/20":].to_numpy())

# global geo matcher
global_epi["UID"] = global_epi.index.values
global_data = pd.merge(global_epi, loc_lookup[["UID", "Population", "FIPS"]], on="UID", how="left")
global_data = global_data.rename(columns={"Country/Region":"country_name","Province/State":"admin1","Lat":"lat","Long":"long","Population":"population"})
global_data.loc[:, "country_name"] = global_data["country_name"].str.replace("*","", regex=False)
global_data["name"] = global_data["country_name"]
global_data.loc[~pd.isnull(global_data["admin1"]), "name"] = global_data["admin1"]
global_data = pd.merge(global_data, countries_geo[["name", "ID"]].rename(columns={"name":"country_name"}), on="country_name", how="left")
global_data = pd.merge(global_data, countries_geo[["name", "geometry"]], on="name", how="left")
global_data["iso3"] = global_data["ID"]
global_data["ID"] = None
global_data["admin2"] = None
global_data["dates_list"] = [dates_list_global]*len(global_data)
global_data = global_data[metakeys]
#missing_filter = global_data.drop(columns=["admin1", "population"]).isna().any(axis=1)
#missing_provinces = global_data[missing_filter]
#global_data = global_data[~missing_filter]

# US data read-in & numpy conversion
print("reading in county-level epi data", file=sys.stderr)
counties_epi = pd.read_csv(workdir+"/time_series_covid19_confirmed_US.csv")
dates_list_counties = pd.to_datetime(list(counties_epi.loc[:, "1/22/20":]), format="%m/%d/%y").strftime("%Y-%m-%d")
counties_cases = list(counties_epi.loc[:, "1/22/20":].to_numpy())
counties_epi = counties_epi.loc[:, :"Long_"]
counties_epi["cases"] = counties_cases
counties_epi["losses"] = list(pd.read_csv(workdir+"/time_series_covid19_deaths_US.csv").loc[:, "1/22/20":].to_numpy())

# US geo matcher
counties_epi["ID"] = "USA_US-" + counties_epi["Province_State"].apply(tostate) + "_" + \
    counties_epi[["FIPS", "Admin2"]].apply(lambda x: "{:05d}".format(int(x["FIPS"])) if pd.notna(x["FIPS"]) else x["Admin2"], axis=1)
county_data = pd.merge(counties_epi, counties_geo, on="ID", how="left").drop(columns=["FIPS"])
county_data = pd.merge(county_data, loc_lookup[["UID", "Population", "FIPS"]], on="UID", how="left")
#print(counties_epi["FIPS"], file=sys.stderr)
county_data = county_data.rename(columns={"Admin2":"admin2","Province_State":"admin1","Lat":"lat","Long_":"long","Population":"population"})
county_data["dates_list"] = [dates_list_counties]*len(county_data)
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
all_data.to_json("/dev/stdout", orient="records", lines=True)

