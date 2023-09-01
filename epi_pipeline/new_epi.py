from dask.distributed import Client

client = Client('scheduler:8786')
client

import dask
import dask.dataframe as dd
import pandas as pd
import numpy as np
import json
import sys
from datetime import date, timedelta

loc_lookup = pd.read_csv(workdir+"/UID_ISO_FIPS_LookUp_Table.csv")
ll = loc_lookup.set_index('iso2').sort_values(by='Population', ascending=False)
ll = ll[~ll.index.duplicated(keep='first')].reset_index()
countries_geo = pd.read_csv(rundir+"/countries_geo.tsv", sep="\t", header=None, names=["ID", "name", "geometry"])
countries_geo.loc[:, "geometry"] = countries_geo["geometry"].apply(json.loads)
cg = countries_geo.set_index('ID')
cg = cg[~cg.index.duplicated(keep='first')].reset_index()
# ll = ll[~ll.index.duplicated(keep='first')].reset_index()
wb_region_table = {"AFG":"South Asia","ALB":"Europe & Central Asia","DZA":"Middle East & North Africa","ASM":"East Asia & Pacific","AND":"Europe & Central Asia","AGO":"Sub-Saharan Africa","ATG":"Latin America & Caribbean","ARG":"Latin America & Caribbean","ARM":"Europe & Central Asia","ABW":"Latin America & Caribbean","AUS":"East Asia & Pacific","AUT":"Europe & Central Asia","AZE":"Europe & Central Asia","BHS":"Latin America & Caribbean","BHR":"Middle East & North Africa","BGD":"South Asia","BRB":"Latin America & Caribbean","BLR":"Europe & Central Asia","BEL":"Europe & Central Asia","BLZ":"Latin America & Caribbean","BEN":"Sub-Saharan Africa","BMU":"North America","BTN":"South Asia","BOL":"Latin America & Caribbean","BIH":"Europe & Central Asia","BWA":"Sub-Saharan Africa","BRA":"Latin America & Caribbean","VGB":"Latin America & Caribbean","BRN":"East Asia & Pacific","BGR":"Europe & Central Asia","BFA":"Sub-Saharan Africa","BDI":"Sub-Saharan Africa","CPV":"Sub-Saharan Africa","KHM":"East Asia & Pacific","CMR":"Sub-Saharan Africa","CAN":"North America","CYM":"Latin America & Caribbean","CAF":"Sub-Saharan Africa","TCD":"Sub-Saharan Africa","CHI":"Europe & Central Asia","CHL":"Latin America & Caribbean","CHN":"East Asia & Pacific","COL":"Latin America & Caribbean","COM":"Sub-Saharan Africa","COD":"Sub-Saharan Africa","COG":"Sub-Saharan Africa","CRI":"Latin America & Caribbean","CIV":"Sub-Saharan Africa","HRV":"Europe & Central Asia","CUB":"Latin America & Caribbean","CUW":"Latin America & Caribbean","CYP":"Europe & Central Asia","CZE":"Europe & Central Asia","DNK":"Europe & Central Asia","DJI":"Middle East & North Africa","DMA":"Latin America & Caribbean","DOM":"Latin America & Caribbean","ECU":"Latin America & Caribbean","EGY":"Middle East & North Africa","SLV":"Latin America & Caribbean","GNQ":"Sub-Saharan Africa","ERI":"Sub-Saharan Africa","EST":"Europe & Central Asia","SWZ":"Sub-Saharan Africa","ETH":"Sub-Saharan Africa","FRO":"Europe & Central Asia","FJI":"East Asia & Pacific","FIN":"Europe & Central Asia","FRA":"Europe & Central Asia","PYF":"East Asia & Pacific","GAB":"Sub-Saharan Africa","GMB":"Sub-Saharan Africa","GEO":"Europe & Central Asia","DEU":"Europe & Central Asia","GHA":"Sub-Saharan Africa","GIB":"Europe & Central Asia","GRC":"Europe & Central Asia","GRL":"Europe & Central Asia","GRD":"Latin America & Caribbean","GUM":"East Asia & Pacific","GTM":"Latin America & Caribbean","GIN":"Sub-Saharan Africa","GNB":"Sub-Saharan Africa","GUY":"Latin America & Caribbean","HTI":"Latin America & Caribbean","HND":"Latin America & Caribbean","HKG":"East Asia & Pacific","HUN":"Europe & Central Asia","ISL":"Europe & Central Asia","IND":"South Asia","IDN":"East Asia & Pacific","IRN":"Middle East & North Africa","IRQ":"Middle East & North Africa","IRL":"Europe & Central Asia","IMN":"Europe & Central Asia","ISR":"Middle East & North Africa","ITA":"Europe & Central Asia","JAM":"Latin America & Caribbean","JPN":"East Asia & Pacific","JOR":"Middle East & North Africa","KAZ":"Europe & Central Asia","KEN":"Sub-Saharan Africa","KIR":"East Asia & Pacific","PRK":"East Asia & Pacific","KOR":"East Asia & Pacific","XKX":"Europe & Central Asia","KWT":"Middle East & North Africa","KGZ":"Europe & Central Asia","LAO":"East Asia & Pacific","LVA":"Europe & Central Asia","LBN":"Middle East & North Africa","LSO":"Sub-Saharan Africa","LBR":"Sub-Saharan Africa","LBY":"Middle East & North Africa","LIE":"Europe & Central Asia","LTU":"Europe & Central Asia","LUX":"Europe & Central Asia","MAC":"East Asia & Pacific","MDG":"Sub-Saharan Africa","MWI":"Sub-Saharan Africa","MYS":"East Asia & Pacific","MDV":"South Asia","MLI":"Sub-Saharan Africa","MLT":"Middle East & North Africa","MHL":"East Asia & Pacific","MRT":"Sub-Saharan Africa","MUS":"Sub-Saharan Africa","MEX":"Latin America & Caribbean","FSM":"East Asia & Pacific","MDA":"Europe & Central Asia","MCO":"Europe & Central Asia","MNG":"East Asia & Pacific","MNE":"Europe & Central Asia","MAR":"Middle East & North Africa","MOZ":"Sub-Saharan Africa","MMR":"East Asia & Pacific","NAM":"Sub-Saharan Africa","NRU":"East Asia & Pacific","NPL":"South Asia","NLD":"Europe & Central Asia","NCL":"East Asia & Pacific","NZL":"East Asia & Pacific","NIC":"Latin America & Caribbean","NER":"Sub-Saharan Africa","NGA":"Sub-Saharan Africa","MKD":"Europe & Central Asia","MNP":"East Asia & Pacific","NOR":"Europe & Central Asia","OMN":"Middle East & North Africa","PAK":"South Asia","PLW":"East Asia & Pacific","PAN":"Latin America & Caribbean","PNG":"East Asia & Pacific","PRY":"Latin America & Caribbean","PER":"Latin America & Caribbean","PHL":"East Asia & Pacific","POL":"Europe & Central Asia","PRT":"Europe & Central Asia","PRI":"Latin America & Caribbean","QAT":"Middle East & North Africa","ROU":"Europe & Central Asia","RUS":"Europe & Central Asia","RWA":"Sub-Saharan Africa","WSM":"East Asia & Pacific","SMR":"Europe & Central Asia","STP":"Sub-Saharan Africa","SAU":"Middle East & North Africa","SEN":"Sub-Saharan Africa","SRB":"Europe & Central Asia","SYC":"Sub-Saharan Africa","SLE":"Sub-Saharan Africa","SGP":"East Asia & Pacific","SXM":"Latin America & Caribbean","SVK":"Europe & Central Asia","SVN":"Europe & Central Asia","SLB":"East Asia & Pacific","SOM":"Sub-Saharan Africa","ZAF":"Sub-Saharan Africa","SSD":"Sub-Saharan Africa","ESP":"Europe & Central Asia","LKA":"South Asia","KNA":"Latin America & Caribbean","LCA":"Latin America & Caribbean","MAF":"Latin America & Caribbean","VCT":"Latin America & Caribbean","SDN":"Sub-Saharan Africa","SUR":"Latin America & Caribbean","SWE":"Europe & Central Asia","CHE":"Europe & Central Asia","SYR":"Middle East & North Africa","TWN":"East Asia & Pacific","TJK":"Europe & Central Asia","TZA":"Sub-Saharan Africa","THA":"East Asia & Pacific","TLS":"East Asia & Pacific","TGO":"Sub-Saharan Africa","TON":"East Asia & Pacific","TTO":"Latin America & Caribbean","TUN":"Middle East & North Africa","TUR":"Europe & Central Asia","TKM":"Europe & Central Asia","TCA":"Latin America & Caribbean","TUV":"East Asia & Pacific","UGA":"Sub-Saharan Africa","UKR":"Europe & Central Asia","ARE":"Middle East & North Africa","GBR":"Europe & Central Asia","USA":"North America","URY":"Latin America & Caribbean","UZB":"Europe & Central Asia","VUT":"East Asia & Pacific","VEN":"Latin America & Caribbean","VNM":"East Asia & Pacific","VIR":"Latin America & Caribbean","PSE":"Middle East & North Africa","YEM":"Middle East & North Africa","ZMB":"Sub-Saharan Africa","ZWE":"Sub-Saharan Africa"}
iso3_name_table = dict(zip(countries_geo["ID"], countries_geo["name"]))


who_cases = dd.read_csv('https://covid19.who.int/WHO-COVID-19-global-data.csv')
cdc_cases = dd.read_csv('https://data.cdc.gov/api/views/vbim-akqf/rows.csv?accessType=DOWNLOAD', blocksize=None, assume_missing=True)
cdc_case_counts = cdc_cases[['cdc_case_earliest_dt ', 'current_status']].groupby('cdc_case_earliest_dt ').count()

def group_dates_who(x):
    b = x.iloc[0].copy().drop(['New_cases', 'Cumulative_cases', 'New_deaths', 'Cumulative_deaths', 'Date_reported'])
    b['dates_list'] = x.loc[:, 'Date_reported'].values
    b['cases'] = x.loc[:, 'New_cases'].fillna(0).values
    b['losses'] = x.loc[:, 'New_deaths'].fillna(0).values
    return b

c=who_cases.groupby('Country_code').apply(group_dates_who).compute()

cdc_case_counts = cdc_cases[['cdc_case_earliest_dt ', 'current_status']].groupby('cdc_case_earliest_dt ').count()

cdc_case_counts = cdc_case_counts.sort_index()

usa_dates = cdc_case_counts.index.to_numpy()
usa_cases = cdc_case_counts['current_status'].to_numpy()

c2 = pd.concat([c[c['Country_code'] != 'US'], pd.DataFrame([['US', 'United States of America', 'AMRO', usa_dates, usa_cases, np.zeros_like(usa_cases)]], columns=c.columns)])


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

global_epi = c2

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


all_data = global_data

# # gather global and US data
# print("gathering and patching", file=sys.stderr)
# all_data = global_data.append(county_data)
all_data = all_data.loc[~pd.isnull(all_data["iso3"])]
all_data = all_data.fillna(np.nan).replace([np.nan], [None])
set_loc_id = lambda x: "{}_{}_{}".format(x["iso3"], x["admin1"], x["admin2"])
all_data["location_id"] = all_data.apply(set_loc_id, axis=1)
# all_data.loc[all_data["iso3"] == "USA", "country_name"] = "United States of America"
all_data["admin_level"] = 2
all_data["FIPS"] = pd.to_numeric(all_data["FIPS"])
all_data.loc[pd.isna(all_data["admin2"]), "admin_level"] = 1
all_data.loc[pd.isna(all_data["admin1"]), "admin_level"] = 0

all_data["country_name"] = all_data["iso3"].apply(lambda x: iso3_name_table[x] if x in iso3_name_table else "")

all_data["wb_region"] = all_data["iso3"].apply(lambda x: wb_region_table[x] if x in wb_region_table else "")

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


all_data["num_subnational"] = 0
def get_sub(level):
    def f(x):
        x.loc[pd.isna(x[level]), "num_subnational"] = len(x) - 1
        return x
    return f
all_data = all_data.groupby(["iso3"]).apply(get_sub("admin1"))
all_data = all_data.groupby(["iso3", "admin1"]).apply(get_sub("admin2"))

# fill blank with None and export to next step
all_data = all_data.fillna(np.nan).replace([np.nan], [None])
all_data = all_data.rename(columns={"ID":"location_id", "location_id":"ID"})

all_data['dates_list'] = all_data['dates_list'].map(lambda dateslist: [pd.to_datetime(d).strftime('%Y-%m-%d') for d in dateslist])

all_data.to_json("epi.out.jsonl", orient="records", lines=True)