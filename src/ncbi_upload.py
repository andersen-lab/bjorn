"""
Use these utils to upload sequences to NCBI
"""

import datetime
import json
import os
import shutil
from typing import Collection, Dict, List, Tuple

import pandas as pd

#from gsheet_interact import gisaid_interactor

def merge_metadata(repo_metadata: pd.DataFrame, online_metadata: pd.DataFrame) -> pd.DataFrame:
    """
    Takes the online metadata and remerges it with that published in the HCoV-19 repo
    to create a combined file with all needed columns
    """
    #TODO: Add code here to get the gsheets from online and do this in google cloud instead
    merged_data = repo_metadata.merge(online_metadata, how="left", left_on="fasta_hdr", right_on="Virus name")
    return merged_data


def convert_metadata(
    metadata: pd.DataFrame,  
    column_config_file: str,
    constant_config_file: str,
    author_conversion_file: str,
    wastewater: bool = False,
) -> pd.DataFrame:
    """
    Takes the metadata given and converts it into the ncbi
    format
    """
    # load the config file
    if wastewater:
        with open(column_config_file, 'r') as infile:
            column_mapping = json.load(infile)['wastewater']
        with open(constant_config_file, 'r') as infile:
            constant_mapping = json.load(infile)["wastewater"]
    else:
        with open(column_config_file, 'r') as infile:
            column_mapping = json.load(infile)["normal"]
        with open(constant_config_file, 'r') as infile:
            constant_mapping = json.load(infile)["normal"]

    # convert columns as per config file
    converted_metadata = pd.DataFrame()
    for key in column_mapping:
        converted_metadata[key] = metadata[column_mapping[key]]
    for key in constant_mapping:
        converted_metadata[key] = constant_mapping[key]
    
    # fix author and originating lab fields
    # map old author names to the new version of those names
    author_conversions = pd.read_csv(author_conversion_file)
    converted_metadata = converted_metadata.merge(author_conversions, how="left", left_on="collected_by_1", right_on="authors_original")
    converted_metadata = converted_metadata.drop(columns=["collected_by_1", "authors_original"])
    converted_metadata = converted_metadata.rename(columns={"authors_new": "collected_by_1"})

    # split the dataframe by if it has both fields, one, or the other
    has_both_fields = converted_metadata[(~converted_metadata["collected_by_1"].isna()) & (~converted_metadata["collected_by_2"].isna())]
    has_both_fields["collected_by"] = has_both_fields["collected_by_1"] + " with the help of " + has_both_fields["collected_by_2"]
    
    # has the first field and not the second
    has_first_field = converted_metadata[(~converted_metadata["collected_by_1"].isna()) & (converted_metadata["collected_by_2"].isna())]
    has_first_field["collected_by"] = has_first_field["collected_by_1"]
    
    # has the second field and not the first
    has_second_field = converted_metadata[(converted_metadata["collected_by_1"].isna()) & (~converted_metadata["collected_by_2"].isna())]
    has_second_field["collected_by"] = has_second_field["collected_by_2"]

    # if neither
    has_neither_field = converted_metadata[(converted_metadata["collected_by_1"].isna()) & (converted_metadata["collected_by_2"].isna())]
    has_neither_field["collected_by"] = "Unknown"

    # combine the 4 dataframes above 
    converted_metadata_2 = pd.concat([has_both_fields, has_first_field, has_second_field, has_neither_field])

    # convert location to have ':' based separation
    converted_metadata_2["geo_loc_name"] = converted_metadata_2["geo_loc_name"].str.replace("/", ":")

    # drop defunct columns
    converted_metadata_3 = converted_metadata_2.drop(columns=["collected_by_1", "collected_by_2"])

    return converted_metadata_3

def recover_vaccine_date(metadata: pd.DataFrame) -> pd.DataFrame:
    """
    Use this separate function to get the actual date out so we can clean up the merged file for errors that would break this
    """
    needs_to_be_vaccine_timestamped = metadata[metadata["last_vaccinated_raw"].str.contains("~", na=False)]
    needs_to_be_vaccine_timestamped["date_of_sars_cov_2_vaccination"], needs_to_be_vaccine_timestamped["vaccine_received"] = _convert_str_date_to_timestamp(needs_to_be_vaccine_timestamped["last_vaccinated_raw"], needs_to_be_vaccine_timestamped["collection_date"])
    needs_to_be_vaccine_timestamped["date_of_sars_cov_2_vaccination"] = needs_to_be_vaccine_timestamped["date_of_sars_cov_2_vaccination"].dt.strftime('%Y-%m-%d')
    does_not_need_to_be_vaccine_timestamped = metadata[~metadata["last_vaccinated_raw"].str.contains("~", na=False)]
    timestamped_metadata = pd.concat([needs_to_be_vaccine_timestamped, does_not_need_to_be_vaccine_timestamped])
    #timestamped_metadata = timestamped_metadata.drop(columns=["last_vaccinated_raw"])
    return timestamped_metadata

def _convert_str_date_to_timestamp(text: pd.Series, collection_date: pd.Series) -> Tuple[pd.Series, pd.Series]:
    """
    Take a text field which says how many weeks or months ago someone was vaccinated and transform that into a timestamp
    and then render that timestamp as text
    """
    # make pairs out of all the values
    collection_timestamps = pd.to_datetime(collection_date).to_list()

    split_series = text.str.split(" ")
    # number of weeks or months
    time_change_numeral = [float(item[0].replace("~", "")) for item in split_series]

    # the weeks or months
    period = [item[1] for item in split_series]
    
    # first or second dose
    dose = [item[3] for item in split_series]

    # vaccine type 
    vaccine = [item[4] if item[4] != "&" else "Johnson & Johnson" for item in split_series]
    
    list_of_vaccine_dates = []
    for entry in zip(collection_timestamps, time_change_numeral, period, dose, vaccine):
        
        # get if week or months
        if "week" in entry[2]:
            multiplier = 7
        else:
            multiplier = 30
        
        # get vaccine date
        if entry[3] == "second":
            vaccine_date = entry[0] - pd.Timedelta(days=(entry[1]*multiplier)+30)
        else:
            vaccine_date = entry[0] - pd.Timedelta(days=entry[1]*multiplier)    
        
        list_of_vaccine_dates.append(vaccine_date)
        
    return list_of_vaccine_dates, vaccine



def batch_by_author(metadata: pd.DataFrame) -> List[pd.DataFrame]:
    """
    Takes the given metadata and returns it as a list of dataframes, each with one set of
    authors
    """
    pass