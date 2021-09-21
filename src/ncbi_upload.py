"""
Use these utils to upload sequences to NCBI
"""

import os
import shutil
from typing import Collection, Dict, List
import time
import json

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
    converted_metadata = convert_metadata.rename(columns={"authors_new": "collected_by_1"})

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

    # convert vaccination text to timestamp
    # converted_metadata_2["last_vaccinated"] = _convert_str_date_to_timestamp(converted_metadata_2["collection_date"], converted_metadata_2["last_vaccinated_raw"])
    
    # drop defunct columns
    converted_metadata_3 = converted_metadata_2.drop(columns=["collected_by_1", "collected_by_2"])

    return converted_metadata_3


def _convert_str_date_to_timestamp(text: pd.Series, collection_date: pd.Series) -> pd.Series:
    """
    Take a text field which says how many weeks or months ago someone was vaccinated and transform that into a timestamp
    and then render that timestamp as text
    """
    #TODO: Needs to handle errors in the collection_date
    #TODO: Needs to handle errors in vaccination date
    collection_timestamps = time.mktime(time.strptime(collection_date, '%d-%m-%Y'))
    pass


def batch_by_author(metadata: pd.DataFrame) -> List[pd.DataFrame]:
    """
    Takes the given metadata and returns it as a list of dataframes, each with one set of
    authors
    """
    pass