"""
Use these utils to upload sequences to NCBI
"""

import os
import shutil
from typing import Dict, List
import json

import pandas as pd

from gsheet_interact import gisaid_interactor

def merge_metadata(repo_metadata: pd.DataFrame, online_metadata: pd.DataFrame) -> pd.DataFrame:
    """
    Takes the online metadata and remerges it with that published in the HCoV-19 repo
    to create a combined file with all needed columns
    """
    merged_data = repo_metadata.merge(online_metadata, how="left", left_on="fasta_hdr", right_on="Virus name")
    return merged_data


def convert_metadata(
    metadata: pd.DataFrame,  
    column_config_file: str,
    constant_config_file: str,
    wastewater: bool = False,
) -> pd.DataFrame:
    """
    Takes the metadata given and converts it into the ncbi
    format
    """
    # load the config file
    if wastewater:
        column_mapping = json.loads(column_config_file)["wastewater"]
        constant_mapping = json.loads(constant_config_file)["wastewater"]
    else:
        column_mapping = json.loads(column_config_file)["normal"]
        constant_mapping = json.loads(constant_config_file)["normal"]
    
    # convert columns
    pass


def batch_by_author(metadata: pd.DataFrame) -> List[pd.DataFrame]:
    """
    Takes the given metadata and returns it as a list of dataframes, each with one set of
    authors
    """
    pass


def separate_wastewater(
    wastewater_sample_list_path: str,
    all_sequences_location: str,
    separated_location: str,
) -> None:
    """
    Pull out the relevant wastewater sequences for separate processing
    """
    # load the list of samples from a data file
    with open(wastewater_sample_list_path, "r") as infile:
        lines = infile.readlines()
        wastewater_samples_list = [line.rstrip() for line in lines]

    # in the location specified, search for these files and extract them into a separate location
    if not os.path.exists(separated_location):
        os.mkdir(separated_location)

    for sequence_file in wastewater_samples_list:
        shutil.copy(
            os.path.join(all_sequences_location, sequence_file), separated_location
        )

    return

def separate_normal(
    normal_sample_list_path: str,
    all_sequences_location: str,
    separated_location: str
) -> None:
    """
    Pull out the relevant normal sequences for separate processing
    """
    # load the list of samples from a data file
    with open(normal_sample_list_path, "r") as infile:
        lines = infile.readlines()
        normal_samples_list = [line.rstrip() for line in lines]

    # in the location specified, search for these files and extract them into a separate location
    if not os.path.exists(separated_location):
        os.mkdir(separated_location)

    for sequence_file in normal_samples_list:
        shutil.copy(
            os.path.join(all_sequences_location, sequence_file), separated_location
        )

    return 


def extract_sequences():
    """
    Not quite sure how this would work exactly, but it would need to compare what
    we have already uploaded vs what we have left to upload and then return a
    relevant sheet of the associated metadata with those sequences and then sequences
    / bam files themselves.
    """
    pass


def extract_bam_files():
    """
    Same as above - maybe this is separate from bam files and sequences?
    """
    pass


def generate_bucket_list():
    """
    We're going to need someway to query our google bucket and see
    what new sequences still need to go? Not sure what the ncbi portal looks like from this aspect
    so tough to know what we need to automate for.
    """
    pass


def map_biosample_ids(
    metadata: pd.DataFrame, biosample_data: pd.DataFrame, sequences: str, bam_files: str
) -> None:
    """
    Need to take the biosample data, our metadata, and use all of that to map the sequences and the bam file for upload
    This should ideally happen after batching by author, but doesn't matter a ton?
    This would create a folder structure in a centralized place that we could upload from
    Potentially could interact with gcp to restructe the folders in our storage buckets to make this easy
    """
    pass
