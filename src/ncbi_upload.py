"""
Use these utils to upload sequences to NCBI
"""

import pandas as pd
from gsheet_interact import gisaid_interactor

from typing import Dict, List

def convert_metadata(metadata: pd.DataFrame, output_mapping: Dict[str, str]) -> pd.DataFrame:
    """
    Takes the metadata given and converts it into the ncbi
    format
    """
    pass

def batch_by_author(metadata: pd.DataFrame) -> List[pd.DataFrame]:
    """
    Takes the given metadata and returns it as a list of dataframes, each with one set of
    authors
    """
    pass

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


def map_biosample_ids(metadata: pd.DataFrame, biosample_data: pd.DataFrame, sequences: str, bam_files: str) -> None:
    """
    Need to take the biosample data, our metadata, and use all of that to map the sequences and the bam file for upload
    This should ideally happen after batching by author, but doesn't matter a ton?
    This would create a folder structure in a centralized place that we could upload from
    Potentially could interact with gcp to restructe the folders in our storage buckets to make this easy
    """
    pass