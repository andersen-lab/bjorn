"""
Interact with gsheets for automation sake
"""

import pandas as pd
import gspread
from configparser import ConfigParser
from typing import Dict
import sys

def _get_config(config_file_path: str = "../bjorn.ini") -> Dict[str, str]:
    """
    gets the config file and loads the paths of relevant information
    """
    config = ConfigParser()
    config.read(config_file_path)
    return {"gisaid_key": config["gsheets"]["gisaid_key"],
            "gisaid_wksht_num": int(config["gsheets"]["gisaid_wksht_num"]),
            "gsheet_key_path": config["gsheets"]["gsheet_key_path"]}

def gisaid_interactor(config_file_path: str = "../bjorn.ini") -> pd.DataFrame:
    """
    Interact with a metadata file and get the appropriate results
    Split this out into separate file when we're done
    """
    config = _get_config(config_file_path)
    metadata = _get_gsheet(config['gisaid_key'], config['gisaid_wksht_num'], config['gsheet_key_path'])
    return metadata

def _get_gsheet(file_key: str, worksheet_num: int, service_account_json: str) -> pd.DataFrame:
    """
    get from gsheet
    """
    gc = gspread.service_account(filename = service_account_json)
    worksheet = gc.open_by_key(file_key).get_worksheet(worksheet_num)
    return pd.DataFrame(worksheet.get_all_records())

def _push_gsheet():
    """
    push data to gsheet
    """
    pass

if __name__=="__main__":
    data = gisaid_interactor(sys.argv[1])
    data.to_csv(sys.argv[2])