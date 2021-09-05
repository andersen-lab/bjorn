"""
Interact with gsheets for automation sake
"""

import pandas as pd
import gspread

def gisaid_interactor():
    """
    Interact with a metadata file and get the appropriate results
    Split this out into separate file when we're done
    """
    gisaid_key = "1alCHJC5WZEzS8wKCzM_-CgvOO67K3piJ5wKhOQ9KXOM"
    gisaid_wksht_num = 1
    metadata = get_gsheet(gisaid_key, gisaid_wksht_num)
    return metadata

def get_gsheet(file_key: str, worksheet_num: int):
    """
    get from gsheet
    """
    gc = gspread.service_account()
    worksheet = gc.open_by_key(key).get_worksheet(worksheet_num)
    df = pd.DataFrame(worksheet.get_all_records())
    return df

def push_gsheet():
    """
    push data to gsheet
    """
    pass