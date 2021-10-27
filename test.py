"""
Unit tests for location information.
"""

import os
import sys
import unittest
import pandas as pd

sys.path.insert(1, './python/')
from  normalize_and_merge import nextstrain_replacement, parse_jsonl, off_by_one_location

class TestMetadataMethods(unittest.TestCase):
    def test_locstring(self):
        metadata_filename = os.path.join("unit_test_files", "potential_locstrings.csv")
        meta = pd.read_csv(metadata_filename)
        meta = nextstrain_replacement(os.path.join("data", "gisaid_geoLocationRules.tsv"), meta)
        std_locs = parse_jsonl(os.path.join("data", "gadm_transformed.jsonl"))
        meta = off_by_one_location(meta, std_locs)
        meta.to_csv(os.path.join("unit_test_files", "unit_test_locations.csv"))

if __name__ == '__main__':
    object_1 = TestMetadataMethods()
    object_1.test_locstring()

