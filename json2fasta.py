import gc
import re
import json
import fiona
import pandas as pd
import geopandas as gpd
import bjorn_support as bs
import visualize as bv
import data as bd


with open('config.json', 'r') as f:
    config = json.load(f)

in_fp = config['gisaid_feed']
out_fp = config['gisaid_fasta']
metacols = ['covv_virus_name', 'covsurver_prot_mutations', 'covv_location',
             'covv_lineage', 'covv_collection_date', 'covv_accession_id',
             'pangolin_lineages_version', 'covv_clade', 'covv_subm_date']
data = [json.loads(line) for line in open(in_fp, 'r')]
print(f"Total number of sequences: {len(data)}")
print(f"Converting to dict...")
regex = re.compile('[^a-zA-Z]')
seqs_dict = {sample['covv_virus_name'].replace('hCoV-19/', '').replace(' ', ''): 
             regex.sub('', sample['sequence'].replace('\n', '')) for sample in data}
print(f"Converting to FASTA...")
bs.dict2fasta(seqs_dict, out_fp)
print(f"FASTA output generated and saved in {out_fp}")