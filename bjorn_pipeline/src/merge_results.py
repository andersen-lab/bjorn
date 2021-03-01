import os
import sys
import glob
import argparse
import time
import json
import pandas as pd
from path import Path


# COLLECTING USER PARAMETERS
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputdir",
                        type=str,
                        required=True,
                        help="Input directory containing chunked mutation results in csv files")
parser.add_argument("-o", "--outfp",
                        type=str,
                        required=True,
                        help="Output filepath")

args = parser.parse_args()

# whether or not to include bam files in the release
input_dir = args.inputdir
out_fp = args.outfp
# fetch list of mutation csv filepaths
mut_fps = glob.glob(f"{input_dir}/*.mutations.csv")
# concat with pd
mut_df = pd.concat((pd.read_csv(fp) for fp in mut_fps))
# generate json
meta_info = [
        'strain', 'accession_id',
        'date_modified', 'date_collected', 'date_submitted',
        'country_id', 'country', 'country_original', 'country_lower', 'country_original_lower',
        'division_id', 'division', 'division_original', 'division_lower', 'division_original_lower',
        'location_id', 'location', 'location_original', 'location_lower', 'location_original_lower',
#         'submitting_lab', 'originating_lab',
#         'authors', 
        'pangolin_lineage', 'pangolin_version',
        'clade', 
#         'nextstrain_clade',
#         'gisaid_epi_isl', 'genbank_accession',
#         'purpose_of_sequencing',
            ]

muts_info = ['type', 'mutation', 'gene', 
             'ref_codon', 'pos', 'alt_codon', 
             'is_synonymous', 
             'ref_aa', 'codon_num', 'alt_aa', 
             'absolute_coords', 
             'change_length_nt', 'is_frameshift',
             'deletion_codon_coords']
mut_df = mut_df[~(mut_df['gene'].isin(['5UTR', '3UTR']))]
# GENERATE JSON DATA MODEL
start = time.time()
# (mut_df.groupby(meta_info, as_index=True)
#              .apply(lambda x: x[muts_info].to_dict('records'))
#              .reset_index()
#              .rename(columns={0:'mutations'})
#              .to_json(out_fp,
#                       orient='records',
#                       compression='gzip'))
mut_df.to_csv(out_fp, index=False)
end = time.time()
print(f'Execution time: {end - start} seconds')
# upload to gcloud

# update gcloud access rights

# send auto-slack message about it? (nah, too much)