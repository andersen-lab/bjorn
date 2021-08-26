import time
import hyperjson
import pandas as pd 


fp = '/valhalla/gisaid/workflow/mutations_2021-03-12.csv'
out_fp = '/home/al/analysis/gisaid/api_data.json'
print(f"Reading data at {fp}")
muts = pd.read_csv(fp)
# generate json
meta_info = [
        'strain', 'accession_id',
        'date_modified', 'date_collected', 'date_submitted',
        'country_id', 'country', 'country_original', 'country_lower', 'country_original_lower',
        'division_id', 'division', 'division_original', 'division_lower', 'division_original_lower',
        'location_id', 'location', 'location_original', 'location_lower', 'location_original_lower',
        'pangolin_lineage', 'pangolin_version',
        'clade', 
            ]

muts_info = ['type', 'mutation', 'gene', 
             'ref_codon', 'pos', 'alt_codon', 
             'is_synonymous', 
             'ref_aa', 'codon_num', 'alt_aa', 
             'absolute_coords', 
             'change_length_nt', 'is_frameshift',
             'deletion_codon_coords']
print(f"Generating JSON...")
start = time.time()
json_data = (muts.groupby(meta_info, as_index=True)
             .apply(lambda x: x[muts_info].to_dict('records'))
             .reset_index()
             .rename(columns={0:'mutations'})
             .to_json(orient='records'))
json_time = time.time() - start
print(f"Writing JSON to {out_fp}")
start = time.time()
with open(out_fp, 'w') as f:
    hyperjson.dump(json_data, f)
io_time = time.time() - start
print(f"JSON Execution time: {json_time} seconds")
print(f"IO Execution time: {io_time} seconds")