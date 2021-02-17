import time
import json
from path import Path


import bjorn_support as bs
import onion_trees as ot
import mutations as bm
import visualize as bv
import reports as br
import data as bd


with open('config.json', 'r') as f:
    config = json.load(f)


out_dir = Path(config['out_dir'])
date = config['date']
gisaid_fasta = config['gisaid_fasta']
gisaid_meta = config['gisaid_meta']
ref_fasta = config['ref_fasta']

# output filenames
fasta_filepath = out_dir/f'sequences_{date}.fasta'
sam_filepath = out_dir/f'sequences_{date}.sam'
alignment_filepath = out_dir/f'sequences_{date}_aligned.fasta'
subs_fp = out_dir/f'subs_long_{date}.csv.gz'
dels_fp = out_dir/f'dels_long_{date}.csv.gz'

# extra configs
data_src = config['data_source']
num_cpus = int(config['num_cpus'])
is_gzip = config['is_gzip']
is_test = config['is_test']
test_fasta = out_dir/'test.fasta'
test_size = 100

if is_gzip:
    cmd = f"gunzip {gisaid_fasta}"
    bs.run_command(cmd)
    if '.gz' in gisaid_fasta:
        gisaid_fasta = gisaid_fasta[:-3]
if is_test:
    gisaid_fasta = bs.sample_fasta(gisaid_fasta, test_fasta, sample_size=test_size)
print(f"STEP 1: Aligning sequences...")
if not Path.isfile(fasta_filepath):
    fasta_filepath = bs.concat_fasta_2([gisaid_fasta, ref_fasta], fasta_filepath)
    print(f"Reference sequence added to input sequences and saved at {fasta_filepath}")
t0 = time.time()
if not Path.isfile(sam_filepath):
    t0 = time.time()
    sam_filepath = bs.run_minimap2(fasta_filepath, sam_filepath, ref_fasta, num_cpus=num_cpus)
    print(f"SAM data generated from consensus sequences and saved at {sam_filepath}")
minimap_time = time.time() - t0
print(f"Generating alignment from SAM data...")
t0 = time.time()
if not Path.isfile(alignment_filepath+'.gz'):
    alignment_filepath = bs.run_datafunk(sam_filepath, ref_fasta, alignment_filepath)
    print(f"Alignment generated and saved at {alignment_filepath} \n")
datafunk_time = time.time() - t0
print(f"STEP 2: Counting variants...")
print(f"Loading alignment file at {alignment_filepath}")
t0 = time.time()
# cmd = f"gzip {alignment_filepath}"
# bs.run_command(cmd)
alignment_filepath += '.gz'
msa_data = bs.load_fasta(alignment_filepath, is_aligned=True, is_gzip=True)
msa_load_time = time.time() - t0
print(f"Identifying substitution-based mutations...")
t0 = time.time()
subs, _ = bm.identify_replacements_per_sample(msa_data, 
                                              gisaid_meta,  
                                              bd.GENE2POS, 
                                              data_src=data_src,
                                              test=is_test)
subs_time = time.time() - t0
subs.to_csv(subs_fp, index=False, compression='gzip')
print(f"Identifying deletion-based mutations...")
t0 = time.time()
dels, _ = bm.identify_deletions_per_sample(msa_data, 
                                           gisaid_meta,  
                                           bd.GENE2POS, 
                                           data_src=data_src, 
                                           min_del_len=1,
                                           test=is_test)
dels_time = time.time() - t0
dels.to_csv(dels_fp, index=False, compression='gzip')
print(f"Substitutions saved in {subs_fp}")
print(f"Deletions saved in {dels_fp}")
total_time = minimap_time + datafunk_time + msa_load_time + subs_time + dels_time
# Data logging
with open(out_dir/'log.txt', 'w') as f:
    f.write(f"minimap2 Execution Time: {(minimap_time / 60):.2f} minutes\n")
    f.write(f"datafunk Execution Time: {(datafunk_time / 60):.2f} minutes\n")
    f.write(f"Alignment loading Execution Time: {(msa_load_time / 60):.2f} minutes\n")
    f.write(f"Counting Substitutions Execution Time: {(subs_time / 60):.2f} minutes\n")
    f.write(f"Counting Deletions Execution Time: {(dels_time / 60):.2f} minutes\n")
    f.write(f"Total Execution Time: {(total_time / 60):.2f} minutes\n")
print(f"""END:\tprocess complete. 
        \n\tlogging information saved at {out_dir/'log.txt'}
        \n\tcontact gkarthik@scripps.edu for any issues ;) """)