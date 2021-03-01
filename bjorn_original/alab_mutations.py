import subprocess
import shlex
import json
from path import Path
import pandas as pd
import bjorn_support as bs
import mutations as bm
import data as bd


with open('config.json', 'r') as f:
    config = json.load(f)

date = config['date']
out_dir = Path(config['alab_out_dir'])
ref_fp = Path(config['ref_fasta'])
patient_zero = config['patient_zero']
num_cpus = config['num_cpus']
in_alab_seqs = Path(config['alab_sequences'])
in_alab_meta = Path(config['alab_meta'])
if not Path.isdir(out_dir):
    Path.mkdir(out_dir)
    print(f"Created results directory: {out_dir}")
else:
    print(f"Results directory {out_dir} already exists...Continuing...")
# concatenate all consensus sequences
fa_fp = out_dir/'alab_seqs.fa'
if not Path.isfile(fa_fp):
    fa_fp = bs.concat_fasta(in_alab_seqs, out_dir/'alab_seqs')
print(f"Concatenated all sequences and wrote to {fa_fp}")
# align consensus sequences
msa_fp = Path(fa_fp.split('.')[0] + '_aligned.fa')
if not Path.isfile(msa_fp):
    print(f"Aligning sequences with reference...")
    msa_fp = bs.align_fasta_reference(fa_fp, msa_fp, ref_fp=ref_fp, num_cpus=num_cpus)
print(f"Multiple sequence alignment of A-lab samples with reference saved in {msa_fp}")
# msa2_fp = Path(fa_fp.split('.')[0] + '_aligned_absolute.fa')
# if not Path.isfile(msa2_fp):
#     print(f"Aligning sequences without reference...")
#     msa2_fp = bs.align_fasta(fa_fp, msa2_fp, num_cpus=num_cpus)
# print(f"Multiple sequence alignment of A-lab samples without reference saved in {msa2_fp}")
# Identify substitutions and deletions 
msa_data = bs.load_fasta(msa_fp, is_aligned=True)
subs_wide = bm.identify_replacements(msa_data, in_alab_meta, data_src='alab')
subs_wide_fp = out_dir/f'alab_substitutions_wide_{date}.csv'
subs_wide.sort_values('num_samples', ascending=False).to_csv(subs_wide_fp, index=False)
print(f"Substitution-based mutations of A-lab samples saved in {subs_wide_fp}")
dels_wide = bm.identify_deletions(msa_data, in_alab_meta, data_src='alab')
dels_wide_fp = out_dir/f'alab_deletions_wide_{date}.csv'
dels_wide.sort_values('num_samples', ascending=False).to_csv(dels_wide_fp, index=False)
print(f"Deletion-based mutations of A-lab samples saved in {dels_wide_fp}")