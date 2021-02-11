# JSON 2 FASTA
import re
import json
import gzip
from Bio import SeqIO


with open('config.json', 'r') as f:
    config = json.load(f)

in_fp = config['gisaid_feed']
out_fp = config['gisaid_fasta']
regex = re.compile('[^a-zA-Z]')
print(f"Loading JSON...")
data = [json.loads(line) for line in open(in_fp, 'r')]
print(f"Converting to dict...")
seqs_dict = {sample['covv_virus_name'].replace('hCoV-19/', '').replace(' ', ''): 
             regex.sub('', sample['sequence'].replace('\n', '')) for sample in data}
print(f"Converting to FASTA...")
with open(out_fp, 'w') as f:
    f.write(''.join(f'>{idx}\n{seq}\n' for idx, seq in seqs_dict.items()))
print(f"FASTA output generated and saved in {out_fp}")