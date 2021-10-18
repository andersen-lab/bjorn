#!/usr/bin/env python
import os
import gzip
import json
import argparse
import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--inputmetadata",
                        type=str,
                        required=True,
                        help="Input filepath metadata csv")
parser.add_argument("-o", "--outfp",
                        type=str,
                        required=True,
                        help="Output filepath")
parser.add_argument("-f", "--filenames",
                        required=True,
                        nargs='+',
                        help="Filenames of interest")


parser.add_argument("-j", "--jobid",
                        required=True,
                        help="Filenames of interest")



print("Opening metadata")
args = parser.parse_args()
metadata_filepath = args.inputmetadata
output_filepath = args.outfp
#read in the overall metadata file
meta_df = pd.read_csv(metadata_filepath)
all_filenames = args.filenames

#get all .fasta files names
all_fasta = [os.path.basename(os.path.splitext(filename)[0]).replace('.fasta','') for filename in all_filenames]

#get all header strains
job_num =args.jobid
strains = ['None'] * len(all_fasta)
for i,fasta in enumerate(all_fasta):
    with gzip.open(os.path.join(output_filepath,fasta +'.fasta.gz'), "r") as file:
        first_line = file.readline().decode("utf-8").strip().replace('>','')
    strains[i] = str(first_line).strip()

meta_df = meta_df[meta_df['fasta_hdr'].isin(strains)]

c = []
d = []
l = []

#drop and add columns
meta_df.drop(['zipcode', 'authors', 'originating_lab', 'ID','gb_accession', \
'percent_coverage_cds', 'avg_depth'], axis=1, inplace=True)
meta_df['date_submitted'] = ['']*len(meta_df)
meta_df['pangolin_lineage'] = ['']*len(meta_df)
meta_df['pangolin_version'] = ['']*len(meta_df)
meta_df['clade'] = ['']*len(meta_df)
meta_df.rename(columns={'collection_date':'date_collected', \
    'gisaid_accession':'accession_id', 'fasta_hdr':'strain', 'location':'locstring'}, inplace=True)
print(meta_df.columns)
meta_df.to_csv(os.path.join(output_filepath, "%s.tsv.gz" %job_num), "\t")
#concatenate all these things together
all_fasta_files = [os.path.join(output_filepath, fasta+'.fasta.gz') for fasta in all_fasta]
with gzip.open(os.path.join(output_filepath,str(job_num) + '.fasta.gz'), 'w') as outfile:
    for fname in all_fasta_files:
        with gzip.open(fname, "r") as infile:
            for line in infile:
                outfile.write(line)
        outfile.write(b'\n')
        os.system("rm %s" %fname)
