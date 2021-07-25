import pandas as pd
from path import Path
from shutil import copy, move
from Bio import SeqIO, AlignIO, Phylo
import argparse
import glob
import subprocess
from multiprocessing import Pool
from itertools import repeat
import os
from datetime import datetime as dt
import bjorn_support as bs
import mutations as bm
# from bjorn_support import concat_fasta, align_fasta, compute_tree, map_gene_to_pos, load_fasta
# from mutations import identify_replacements, identify_deletions, identify_insertions
import data as bd
import json


## MAIN

if __name__=="__main__":
    with open('biolab_config.json', 'r') as f:
        config = json.load(f)
    # out_dir = Path(config['release_outdir'])
    # date = config['biolabs_date']
    # meta_fp = config['biolabs_meta']
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--out-dir",
                        type=str,
                        help="Path to folder where results are to be saved")
    parser.add_argument("--date",
                        type=str,
                        help="Date assigned to the sequencing run")
    args = parser.parse_args()
    out_dir = Path(args.out_dir)
    date = args.date
    fasta_hub = config['raw_fasta_hub']
    results_hub = config['mutations_hub']
    num_cpus = config['num_cpus']
    ref_path = config['reference_filepath']
    patient_zero = config['patient_zero']
    # create folders for all results
    if Path.isdir(out_dir):
        rm_cmd = f"rm -r {out_dir}"
        bs.run_command(rm_cmd)
    Path.mkdir(out_dir);
    msa_dir = out_dir/'msa'
    if not Path.isdir(msa_dir):
        Path.mkdir(msa_dir);
    seqs_dir = out_dir/'fa'
    if not Path.isdir(seqs_dir):
        Path.mkdir(seqs_dir);
    print(f"Transferring FASTA from Windows to Linux subsystem")
    transfer_fasta_cmd = f"cp {fasta_hub}/* {seqs_dir}/."
    bs.run_command(transfer_fasta_cmd)
    out_fasta_fp = f"{msa_dir}/{date}_release.fa"
    all_sequences = []
    ref_seq = SeqIO.read(ref_path, 'fasta')
    all_sequences.append(ref_seq)
    fasta_fps = glob.glob(f"{seqs_dir}/*.fasta")
    for fp in fasta_fps:
        rec = SeqIO.read(fp, 'fasta')
        all_sequences.append(rec)
    SeqIO.write(all_sequences, out_fasta_fp, 'fasta')
    msa_fp = out_fasta_fp.split('.')[0] + '_aligned.fa'
    if not Path.isfile(Path(msa_fp)):
        msa_fp = bs.align_fasta(out_fasta_fp, msa_fp, num_cpus=num_cpus);
    # load multiple sequence alignment
    msa_data = bs.load_fasta(msa_fp, is_aligned=True)
    # identify insertions
    insertions = bm.identify_insertions(msa_data, 
                                        # meta_fp=meta_fp, 
                                        patient_zero=patient_zero, 
                                        min_ins_len=1)
    # save insertion results to file
    insertions.to_csv(out_dir/'insertions.csv', index=False)
    # identify substitution mutations
    substitutions = bm.identify_replacements(msa_data,
                                patient_zero=patient_zero)
    # save substitution results to file
    substitutions.to_csv(out_dir/'replacements.csv', index=False)
    # identify deletions
    deletions = bm.identify_deletions(msa_data,
                                    patient_zero=patient_zero,
                                    min_del_len=1)
    # save deletion results to file
    deletions.to_csv(out_dir/'deletions.csv', index=False)
    print(msa_fp)
    print(f"Transferring results from Linux subsystem to Windows results folder")
    transfer_results_cmd = f"cp -r {out_dir} {results_hub}/."
    bs.run_command(transfer_results_cmd)
    print(f"Processing Complete. Please proceed to the Windows mutations hub at {results_hub} to view the mutation files")
