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
    parser.add_argument("-m", "--metadata",
                        type=str,
                        help="Path to file containing biolabs metadata")
    parser.add_argument("-c", "--coverage", 
                        default=80., type=float,
                        help="Minimum threshold for percentage coverage to accept/reject each sequence")
    parser.add_argument("-d", "--depth", 
                        default=100., type=float,
                        help="Minimum threshold for average depth per nucleotide position to accept/reject each sequence")
    args = parser.parse_args()
    out_dir = Path(args.out_dir)
    date = args.date
    meta_fn = args.metadata
    min_coverage = args.coverage
    min_depth = args.depth
    submitter = config['submitter']
    fasta_hub = config['fasta_hub']
    meta_hub = config['meta_hub']
    results_hub = config['results_hub']
    num_cpus = config['num_cpus']
    ref_path = config['reference_filepath']
    patient_zero = config['patient_zero']
    nonconcerning_genes = config['nonconcerning_genes']
    nonconcerning_mutations = config['nonconcerning_mutations']
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
    print(f"Transferring metadata from Windows to Linux subsystem")
    transfer_meta_cmd = f"cp {meta_hub}/* {out_dir}/{date}_gisaid_metadata_release.xls"
    bs.run_command(transfer_meta_cmd)
    accepted_seqs_dir = Path(out_dir/'fa_accepted')
    if not Path.isdir(accepted_seqs_dir):
        Path.mkdir(accepted_seqs_dir);
    intro = pd.read_excel(f'{out_dir}/{date}_gisaid_metadata_release.xls', sheet_name='Instructions')
    cov = pd.read_excel(f'{out_dir}/{date}_gisaid_metadata_release.xls', sheet_name='Coverage')
    meta = pd.read_excel(f'{out_dir}/{date}_gisaid_metadata_release.xls', sheet_name='Submissions', skiprows=1)
    meta = meta.dropna(subset=["FASTA filename"])
    qc_filter = (cov['pct_coverage']>min_coverage) & (cov['avg_depth']>min_depth)
    accepted_samples = cov.loc[qc_filter, 'Biolab Trans. #'].tolist()
    meta['sample_id'] = meta['FASTA filename'].apply(lambda x : x.split('.')[0]).astype(int)
    accepted_sample_filenames = meta.loc[meta['sample_id'].isin(accepted_samples), 'FASTA filename'].tolist()
    meta = meta.loc[meta['sample_id'].isin(accepted_samples)].drop(columns=['sample_id'])
    meta[['Gender', 'Patient age', 'Patient status']] = 'N/A'
    meta['Submitter'] = submitter
    meta.to_csv(f'{out_dir}/{date}_gisaid_metadata_raw.csv', index=False)
    for sample_filename in accepted_sample_filenames:
        copy(f'{seqs_dir}/{sample_filename}', accepted_seqs_dir)
    # copy(ref_path, accepted_seqs_dir)
    fasta_fps = glob.glob(f"{accepted_seqs_dir}/*.fasta")
    out_fasta_fp = f"{msa_dir}/{date}_release.fa"
    all_sequences = []
    ref_seq = SeqIO.read(ref_path, 'fasta')
    all_sequences.append(ref_seq)
    for fp in fasta_fps:
        rec = SeqIO.read(fp, 'fasta')
        rec.id = meta.loc[meta['FASTA filename']==os.path.basename(fp), 'Virus name'].values[0]
        rec.description = ""
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
                                # meta_fp=meta_fp,
                                patient_zero=patient_zero)
    # save substitution results to file
    substitutions.to_csv(out_dir/'substitutions.csv', index=False)
    # identify deletions
    deletions = bm.identify_deletions(msa_data,
                                    # meta_fp=meta_fp,
                                    patient_zero=patient_zero,
                                    min_del_len=1)
    # save deletion results to file
    deletions.to_csv(out_dir/'deletions.csv', index=False)
    # identify samples with suspicious INDELs and/or substitutions
    sus_ids, sus_muts = bm.identify_samples_with_suspicious_mutations(substitutions, 
                                                                        deletions, 
                                                                        pd.DataFrame(),
                                                                        nonconcerning_genes,
                                                                        nonconcerning_mutations)
    sus_muts.to_csv(out_dir/'suspicious_mutations.csv', index=False)
    print(msa_fp)
    print(f"Transferring metadata from Windows to Linux subsystem")
    transfer_results_cmd = f"cp -r {out_dir} {results_hub}/."
    bs.run_command(transfer_results_cmd)
    results_dirname = out_dir.split('/')[-1]
    while True:
        clean_msa_fp_str = f"{results_hub}/{results_dirname}/msa/{date}_release_aligned_clean.fasta"
        data = input(f"Please clean the alignment file and save as {clean_msa_fp_str} inside the msa folder")
        clean_msa_filepaths = glob.glob(f"{results_hub}/{results_dirname}/msa/{date}_release_aligned_clean*")
        if clean_msa_filepaths:
            clean_msa_fp = clean_msa_filepaths[0]
            aln_cmd = f"bash scripts/convert_to_unaligned_fasta.sh {clean_msa_fp} > {clean_msa_fp.replace('_aligned', '')}"
            bs.run_command((aln_cmd))
            break
        else:
            print("I was not able to find the file containing the cleaned alignment. Please ensure that it exists.")
            continue
            
    
    transfer_clean_fasta_cmd = f"cp {clean_msa_fp} {results_hub}/{results_dirname}/msa/."
    bs.run_command(transfer_clean_fasta_cmd)
    print(f"Processing Complete. Please process to the final manual upload steps")
