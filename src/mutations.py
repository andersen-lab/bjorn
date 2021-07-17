import os
import gc
import math
import re
import gzip
import numpy as np
import pandas as pd

import more_itertools as mit
from Bio import Seq, SeqIO, AlignIO, Phylo, Align
from .bjorn_support import map_gene_to_pos
from .data import GENE2POS



def identify_samples_with_suspicious_mutations(substitutions: pd.DataFrame, 
                                               deletions: pd.DataFrame, 
                                               insertions: pd.DataFrame,
                                               nonconcerning_genes: list,
                                               nonconcerning_mutations: list):
    """Returns list of sample IDs that have triggered an insertion, deletion, and/or 
    substitution-based flag, indicating that they require further manual inspection 
    before public release"""
    # nonconcerning_genes = ['5UTR', 'ORF7a', 'ORF7b', 'ORF8', 'ORF10', 'Non-coding region']
    cols = ['type', 'mutation', 'absolute_coords', 'gene', 'ref_codon', 'alt_codon', 'pos', 'ref_aa',
       'codon_num', 'alt_aa', 'num_samples',
       'is_frameshift', 'indel_len', 'indel_seq', 'relative_coords',
       'prev_5nts', 'next_5nts', 'samples']
    try:
        subs_flag = ((substitutions['alt_aa']=='*') & (~substitutions['gene'].isin(nonconcerning_genes)) & (~substitutions['mutation'].isin(nonconcerning_mutations))) 
        sus_subs_ids = substitutions.loc[subs_flag, 'samples'].str.split(',').explode().unique().tolist()
        sus_subs = substitutions.loc[subs_flag]
    except:
        sus_subs_ids = []
        sus_subs = pd.DataFrame()
    try:
        dels_flag = ((deletions['is_frameshift']==True) & (~deletions['gene'].isin(nonconcerning_genes)) & (~deletions['mutation'].isin(nonconcerning_mutations)))
        sus_dels_ids = deletions.loc[dels_flag, 'samples'].str.split(',').explode().unique().tolist()
        sus_dels = deletions.loc[dels_flag]
    except:
        sus_dels_ids = []
        sus_dels = pd.DataFrame()
    try:
        ins_flag = ((insertions['is_frameshift']==True) & (~insertions['gene'].isin(nonconcerning_genes)) & (~insertions['mutation'].isin(nonconcerning_mutations)))
        sus_ins_ids = insertions.loc[ins_flag, 'samples'].str.split(',').explode().unique().tolist()
        sus_inss = insertions.loc[ins_flag]
    except:
        sus_ins_ids = []
        sus_inss = pd.DataFrame()
    sus_subs = sus_subs.loc[:,~sus_subs.columns.duplicated()]
    sus_dels = sus_dels.loc[:,~sus_dels.columns.duplicated()]
    sus_inss = sus_inss.loc[:,~sus_inss.columns.duplicated()]
    sus_mutations = pd.concat([sus_subs, sus_dels, sus_inss])
    try:
        sus_mutations = sus_mutations[cols]
    except:
        pass
    sus_ids = list(set(sus_subs_ids + sus_dels_ids + sus_ins_ids))
    return sus_ids, sus_mutations

def aggregate_replacements(subs: pd.DataFrame, 
                           date: str,
                           data_src: str):
    """Return a dataframe containing all substitution-based mutations (inside open reading frames)
    The data is aggregated on each type of mutation e.g. S:N501Y"""
    if data_src=='gisaid':
        subs.rename(columns={'date': 'date_collected'}, inplace=True)
    subs.loc[:, 'tmp'] = subs['date_collected'].str.split('-').copy()
    subs = subs[subs['tmp'].str.len()>=2].copy()
    subs['date_collected'] = subs['date_collected'].astype(str)
    subs.loc[subs['tmp'].str.len()==2, 'date_collected'] += '-15'
    subs = subs[subs['date_collected']<date].copy()
    subs.loc[:, 'date_collected'] = pd.to_datetime(subs['date_collected'], errors='coerce')
    # aggregate on mutation
    subs_agg = (subs.groupby(['mutation', 'gene', 'ref_codon', 'pos', 'alt_codon', 'ref_aa', 'codon_num', 'alt_aa'])
                .agg(
                 num_samples=('strain', 'nunique'),
                 first_detected=('date_collected', 'min'),
                 last_detected=('date_collected', 'max'),
                 num_locations=('location_normed', 'nunique'),
                 location_counts=('location_normed', 
                                  lambda x: np.unique(x, 
                                                      return_counts=True)),
                 num_divisions=('division_normed', 'nunique'),
                 division_counts=('division_normed', 
                                  lambda x: np.unique(x, 
                                                      return_counts=True)),
                 num_countries=('country_normed', 'nunique'),
                 country_counts=('country_normed', 
                                 lambda x: np.unique(x, 
                                                     return_counts=True))
                )
                .reset_index())
    subs_agg.loc[:, 'divisions'] = subs_agg['division_counts'].apply(lambda x: x[0]).apply(lambda x: ','.join(x))
    subs_agg.loc[:, 'division_counts'] = subs_agg['division_counts'].apply(lambda x: x[1]).apply(lambda x: ','.join(map(str, x)))
    subs_agg.loc[:, 'countries'] = subs_agg['country_counts'].apply(lambda x: x[0]).apply(lambda x: ','.join(x))#.astype(str)
    subs_agg.loc[:, 'country_counts'] = subs_agg['country_counts'].apply(lambda x: x[1]).apply(lambda x: ','.join(map(str, x)))
    return subs_agg


def aggregate_deletions(dels: pd.DataFrame, 
                        date: str,
                        data_src: str):
    """Return a dataframe containing all deletion-based mutations (inside open reading frames)
    The data is aggregated on each type of mutation e.g. S:DEL69/70"""
    if data_src=='gisaid':
        dels.rename(columns={'date': 'date_collected'}, inplace=True)
    dels.loc[:, 'tmp'] = dels['date_collected'].str.split('-').copy()
    dels['date_collected'] = dels['date_collected'].astype(str)
    dels = dels[dels['tmp'].str.len()>=2].copy()
    dels.loc[dels['tmp'].str.len()==2, 'date_collected'] += '-15'
    dels = dels[dels['date_collected']<date].copy()
    dels.loc[:, 'date_collected'] = pd.to_datetime(dels['date_collected'], errors='coerce')
    # aggregate on mutation
    dels_agg = (dels.groupby(['mutation', 'relative_coords', 'del_len'])
                .agg(
                 num_samples=('strain', 'nunique'),
                 first_detected=('date_collected', 'min'),
                 last_detected=('date_collected', 'max'),
                 num_locations=('location_normed', 'nunique'),
                 location_counts=('location_normed', 
                                  lambda x: np.unique(x, 
                                                      return_counts=True)),
                 num_divisions=('division_normed', 'nunique'),
                 division_counts=('division_normed', 
                                  lambda x: np.unique(x, 
                                                      return_counts=True)),
                 num_countries=('country_normed', 'nunique'),
                 country_counts=('country_normed', 
                                 lambda x: np.unique(x, 
                                                     return_counts=True))
                )
                .reset_index())
    dels_agg.loc[:, 'divisions'] = dels_agg['division_counts'].apply(lambda x: x[0]).apply(lambda x: ','.join(x))
    dels_agg.loc[:, 'division_counts'] = dels_agg['division_counts'].apply(lambda x: x[1]).apply(lambda x: ','.join(map(str, x)))
    dels_agg.loc[:, 'countries'] = dels_agg['country_counts'].apply(lambda x: x[0]).apply(lambda x: ','.join(x))
    dels_agg.loc[:, 'country_counts'] = dels_agg['country_counts'].apply(lambda x: x[1]).apply(lambda x: ','.join(map(str, x)))
    return dels_agg


def process_samples(x):
    """Support function for processing SARS-CoV-2 sequence sample names
    Expects naming conventions followed by the Andersen Lab"""
    x = ','.join(x)
    return x


def identify_replacements_per_sample(cns, 
                                     meta_fp=None,
                                     gene2pos: dict=GENE2POS,
                                     data_src: str='gisaid',
                                     min_seq_len=20000,
                                     max_num_subs=5000,
                                     patient_zero: str='NC_045512.2',
                                     ref_path: str='/home/al/data/hcov19/NC045512.fasta',
                                     test: bool=False):
    """Returns dataframe of all substitution-based mutations from a pre-loaded multiple sequence alignment, 
    containing the reference sequence (default: NC_045512.2)
    The data is NOT aggregated, meaning that there will be a record for each observed substitution for each sample"""
    print(f"Initial cleaning...")
    seqs, ref_seq = process_cns_seqs(cns, patient_zero,
                                     start_pos=0, end_pos=29674)
#     ref_seq = get_seq_from_fasta(ref_path)
    seqsdf = (pd.DataFrame(index=seqs.keys(), 
                           data=seqs.values(), 
                           columns=['sequence'])
                .reset_index()
                .rename(columns={'index': 'idx'}))
    if test:
        seqsdf = seqsdf.sample(100)
    try:
        # compute length of each sequence
        seqsdf['seq_len'] = seqsdf['sequence'].str.len()
        # filter out seqs that are too short
        seqsdf = seqsdf[seqsdf['seq_len']>min_seq_len]
        print(f"Identifying mutations...")
        # for each sample, identify list of substitutions (position:alt)
        seqsdf['replacements'] = seqsdf['sequence'].apply(find_replacements, 
                                                        args=(ref_seq,))
        # sequences with one or more substitutions
        seqsdf = seqsdf.loc[seqsdf['replacements'].str.len() > 0]
        seqsdf = seqsdf.loc[seqsdf['replacements'].str.len() < max_num_subs]
        seqs = dict(zip(seqsdf['idx'], seqsdf['sequence']))
        # drop the actual sequences to save mem
        seqsdf.drop(columns=['sequence'], inplace=True)
        gc.collect();
        # wide-to-long data manipulation
        seqsdf = seqsdf.explode('replacements')
        # initialize position column
        seqsdf['pos'] = -1
        # populate position column
        seqsdf.loc[~seqsdf['replacements'].isna(), 'pos'] = (seqsdf.loc[~seqsdf['replacements'].isna(), 'replacements']
        .apply(lambda x: int(x.split(':')[0])))
        # filter out non-substitutions
        seqsdf = seqsdf.loc[seqsdf['pos']!=-1]
        print(f"Mapping Genes to mutations...")
        # identify gene of each substitution
        seqsdf['gene'] = seqsdf['pos'].apply(map_gene_to_pos)
        seqsdf.loc[seqsdf['gene'].isna(), 'gene'] = 'Non-coding region'
        seqsdf.loc[seqsdf['gene']=='nan', 'gene'] = 'Non-coding region'
        # filter our substitutions in non-gene positions
        seqsdf = seqsdf.loc[seqsdf['gene']!='nan']
        print(f"Computing codon numbers...")
        # compute codon number of each substitution
        seqsdf['gene_start_pos'] = seqsdf['gene'].apply(lambda x: gene2pos.get(x, {}).get('start', 0))
        seqsdf['codon_num'] = np.ceil((seqsdf['pos'] - seqsdf['gene_start_pos'] + 1) / 3).astype(int)
        print(f"Fetching reference codon...")
        # fetch the reference codon for each substitution
        seqsdf['codon_start'] = seqsdf['gene_start_pos'] + (3*(seqsdf['codon_num'] - 1))
        seqsdf['ref_codon'] = seqsdf['codon_start'].apply(lambda x: ref_seq[x:x+3].upper())
        print(f"Fetching alternative codon...")
        # fetch the alternative codon for each substitution
        seqsdf['alt_codon'] = seqsdf[['idx', 'codon_start']].apply(get_alt_codon, args=(seqs,), axis=1)
        del seqs
        gc.collect();
        print(f"Mapping amino acids...")
        # fetch the reference and alternative amino acids
        seqsdf['ref_aa'] = seqsdf['ref_codon'].apply(get_aa)
        seqsdf['alt_aa'] = seqsdf['alt_codon'].apply(get_aa)
        # filter out substitutions with non-amino acid alternates (bad consensus calls)
        seqsdf = seqsdf.loc[seqsdf['alt_aa']!='nan']
        print("Naming substitutions")
        seqsdf['mutation'] = seqsdf['gene'] + ':' + seqsdf['ref_aa'] + seqsdf['codon_num'].astype(str) + seqsdf['alt_aa']
        seqsdf['type'] = 'substitution'
        print(f"Fusing with metadata...")
        # load and join metadata
        if meta_fp:
            if data_src=='alab':
                meta = pd.read_csv(meta_fp)
                seqsdf = pd.merge(seqsdf, meta, left_on='idx', right_on='fasta_hdr')
                # clean and process sample collection dates
                seqsdf = seqsdf.loc[(seqsdf['collection_date']!='Unknown') 
                            & (seqsdf['collection_date']!='1900-01-00')]
                seqsdf.loc[seqsdf['collection_date'].str.contains('/'), 'collection_date'] = seqsdf['collection_date'].apply(lambda x: x.split('/')[0])
                seqsdf['date'] = pd.to_datetime(seqsdf['collection_date'])
            elif data_src=='gisaid':
                meta = pd.read_csv(meta_fp, sep='\t', compression='gzip')
                seqsdf = pd.merge(seqsdf, meta, left_on='idx', right_on='strain')
                seqsdf.loc[seqsdf['location'].isna(), 'location'] = 'unk'
                seqsdf = seqsdf[seqsdf['host']=='Human']
                seqsdf.loc[seqsdf['country']=='USA', 'country'] = 'United States of America'
            elif data_src=='gisaid_feed':
                meta = pd.read_csv(meta_fp, sep='\t', compression='gzip')
                seqsdf = pd.merge(seqsdf, meta, left_on='idx', right_on='strain')
                # seqsdf.loc[seqsdf['country']=='USA', 'country'] = 'United States of America'
            else:
                raise ValueError(f"user-specified data source {data_src} not recognized. Aborting.")
    except:
        print(f"No substitutions found in any of the sequences in the alignment. Output will contain an empty dataframe")
        seqsdf = pd.DataFrame()
    return seqsdf, ref_seq


# def identify_out_frame_replacements()


def find_replacements(x, ref):
    "Support function for enumerating nucleotide substitutions"
    # return [f'{i}:{n}' for i, n in enumerate(x) 
    #         if n!=ref[i] and n!='-' and n!='n']
    return [f'{i}:{n}' for i, n in enumerate(x) 
            if n!=ref[i] and n!='n']


def get_ref_codon(x, ref_seq, gene2pos: dict):
    "Support function for fetching the reference codon"
    try:
        ref_pos = gene2pos[x['gene']]['start']
        codon_start = ref_pos + ((x['codon_num'] - 1) * 3)
        return ref_seq[codon_start: codon_start+3].upper()
    except:
        return 'NA'


def get_alt_codon(x, seqs: dict):
    "Support function for fetching the alternative codon"
    seq = seqs[x['idx']]
    codon_start = x['codon_start']
    return seq[codon_start:].replace('-', '')[:3].upper()
    # return seq[codon_start:codon_start+3].upper()
    # except:
    #     return 'NA'


def identify_deletions_per_sample(cns, 
                                  meta_fp=None, 
                                  gene2pos: dict=GENE2POS, 
                                  data_src='gisaid',
                                  min_del_len=1, 
                                  max_del_len=500,
                                  min_seq_len=20000,
                                  start_pos=265,
                                  end_pos=29674,
                                  patient_zero: str='NC_045512.2',
                                  test=False):
    """Returns dataframe of all deletion-based mutations from a pre-loaded multiple sequence alignment, 
    containing the reference sequence (default: NC_045512.2)
    The data is NOT aggregated, meaning that there will be a record for each observed deletion for each sample"""
    seqs, ref_seq = process_cns_seqs(cns, patient_zero, start_pos, end_pos)
    print(f"Initial cleaning...")
    # load into dataframe
    seqsdf = (pd.DataFrame(index=seqs.keys(), data=seqs.values(), 
                           columns=['sequence'])
                .reset_index().rename(columns={'index': 'idx'}))
    if test:
        seqsdf = seqsdf.sample(100)
    try:
        # compute length of each sequence
        seqsdf['seq_len'] = seqsdf['sequence'].str.len()
        seqsdf = seqsdf[seqsdf['seq_len']>min_seq_len]
        print(f"Identifying deletions...")
        # identify deletion positions
        seqsdf['del_positions'] = seqsdf['sequence'].apply(find_deletions)
        # dump sequences to save mem, boost speed
        seqsdf.drop(columns=['sequence'], inplace=True)
        gc.collect();
        # sequences with one or more deletions
        seqsdf = seqsdf.loc[seqsdf['del_positions'].str.len() > 0]
        # sequences with less than 500 deletions
        seqsdf = seqsdf.loc[seqsdf['del_positions'].str.len() < max_del_len]
        seqsdf = seqsdf.explode('del_positions')
        # compute length of each deletion
        seqsdf['del_len'] = seqsdf['del_positions'].apply(len)
        # only consider deletions longer than 1nts
        seqsdf = seqsdf[seqsdf['del_len'] >= min_del_len]
        # only consider deletions shorter than 500nts
        seqsdf = seqsdf[seqsdf['del_len'] < max_del_len]
        # fetch coordinates of each deletion
        seqsdf['relative_coords'] = seqsdf['del_positions'].apply(get_indel_coords)
        seqsdf['type'] = 'deletion'
        # adjust coordinates to account for the nts trimmed from beginning e.g. 265nts
        seqsdf['absolute_coords'] = seqsdf['relative_coords'].apply(adjust_coords, args=(start_pos+1,))
        seqsdf['pos'] = seqsdf['absolute_coords'].apply(lambda x: int(x.split(':')[0])+1)
        print(f"Mapping Genes to mutations...")
        # approximate the gene where each deletion was identified
        seqsdf['gene'] = seqsdf['pos'].apply(map_gene_to_pos)
        seqsdf.loc[seqsdf['gene'].isna(), 'gene'] = 'Non-coding region'
        # filter our substitutions in non-gene positions
        seqsdf.loc[seqsdf['gene']=='nan', 'gene'] = 'Non-coding region'
        print(f"Computing codon numbers...")
        seqsdf['codon_num'] = seqsdf.apply(compute_codon_num, args=(gene2pos,), axis=1)
        print(f"Fetching reference codon...")
        # fetch the reference codon for each substitution
        print(f"Mapping amino acids...")
        # fetch the reference and alternative amino acids
        # record the deletion subsequence
        seqsdf['del_seq'] = seqsdf['absolute_coords'].apply(get_deletion, args=(ref_seq,))
        # record the 5 nts before each deletion (based on reference seq)
        seqsdf['prev_5nts'] = seqsdf['absolute_coords'].apply(lambda x: ref_seq[int(x.split(':')[0])-5:int(x.split(':')[0])])
        # record the 5 nts after each deletion (based on reference seq)
        seqsdf['next_5nts'] = seqsdf['absolute_coords'].apply(lambda x: ref_seq[int(x.split(':')[1])+1:int(x.split(':')[1])+6])
        print("Naming deletions")
        seqsdf['pos'] = seqsdf['absolute_coords'].apply(lambda x: int(x.split(':')[0]))
        seqsdf['ref_codon'] = seqsdf['del_seq'].copy()
        seqsdf['gene_start_pos'] = seqsdf['gene'].apply(lambda x: gene2pos.get(x, {}).get('start', -2)+2)
        seqsdf['pos_in_codon'] = (seqsdf['pos'] - seqsdf['gene_start_pos']) % 3
        seqsdf['mutation'] = seqsdf[['pos_in_codon', 'gene', 'codon_num', 'del_len']].apply(assign_deletion_v2, axis=1)
        seqsdf['deletion_codon_coords'] = seqsdf[['pos_in_codon', 'gene', 'codon_num', 'del_len']].apply(assign_deletion_codon_coords, axis=1)
        seqsdf['is_frameshift'] = seqsdf['del_len'].apply(is_frameshift)
        del seqs
        gc.collect();
        print(f"Fuse with metadata...")
        # load and join metadata
        if meta_fp:
            if data_src=='alab':
                meta = pd.read_csv(meta_fp)
                seqsdf = pd.merge(seqsdf, meta, left_on='idx', right_on='fasta_hdr')
                # clean and process sample collection dates
                seqsdf = seqsdf.loc[(seqsdf['collection_date']!='Unknown') 
                            & (seqsdf['collection_date']!='1900-01-00')]
                seqsdf.loc[seqsdf['collection_date'].str.contains('/'), 'collection_date'] = seqsdf['collection_date'].apply(lambda x: x.split('/')[0])
                seqsdf['date'] = pd.to_datetime(seqsdf['collection_date'])
            elif data_src=='gisaid':
                meta = pd.read_csv(meta_fp, sep='\t', compression='gzip')
                seqsdf = pd.merge(seqsdf, meta, left_on='idx', right_on='strain')
                seqsdf.loc[seqsdf['location'].isna(), 'location'] = 'unk'
                seqsdf = seqsdf[seqsdf['host']=='Human']
                seqsdf.loc[seqsdf['country']=='USA', 'country'] = 'United States of America'
            elif data_src=='gisaid_feed':
                meta = pd.read_csv(meta_fp, sep='\t', compression='gzip')
                seqsdf = pd.merge(seqsdf, meta, left_on='idx', right_on='strain')
                # seqsdf.loc[seqsdf['country']=='USA', 'country'] = 'United States of America'
            else:
                raise ValueError(f"user-specified data source {data_src} not recognized. Aborting.")
    except:
        print(f"No deletions found in any of the sequences in the alignment. Output will contain an empty dataframe")
        seqsdf = pd.DataFrame()
    return seqsdf, ref_seq


def get_deletion(x, ref_seq):
    start, end = int(x.split(':')[0]), int(x.split(':')[1])
    return ref_seq[start+1:end+2]


def assign_deletion(x):
    deletion = x['gene'] + ':DEL' + str(x['codon_num']) + '/' + str(x['codon_num'] + (x['del_len']/3) - 1)
    return deletion

    
def identify_insertions_per_sample(cns, 
                                   meta_fp=None, 
                                   gene2pos: dict=GENE2POS, 
                                   data_src='gisaid',
                                   min_ins_len=1, 
                                   start_pos=265,
                                   end_pos=29674,
                                   patient_zero: str='NC_045512.2',
                                   test=False):
        """Returns dataframe of all insertion-based mutations from a pre-loaded multiple sequence alignment, 
            containing the reference sequence (default: NC_045512.2)
            The data is NOT aggregated, meaning that there will be a record for each observed insertion for each sample"""
        # load into dataframe
        ref_seq = get_seq(cns, patient_zero)[start_pos:end_pos]
        insert_positions = identify_insertion_positions(ref_seq)
        if insert_positions:
            seqs = get_seqs(cns)
        else:
            print(f"No insertions found in any of the sequences in the alignment. Output will contain an empty dataframe")
            return pd.DataFrame(), ref_seq
        seqsdf = (pd.DataFrame(index=seqs.keys(), 
                               data=seqs.values(), 
                               columns=['sequence'])
                    .reset_index()
                    .rename(columns={'index': 'idx'}))
        seqsdf['seq_len'] = seqsdf['sequence'].str.len()
        # identify contiguous insertions 
        seqsdf['ins_positions'] = seqsdf['sequence'].apply(find_insertions, args=(insert_positions,))
        # keep sequences with one or more insertions
        seqsdf = seqsdf.loc[seqsdf['ins_positions'].str.len() > 0]
        # drop sequences to save mem 
        seqsdf.drop(columns=['sequence'], inplace=True)
        seqsdf = seqsdf.explode('ins_positions')
        # compute length of each insertion
        seqsdf['ins_len'] = seqsdf['ins_positions'].apply(len)
        # only consider insertions longer than 2nts
        seqsdf = seqsdf[seqsdf['ins_len'] >= min_ins_len]
        # fetch coordinates of each insertion
        seqsdf['relative_coords'] = seqsdf['ins_positions'].apply(get_indel_coords)
        seqsdf['absolute_coords'] = seqsdf['relative_coords'].apply(adjust_coords, args=(start_pos,))
        # record the deletion subsequence
        # seqsdf['ins_seq'] = seqsdf['absolute_coords'].apply(get_deletion, args=(ref_seq,))
        seqsdf['pos'] = seqsdf['absolute_coords'].apply(lambda x: int(x.split(':')[0])+1)
        # approximate the gene where each insertion was identified
        seqsdf['gene'] = seqsdf['pos'].apply(map_gene_to_pos)
        seqsdf.loc[seqsdf['gene'].isna(), 'gene'] = 'Non-coding region'
        # seqsdf = seqsdf.loc[~seqsdf['gene'].isna()]
        # filter our substitutions in non-gene positions
        seqsdf.loc[seqsdf['gene']=='nan', 'gene'] = 'Non-coding region'
        # compute codon number of each substitution
        seqsdf['codon_num'] = seqsdf.apply(compute_codon_num, args=(gene2pos,), axis=1)
        # fetch the reference codon for each substitution
        seqsdf['ref_codon'] = seqsdf.apply(get_ref_codon, args=(ref_seq, gene2pos), axis=1)
        # fetch the reference and alternative amino acids
        seqsdf['ref_aa'] = seqsdf['ref_codon'].apply(get_aa)
        # start position of the gene that each insert is found on
        seqsdf['gene_start_pos'] = seqsdf['gene'].apply(lambda x: gene2pos.get(x, {}).get('start', -2)+2)
        # insert position in codon number
        seqsdf['pos_in_codon'] = (seqsdf['pos'] - seqsdf['gene_start_pos']) % 3
        # insert mutation name
        seqsdf['mutation'] = seqsdf[['pos_in_codon', 'gene', 'codon_num', 'ins_len']].apply(assign_insertion_v2, axis=1)
        seqsdf['type'] = 'insertion'
        seqsdf['is_frameshift'] = seqsdf['ins_len'].apply(is_frameshift)
        # record the 5 nts before each deletion (based on reference seq)
        seqsdf['prev_5nts'] = seqsdf['relative_coords'].apply(lambda x: ref_seq[int(x.split(':')[0])-5:int(x.split(':')[0])])
        # record the 5 nts after each deletion (based on reference seq)
        seqsdf['next_5nts'] = seqsdf['relative_coords'].apply(lambda x: ref_seq[int(x.split(':')[1])+1:int(x.split(':')[1])+6])
        if meta_fp:
            # load and join metadata
            meta = pd.read_csv(meta_fp)
            seqsdf = pd.merge(seqsdf, meta, left_on='idx', right_on='fasta_hdr')
            # clean and process sample collection dates
            seqsdf = seqsdf.loc[(seqsdf['collection_date']!='Unknown') 
                        & (seqsdf['collection_date']!='1900-01-00')]
            seqsdf.loc[seqsdf['collection_date'].str.contains('/'), 'collection_date'] = seqsdf['collection_date'].apply(lambda x: x.split('/')[0])
            seqsdf['date'] = pd.to_datetime(seqsdf['collection_date'])
        return seqsdf, ref_seq


def pad_aligned_sequences(in_fp, out_fp):
    """helper function that ensures all sequences in the given alignment have equal length using padding"""
    records = SeqIO.parse(in_fp, 'fasta')
    records = list(records) # make a copy, otherwise our generator
                            # is exhausted after calculating maxlen
    maxlen = max(len(record.seq) for record in records)

    # pad sequences so that they all have the same length
    for record in records:
        if len(record.seq) != maxlen:
            sequence = str(record.seq).ljust(maxlen, '.')
            record.seq = Seq.Seq(sequence)
    assert all(len(record.seq) == maxlen for record in records)
    # write to temporary file and do alignment
    with open(out_fp, 'w') as f:
        SeqIO.write(records, f, 'fasta')
    return out_fp
    

def process_cns_seqs(cns_data: Align.MultipleSeqAlignment, patient_zero: str,
                     start_pos: int, end_pos: int) -> (dict, str):
    """Process aligned consensus sequences to prepare them for identifying deletions. 
    The reference sequence is used to identify insertion positions, 
    which are then removed and position numbers are updated."""
    # sequence for patient zero (before removing pseudo deletions)
    ref_seq = get_seq(cns_data, patient_zero)
    # identify insertions (by identifyin 'fake' deletions in the aligned reference sequence)
    insertion_positions = identify_insertion_positions(ref_seq)
    # remove insertions from each sequence to consolidate correct nt positions
    for rec in cns_data:
        rec.seq = remove_insertions(str(rec.seq), insertion_positions)
    # sanity check: ensure that there are no "fake" deletions in reference sequence
    ref_seq = get_seq(cns_data, patient_zero)
    assert not identify_insertion_positions(ref_seq)
    # grab sequences from MSA
    seqs = get_seqs(cns_data, start_pos, end_pos)
    return seqs, ref_seq


def identify_insertion_positions(ref_seq: str) -> list:
    """helper function to identify positions where '-' was found in a sequence"""
    return [m.start() for m in re.finditer('-', str(ref_seq))]


# support functions
def get_seqs(bio_seqs: Align.MultipleSeqAlignment, min_pos: int=265, max_pos: int=29674) -> dict:
    """Parse aligned sequences from Bio.Align.MultipleSeqAlignment to a dict object.
    The keys are sample names and values are their consensus sequences. 
    Each sequence is trimmed from both ends using `min_pos` and `max_pos`"""
    seqs = {}
    for row in bio_seqs:
        sample_name = str(row.id)
        s = str(row.seq)
        seqs[sample_name] = s[min_pos:max_pos]
    return seqs



def compute_codon_num(x, gene2pos: dict):
    "Support function for computing codon number within a gene"
    pos = x['pos']
    try:
        ref_pos = gene2pos[x['gene']]['start']
        return math.ceil((pos - ref_pos + 1) / 3)
    except:
        return 0


def get_alt_codon_old(x, seqs: dict, gene2pos: dict):
    "[DEPRECATED] Support function for fetching the alternative codon"
    try:
        ref_pos = gene2pos[x['gene']]['start']
        codon_start = ref_pos + ((x['codon_num'] - 1) * 3)
        return seqs[x['idx']][codon_start: codon_start+3].upper()
    except:
        return 'NA'


def get_aa(codon: str):
    "Support function for mapping codon to amino acid/stop"
    CODON2AA = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W', 
    } 
    return CODON2AA.get(codon, 'nan')


def find_deletions(x):
    """Support function for identifying contiguous deletions within a sequence"""
    del_positions = [m.start() for m in re.finditer('-', x)]
    deletions = [list(deletion) for deletion in mit.consecutive_groups(del_positions)]
    return deletions


def find_insertions(x, insert_positions: list):
    """Support function for identifying contiguous insertions within a sequence"""
    ins_positions = [m for m in insert_positions if x[m]!='-' and x[m]!='n']
    insertions = [list(insert) for insert in mit.consecutive_groups(ins_positions)]
    return insertions


def get_seq(all_seqs: Align.MultipleSeqAlignment, sample_name: str) -> str:
    """Fetches the aligned sequence of a specific sample name"""
    seq = ''
    for rec in all_seqs:
        if sample_name in rec.id:
            seq = rec.seq
            break
    if len(str(seq))==0:
        print('WARNING: reference sequence not acquired. Something is off.')
    return str(seq)


def get_seq_from_fasta(fasta_filepath):
    """Takes the path of a fasta file (str) containing a single sequence e.g. the reference sequence
    Returns the sequence as string"""
    seq = SeqIO.read(fasta_filepath, 'fasta')
    return str(seq.seq)


def remove_insertions(seq: str, positions: list) -> str:
    """Support function for removing insertions from a sequence, which is useful for 
    correctly identifying deletions and/or substitution-based mutations"""
    for i, pos in enumerate(positions):
#         seqs = seqs[:, :pos-i] + align[:, pos+1-i]
        seq = seq[:pos-i] + seq[pos+1-i:]
    return seq


def get_indel_coords(x):
    """helper function to get deletion coordinates (start:end) from a Pandas Series containing deletion positions"""
    min_pos = np.min(x)
    max_pos = np.max(x)
    return f'{min_pos}:{max_pos}'
    
    
def adjust_coords(x, start_pos: int):
    """helper function to adjust deletion coordinates by adding the number of nucleotides that were trimmed (265nts)"""
    start = int(x.split(':')[0])
    end = int(x.split(':')[1])
    return f'{start+start_pos}:{end+start_pos}'


def find_del_positions(x):
    """Support function for enumerating deletions"""
    return [m.start() for m in re.finditer('-', x)]
    
    
def find_deletions_old(x):
    """[DEPRECATED] Support function for enumerating deletions"""
    del_positions = [m.start() for m in re.finditer('-', x)]
    return [list(map(itemgetter(1), g)) for k, g in groupby(enumerate(del_positions), 
                                                            lambda x: x[0]-x[1])]

def cross_join(df1: pd.DataFrame, df2: pd.DataFrame) -> pd.DataFrame:
    """helper function to perform a cross-join between two dataframes
    Useful for computing pairwise relationships...etc."""
    df1 = df1.assign(key=0)
    df2 = df2.assign(key=0)
    return pd.merge(df1, df2, on='key').drop(columns='key')


def is_deletion_common(x):
    """Support function for deciding whether two given deletions are common"""
    return x['del_positions_x']==x['del_positions_y']


def assign_deletion_codon_coords(x):
    """Support function for assiging the specific codon coordinates (floats) for a given deletion e.g. 69.0/70.0"""
    if (x['pos_in_codon'] + x['del_len']) <= 3:
        return x['gene'] + ':DEL' + str(x['codon_num'] + (x['pos_in_codon']/3))
    deletion = x['gene'] + ':DEL' + str(x['codon_num'] + (x['pos_in_codon']/3)) + '/' + str(x['codon_num'] + (1 + (x['pos_in_codon']/3))  + (x['del_len']/3) - 1)
    return deletion


def assign_deletion_v2(x):
    """Support function for assigning the non-specific codon coordinates (integers) for a given deletion e.g. 69/70"""
    if (x['pos_in_codon'] + x['del_len']) <= 3:
        return x['gene'] + ':DEL' + str(x['codon_num'])
    deletion = x['gene'] + ':DEL' + str(x['codon_num']) + '/' + str(x['codon_num'] + np.maximum(((x['del_len']//3) - 1), 1))
    return deletion


def assign_insertion_v2(x):
    """Support function for assigning the non-specific codon coordinates (integers) for a given deletion e.g. 69/70"""
    if (x['pos_in_codon'] + x['ins_len']) <= 3:
        return x['gene'] + ':INS' + str(x['codon_num'])
    deletion = x['gene'] + ':INS' + str(x['codon_num']) + '/' + str(x['codon_num'] + np.maximum(((x['ins_len']//3) - 1), 1))
    return deletion


def is_frameshift(x):
    """Support function for deciding whether an INDEL is frame-shifting"""
    if x % 3 == 0:
        return False
    return True


def get_dels_separated(x):
    """Support function for separating deletions cuz karthik said so"""
    c1 = int(x[x.find('DEL')+3:].split('/')[0])
    try:
        c2 = int(x.split('/')[1])
        return np.arange(c1, c2+1)
    except:
        return c1


def identify_replacements(cns, 
                          meta_fp,
                          patient_zero: str='NC_045512.2', 
                          gene2pos: dict=GENE2POS,
                          location: str=None,
                          data_src: str='alab'):
    """Returns dataframe of substitution-based mutations from a pre-loaded multiple sequence alignment (fasta),
    which contains the reference sequence (default: NC_045512.2)"""
    print(f"Creating a dataframe...")
    seqsdf, _ = identify_replacements_per_sample(cns, 
                                                 meta_fp,  
                                                 gene2pos,
                                                 data_src=data_src,
                                                 patient_zero=patient_zero)
    if location:
        seqsdf = seqsdf.loc[seqsdf['location'].str.contains(location)]
    # aggregate on each substitutions, compute number of samples and other attributes
    if data_src=='gisaid':
        subs = (seqsdf.groupby(['gene', 'pos', 'ref_aa', 'codon_num', 'alt_aa'])
                .agg(
                 num_samples=('idx', 'nunique'),
                 first_detected=('date', 'min'),
                 last_detected=('date', 'max'),
                 num_locations=('location', 'nunique'),
                 location_counts=('location', 
                                  lambda x: np.unique(x, 
                                                      return_counts=True)),
                 num_divisions=('division', 'nunique'),
                 division_counts=('division', 
                                  lambda x: np.unique(x, 
                                                      return_counts=True)),
                 num_countries=('country', 'nunique'),
                 country_counts=('country', 
                                 lambda x: np.unique(x, 
                                                     return_counts=True))
                )
                .reset_index())
        subs['divisions'] = subs['division_counts'].apply(lambda x: x[0]).apply(lambda x: ','.join(x))
        subs['division_counts'] = subs['division_counts'].apply(lambda x: x[1]).apply(lambda x: ','.join(map(str, x)))
        subs['countries'] = subs['country_counts'].apply(lambda x: x[0]).apply(lambda x: ','.join(x))
        subs['country_counts'] = subs['country_counts'].apply(lambda x: x[1]).apply(lambda x: ','.join(map(str, x)))
    else:
        subs = (seqsdf.groupby(['type', 'mutation', 'gene', 
                                'ref_codon', 'alt_codon', 
                                'pos', 'ref_aa', 
                                'codon_num', 'alt_aa'])
                    .agg(
                    num_samples=('idx', 'nunique'),
                    samples=('idx', 'unique')
                    )
                    .reset_index())
        if data_src=='alab':
            subs['samples'] = subs['samples'].apply(process_samples).astype(str)
    # 1-based nucleotide position coordinate system
    subs['pos'] = subs['pos'] + 1
    return subs


def identify_deletions(cns, 
                       meta_fp=None,
                       patient_zero: str='NC_045512.2',
                       location: str=None,
                       gene2pos: dict=GENE2POS,
                       min_del_len: int=2,
                       start_pos: int=265, 
                       end_pos: int=29674,
                       data_src: str='alab') -> pd.DataFrame:
    """Identify deletions found in the aligned sequences. 
    input_filepath: path to fasta multiple sequence alignment
    patient_zero: name of the reference sequence in the alignment
    min_del_len: minimum length of deletions to be identified"""
    seqsdf, ref_seq = identify_deletions_per_sample(cns, 
                                           meta_fp, 
                                           gene2pos, 
                                           patient_zero=patient_zero,
                                           min_del_len=min_del_len,
                                           start_pos=start_pos,
                                           end_pos=end_pos,
                                           data_src=data_src)
    seqsdf.rename(columns={'del_len': 'indel_len'}, inplace=True)
    seqsdf.rename(columns={'del_seq': 'indel_seq'}, inplace=True)
    if seqsdf.shape[0]==0:
        return seqsdf
    if location:
        seqsdf = seqsdf.loc[seqsdf['location'].str.contains(location)]
    # group sample by the deletion they share
    if data_src=='gisaid':
        del_seqs = (seqsdf.groupby(['relative_coords', 'indel_len'])
                .agg(
                 num_samples=('idx', 'nunique'),
                )
                .reset_index())
    else:
        del_seqs = (seqsdf.groupby(['type', 'mutation', 'absolute_coords', 'is_frameshift', 'gene', 'indel_len', 'indel_seq', 'relative_coords', 'prev_5nts', 'next_5nts'])
                        .agg(
                            num_samples=('idx', 'nunique'),
                            samples=('idx', 'unique'),
                            )
                        .reset_index()
                        .sort_values('num_samples'))
        if data_src=='alab':
            del_seqs['samples'] = del_seqs['samples'].apply(process_samples).astype(str)
    return del_seqs#[cols]


def identify_insertions(cns,
                        meta_fp: str,
                        data_src: str='alab',
                        patient_zero: str='NC_045512.2',
                        gene2pos: dict=GENE2POS,
                        min_ins_len: int=1,
                        start_pos: int=265, 
                        end_pos: int=29674) -> pd.DataFrame:
    """Identify insertions found in the aligned sequences. 
    input_filepath: path to fasta multiple sequence alignment
    patient_zero: name of the reference sequence in the alignment
    min_ins_len: minimum length of insertions to be identified"""
    seqsdf, _ = identify_insertions_per_sample(cns, meta_fp=meta_fp)
    seqsdf.rename(columns={'ins_len': 'indel_len'}, inplace=True)
    if seqsdf.shape[0]>0:
        if data_src=='alab':
            ins_seqs = (seqsdf.groupby(['type', 'mutation', 'absolute_coords', 'is_frameshift', 'gene', 'indel_len', 'relative_coords', 'prev_5nts', 'next_5nts'])
                                .agg(samples=('idx', 'unique'),       # list of sample IDs with the deletion
                                        num_samples=('ID', 'nunique'),
                                        )
                                .reset_index()
                                .sort_values('num_samples'))
            
            ins_seqs['samples'] = ins_seqs['samples'].apply(process_samples).astype(str)
            return ins_seqs
        else:
            ins_seqs = (seqsdf.groupby(['type', 'is_frameshift', 'gene', 'absolute_coords', 'relative_coords', 'indel_len', 'pos', 'codon_num', 'prev_5nts', 'next_5nts'])
                                .agg(samples=('idx', 'unique'),       # list of sample IDs with the deletion
                                     num_samples=('idx', 'nunique'))  # num of samples with the deletion
                                .reset_index()
                                .sort_values('num_samples'))
            return ins_seqs[['type', 'mutation', 'absolute_coords', 'is_frameshift', 'gene', 'indel_len', 'relative_coords',
                            'pos', 'num_samples',
                            'prev_5nts', 'next_5nts', 'samples'
                        ]]
    return seqsdf