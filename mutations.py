import os
import math
import re
import gzip
import numpy as np
import pandas as pd
import more_itertools as mit
from Bio import Seq, SeqIO, AlignIO, Phylo, Align
from bjorn_support import map_gene_to_pos
from data import GENE2POS


def identify_replacements(input_fasta, 
                          meta_fp,
                          patient_zero: str='NC_045512.2', 
                          gene2pos: dict=GENE2POS,
                          location: str=None,
                          data_src: str='alab'):
    print(f"Creating a dataframe...")
    seqsdf, _ = identify_replacements_per_sample(input_fasta, 
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
        # subs['divisions'] = subs['divisions'].apply(process_list).astype(str)
        subs['division_counts'] = subs['division_counts'].apply(lambda x: x[1]).apply(lambda x: ','.join(map(str, x)))
        # subs['division_counts'] = subs['division_counts'].apply(process_list).astype(str)
        subs['countries'] = subs['country_counts'].apply(lambda x: x[0]).apply(lambda x: ','.join(x))#.astype(str)
        # subs['countries'] = subs['countries'].apply(process_list).astype(str)
        subs['country_counts'] = subs['country_counts'].apply(lambda x: x[1]).apply(lambda x: ','.join(map(str, x)))
        # subs['country_counts'] = subs['country_counts'].apply(process_list).astype(str)
    else:
        subs = (seqsdf.groupby(['gene', 'ref_codon', 'alt_codon', 'pos', 'ref_aa', 
                            'codon_num', 'alt_aa'])
                    .agg(
                    num_samples=('idx', 'nunique'),
                    # first_detected=('date', 'min'),
                    # last_detected=('date', 'max'),
                #      locations=('location', uniq_locs),
                    # location_counts=('location', 
                    #                 lambda x: np.unique(x, return_counts=True)),
                    samples=('idx', 'unique')
                    )
                    .reset_index())
        subs['samples'] = subs['samples'].apply(process_samples).astype(str)
    # subs['locations'] = subs['location_counts'].apply(lambda x: x[0]).apply(lambda x: ','.join(x))
    # subs['location_counts'] = subs['location_counts'].apply(lambda x: x[1]).apply(lambda x: ','.join(map(str, x)))
    # 1-based nucleotide position coordinate system
    subs['pos'] = subs['pos'] + 1
    return subs


def process_samples(x):
    x = ','.join(x)
    return x


def identify_replacements_per_sample(input_fasta, 
                                     meta_fp,
                                     gene2pos,
                                     data_src,
                                     patient_zero: str='NC_045512.2',
                                     is_gzip: bool=False,
                                     test: bool=False):
    print(f"Loading Alignment file at: {input_fasta}")
    if is_gzip:
        with gzip.open(input_fasta, "rt") as handle:
            cns = AlignIO.read(handle, 'fasta')
    else:
        cns = AlignIO.read(input_fasta, 'fasta')
    print(f"Initial cleaning...")
    seqs, ref_seq = process_cns_seqs(cns, patient_zero,
                                     start_pos=0, end_pos=30000)
    seqsdf = (pd.DataFrame(index=seqs.keys(), 
                           data=seqs.values(), 
                           columns=['sequence'])
                .reset_index()
                .rename(columns={'index': 'idx'}))
    if test:
        seqsdf = seqsdf.sample(100)
    print(f"Identifying mutations...")
    # for each sample, identify list of substitutions (position:alt)
    seqsdf['replacements'] = seqsdf['sequence'].apply(find_replacements, 
                                                      args=(ref_seq,))
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
    seqsdf = seqsdf.loc[~seqsdf['gene'].isna()]
    # filter our substitutions in non-gene positions
    seqsdf = seqsdf.loc[seqsdf['gene']!='nan']
    print(f"Computing codon numbers...")
    # compute codon number of each substitution
    seqsdf['codon_num'] = seqsdf.apply(compute_codon_num, args=(gene2pos,), axis=1)
    print(f"Fetching reference codon...")
    # fetch the reference codon for each substitution
    seqsdf['ref_codon'] = seqsdf.apply(get_ref_codon, args=(ref_seq, gene2pos), axis=1)
    print(f"Fetching alternative codon...")
    # fetch the alternative codon for each substitution
    seqsdf['alt_codon'] = seqsdf.apply(get_alt_codon, args=(gene2pos,), axis=1)
    # drop the actual sequences to save mem
    seqsdf.drop(columns=['sequence'], inplace=True)
    print(f"Mapping amino acids...")
    # fetch the reference and alternative amino acids
    seqsdf['ref_aa'] = seqsdf['ref_codon'].apply(get_aa)
    seqsdf['alt_aa'] = seqsdf['alt_codon'].apply(get_aa)
    # filter out substitutions with non-amino acid alternates (bad consensus calls)
    seqsdf = seqsdf.loc[seqsdf['alt_aa']!='nan']
    print("Naming substitutions")
    seqsdf['mutation'] = seqsdf['gene'] + ':' + seqsdf['ref_aa'] + seqsdf['codon_num'].astype(str) + seqsdf['alt_aa']
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
            # filter out improper collection dates
            # meta['tmp'] = meta['date'].str.split('-')
            # meta = meta[meta['tmp'].str.len()>=3]
            seqsdf = pd.merge(seqsdf, meta, left_on='idx', right_on='strain')
            # seqsdf['date'] = pd.to_datetime(seqsdf['date'], errors='coerce')
            # seqsdf['month'] = seqsdf['date'].dt.month
            seqsdf.loc[seqsdf['location'].isna(), 'location'] = 'unk'
            seqsdf = seqsdf[seqsdf['host']=='Human']
            seqsdf.loc[seqsdf['country']=='USA', 'country'] = 'United States of America'
        else:
            raise ValueError(f"user-specified data source {data_src} not recognized. Aborting.")
    return seqsdf, ref_seq


def find_replacements(x, ref):
    return [f'{i}:{n}' for i, n in enumerate(x) 
            if n!=ref[i] and n!='-' and n!='n']


def compute_codon_num(x, gene2pos: dict):
    pos = x['pos']
    ref_pos = gene2pos[x['gene']]['start']
    return math.ceil((pos - ref_pos + 1) / 3)


def get_ref_codon(x, ref_seq, gene2pos: dict):
    ref_pos = gene2pos[x['gene']]['start']
    codon_start = ref_pos + ((x['codon_num'] - 1) * 3)
    return ref_seq[codon_start: codon_start+3].upper()


def get_alt_codon(x, gene2pos: dict):
    ref_pos = gene2pos[x['gene']]['start']
    codon_start = ref_pos + ((x['codon_num'] - 1) * 3)
    return x['sequence'][codon_start: codon_start+3].upper()


def get_aa(codon: str):
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
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
    } 
    return CODON2AA.get(codon, 'nan')


def identify_deletions(input_fasta: str, 
                       meta_fp: str,
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
    seqsdf, ref_seq = identify_deletions_per_sample(input_fasta, 
                                           meta_fp,  
                                           patient_zero,
                                           gene2pos,
                                           min_del_len,
                                           start_pos,
                                           end_pos,
                                           data_src)
    if location:
        seqsdf = seqsdf.loc[seqsdf['location'].str.contains(location)]
    # group sample by the deletion they share
    if data_src=='gisaid':
        del_seqs = (seqsdf.groupby(['relative_coords', 'del_len'])
                .agg(
                 num_samples=('idx', 'nunique'),
                #  first_detected=('date', 'min'),
                #  last_detected=('date', 'max'),
                #  num_locations=('location', 'nunique'),
                #  location_counts=('location', 
                #                   lambda x: np.unique(x, 
                #                                       return_counts=True)),
                #  num_divisions=('division', 'nunique'),
                #  division_counts=('division', 
                #                   lambda x: np.unique(x, 
                #                                       return_counts=True)),
                #  num_countries=('country', 'nunique'),
                #  country_counts=('country', 
                #                  lambda x: np.unique(x, 
                #                                      return_counts=True))
                )
                .reset_index())
        # del_seqs['divisions'] = del_seqs['division_counts'].apply(lambda x: x[0]).apply(lambda x: ','.join(x))
        # del_seqs['division_counts'] = del_seqs['division_counts'].apply(lambda x: x[1]).apply(lambda x: ','.join(map(str, x)))
        # del_seqs['countries'] = del_seqs['country_counts'].apply(lambda x: x[0]).apply(lambda x: ','.join(x))
        # del_seqs['country_counts'] = del_seqs['country_counts'].apply(lambda x: x[1]).apply(lambda x: ','.join(map(str, x)))
    else:
        del_seqs = (seqsdf.groupby(['relative_coords', 'del_len'])
                        .agg(samples=('idx', 'unique'),
                             num_samples=('idx', 'nunique'),
                            #  first_detected=('date', 'min'),
                            #  last_detected=('date', 'max'),
                            #  location_counts=('location', 
                            #                   lambda x: np.unique(x, return_counts=True))
                            )
                        .reset_index()
                        .sort_values('num_samples'))
        del_seqs['samples'] = del_seqs['samples'].apply(process_samples).astype(str)
    # del_seqs['locations'] = del_seqs['location_counts'].apply(lambda x: x[0]).apply(lambda x: ','.join(x))
    # del_seqs['location_counts'] = del_seqs['location_counts'].apply(lambda x: x[1]).apply(lambda x: ','.join(map(str, x)))
    del_seqs['type'] = 'deletion'
    # adjust coordinates to account for the nts trimmed from beginning e.g. 265nts
    del_seqs['absolute_coords'] = del_seqs['relative_coords'].apply(adjust_coords, args=(start_pos+1,))
    del_seqs['pos'] = del_seqs['absolute_coords'].apply(lambda x: int(x.split(':')[0]))
    # approximate the gene where each deletion was identified
    del_seqs['gene'] = del_seqs['pos'].apply(map_gene_to_pos)
    del_seqs = del_seqs.loc[~del_seqs['gene'].isna()]
    # filter our substitutions in non-gene positions
    del_seqs = del_seqs.loc[del_seqs['gene']!='nan']
    # compute codon number of each substitution
    del_seqs['codon_num'] = del_seqs.apply(compute_codon_num, args=(gene2pos,), axis=1)
    # fetch the reference codon for each substitution
    del_seqs['ref_codon'] = del_seqs.apply(get_ref_codon, args=(ref_seq, gene2pos), axis=1)
    # fetch the reference and alternative amino acids
    del_seqs['ref_aa'] = del_seqs['ref_codon'].apply(get_aa)
    # record the 5 nts before each deletion (based on reference seq)
    del_seqs['prev_5nts'] = del_seqs['absolute_coords'].apply(lambda x: ref_seq[int(x.split(':')[0])-5:int(x.split(':')[0])])
    # record the 5 nts after each deletion (based on reference seq)
    del_seqs['next_5nts'] = del_seqs['absolute_coords'].apply(lambda x: ref_seq[int(x.split(':')[1])+1:int(x.split(':')[1])+6])
    if data_src=='alab':
        cols = ['type', 'gene', 'absolute_coords', 'del_len', 'pos', 
                     'ref_aa', 'codon_num', 'num_samples',
                     'first_detected', 'last_detected', 'locations',
                     'location_counts', 'samples',
                     'ref_codon', 'prev_5nts', 'next_5nts']
    elif data_src=='gisaid':
        cols = del_seqs.columns
    return del_seqs#[cols]


def identify_deletions_per_sample(input_fasta, meta_fp, 
                                  gene2pos, 
                                  data_src,
                                  min_del_len=1, 
                                  start_pos=265,
                                  end_pos=29674,
                                  patient_zero: str='NC_045512.2',
                                  is_gzip: bool=False,
                                  test=False):
    # read MSA file
    print(f"Loading Alignment file at: {input_fasta}")
    if is_gzip:
        with gzip.open(input_fasta, "rt") as handle:
            cns = AlignIO.read(handle, 'fasta')
    else:
        cns = AlignIO.read(input_fasta, 'fasta')
    # prcess MSA to remove insertions and fix position coordinate systems
    seqs, ref_seq = process_cns_seqs(cns, patient_zero, start_pos, end_pos)
    print(f"Initial cleaning...")
    # load into dataframe
    seqsdf = (pd.DataFrame(index=seqs.keys(), data=seqs.values(), 
                           columns=['sequence'])
                .reset_index().rename(columns={'index': 'idx'}))
    if test:
        seqsdf = seqsdf.sample(100)
    # compute length of each sequence
    seqsdf['seq_len'] = seqsdf['sequence'].str.len()
    print(f"Identifying deletions...")
    # identify deletion positions
    seqsdf['del_positions'] = seqsdf['sequence'].apply(find_deletions)
    # dump sequences to save mem, boost speed
    seqsdf.drop(columns=['sequence'], inplace=True)
    # sequences with one or more deletions
    seqsdf = seqsdf.loc[seqsdf['del_positions'].str.len() > 0]
    seqsdf = seqsdf.explode('del_positions')
    # compute length of each deletion
    seqsdf['del_len'] = seqsdf['del_positions'].apply(len)
    # only consider deletions longer than 2nts
    seqsdf = seqsdf[seqsdf['del_len'] >= min_del_len]
    # fetch coordinates of each deletion
    seqsdf['relative_coords'] = seqsdf['del_positions'].apply(get_indel_coords)
    seqsdf['type'] = 'deletion'
    # adjust coordinates to account for the nts trimmed from beginning e.g. 265nts
    seqsdf['absolute_coords'] = seqsdf['relative_coords'].apply(adjust_coords, args=(start_pos+1,))
    seqsdf['pos'] = seqsdf['absolute_coords'].apply(lambda x: int(x.split(':')[0])+1)
    print(f"Mapping Genes to mutations...")
    # approximate the gene where each deletion was identified
    seqsdf['gene'] = seqsdf['pos'].apply(map_gene_to_pos)
    seqsdf = seqsdf.loc[~seqsdf['gene'].isna()]
    # filter our substitutions in non-gene positions
    seqsdf = seqsdf.loc[seqsdf['gene']!='nan']
    print(f"Computing codon numbers...")
    # compute codon number of each substitution
    seqsdf['codon_num'] = seqsdf.apply(compute_codon_num, args=(gene2pos,), axis=1)
    print(f"Fetching reference codon...")
    # fetch the reference codon for each substitution
    seqsdf['ref_codon'] = seqsdf.apply(get_ref_codon, args=(ref_seq, gene2pos), axis=1)
    print(f"Mapping amino acids...")
    # fetch the reference and alternative amino acids
    seqsdf['ref_aa'] = seqsdf['ref_codon'].apply(get_aa)
    # record the deletion subsequence
    seqsdf['del_seq'] = seqsdf['absolute_coords'].apply(get_deletion, args=(ref_seq,))
    # record the 5 nts before each deletion (based on reference seq)
    seqsdf['prev_5nts'] = seqsdf['absolute_coords'].apply(lambda x: ref_seq[int(x.split(':')[0])-5:int(x.split(':')[0])])
    # record the 5 nts after each deletion (based on reference seq)
    seqsdf['next_5nts'] = seqsdf['absolute_coords'].apply(lambda x: ref_seq[int(x.split(':')[1])+1:int(x.split(':')[1])+6])
    print("Naming deletions")
    seqsdf['mutation'] = seqsdf[['gene', 'codon_num', 'del_len']].apply(assign_deletion, axis=1)
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
            # filter out improper collection dates
            # meta['tmp'] = meta['date'].str.split('-')
            # meta = meta[meta['tmp'].str.len()>=3]
            seqsdf = pd.merge(seqsdf, meta, left_on='idx', right_on='strain')
            # seqsdf['date'] = pd.to_datetime(seqsdf['date'], errors='coerce')
            # seqsdf['month'] = seqsdf['date'].dt.month
            seqsdf.loc[seqsdf['location'].isna(), 'location'] = 'unk'
            seqsdf = seqsdf[seqsdf['host']=='Human']
            seqsdf.loc[seqsdf['country']=='USA', 'country'] = 'United States of America'
        else:
            raise ValueError(f"user-specified data source {data_src} not recognized. Aborting.")
    return seqsdf, ref_seq


def get_deletion(x, ref_seq):
    start, end = int(x.split(':')[0]), int(x.split(':')[1])
    return ref_seq[start+1:end+2]


def assign_deletion(x):
    deletion = x['gene'] + ':DEL' + str(x['codon_num']) + '/' + str(x['codon_num'] + (x['del_len']/3) - 1)
    return deletion


def identify_insertions(input_filepath: str,
                        meta_fp: str,
                        patient_zero: str,
                        gene2pos: dict=GENE2POS,
                        min_ins_len: int=2,
                        start_pos: int=265, 
                        end_pos: int=29674) -> pd.DataFrame:
    """Identify insertions found in the aligned sequences. 
    input_filepath: path to fasta multiple sequence alignment
    patient_zero: name of the reference sequence in the alignment
    min_ins_len: minimum length of insertions to be identified"""
    # read MSA file
    consensus_data = AlignIO.read(input_filepath, 'fasta')
    # load into dataframe
    ref_seq = get_seq(consensus_data, patient_zero)[start_pos:end_pos]
    insert_positions = identify_insertion_positions(ref_seq)
    if insert_positions:
        seqs = get_seqs(consensus_data)
        seqsdf = identify_insertions_per_sample(seqs, 
                                                meta_fp,
                                                insert_positions)
        # sequences with one or more deletions
        ins_seqs = seqsdf.loc[seqsdf['ins_positions'].str.len() > 0]
        ins_seqs = ins_seqs.explode('ins_positions')
        # compute length of each deletion
        ins_seqs['ins_len'] = ins_seqs['ins_positions'].apply(len)
        # only consider deletions longer than 2nts
        ins_seqs = ins_seqs[ins_seqs['ins_len'] >= min_ins_len]
        # fetch coordinates of each deletion
        ins_seqs['relative_coords'] = ins_seqs['ins_positions'].apply(get_indel_coords)
        # group sample by the deletion they share
        ins_seqs = (ins_seqs.groupby(['relative_coords', 'ins_len'])
                            .agg(samples=('ID', 'unique'),       # list of sample IDs with the deletion
                                 num_samples=('ID', 'nunique'),
                                 first_detected=('date', 'min'),
                                 last_detected=('date', 'max'),
    #                              locations=('location', uniq_locs),
                                 location_counts=('location', lambda x: np.unique(x, return_counts=True)))  # num of samples with the deletion
                            .reset_index()
                            .sort_values('num_samples'))
        
        ins_seqs['locations'] = ins_seqs['location_counts'].apply(lambda x: list(x[0]))
        ins_seqs['location_counts'] = ins_seqs['location_counts'].apply(lambda x: list(x[1]))
        # adjust coordinates to account for the nts trimmed from beginning e.g. 265nts
        ins_seqs['absolute_coords'] = ins_seqs['relative_coords'].apply(adjust_coords, args=(start_pos,))
        ins_seqs['pos'] = ins_seqs['absolute_coords'].apply(lambda x: int(x.split(':')[0]))
        # approximate the gene where each deletion was identifiedins_seqs
        ins_seqs['gene'] = ins_seqs['pos'].apply(map_gene_to_pos)
        ins_seqs = ins_seqs.loc[~ins_seqs['gene'].isna()]
        # filter our substitutions in non-gene positions
        ins_seqs = ins_seqs.loc[ins_seqs['gene']!='nan']
        # compute codon number of each substitution
        ins_seqs['codon_num'] = ins_seqs.apply(compute_codon_num, args=(gene2pos,), axis=1)
        # fetch the reference codon for each substitution
        ins_seqs['ref_codon'] = ins_seqs.apply(get_ref_codon, args=(ref_seq, gene2pos), axis=1)
        # fetch the reference and alternative amino acids
        ins_seqs['ref_aa'] = ins_seqs['ref_codon'].apply(get_aa)
        # record the 5 nts before each deletion (based on reference seq)
        ins_seqs['prev_5nts'] = ins_seqs['absolute_coords'].apply(lambda x: ref_seq[int(x.split(':')[0])-5:int(x.split(':')[0])])
        # record the 5 nts after each deletion (based on reference seq)
        ins_seqs['next_5nts'] = ins_seqs['absolute_coords'].apply(lambda x: ref_seq[int(x.split(':')[1])+1:int(x.split(':')[1])+6])
        ins_seqs['type'] = 'insertion'
        # record the 5 nts before each deletion (based on reference seq)
        ins_seqs['prev_5nts'] = ins_seqs['relative_coords'].apply(lambda x: ref_seq[int(x.split(':')[0])-5:int(x.split(':')[0])])
        # record the 5 nts after each deletion (based on reference seq)
        ins_seqs['next_5nts'] = ins_seqs['relative_coords'].apply(lambda x: ref_seq[int(x.split(':')[1])+1:int(x.split(':')[1])+6])
        return ins_seqs[['type', 'gene', 'absolute_coords', 
                         'ins_len', 'pos', 
                         'ref_aa', 'codon_num', 'num_samples',
                         'first_detected', 'last_detected', 'locations',
                         'location_counts', 'samples',
                         'ref_codon', 'prev_5nts', 'next_5nts'
                       ]]
    else:
        return pd.DataFrame()

    
def identify_insertions_per_sample(seqs, meta_fp, insert_positions):
        seqsdf = (pd.DataFrame(index=seqs.keys(), 
                               data=seqs.values(), 
                               columns=['sequence'])
                    .reset_index()
                    .rename(columns={'index': 'idx'}))
        seqsdf['seq_len'] = seqsdf['sequence'].str.len()
        # load and join metadata
        meta = pd.read_csv(meta_fp)
        seqsdf = pd.merge(seqsdf, meta, left_on='idx', right_on='fasta_hdr')
        # clean and process sample collection dates
        seqsdf = seqsdf.loc[(seqsdf['collection_date']!='Unknown') 
                       & (seqsdf['collection_date']!='1900-01-00')]
        seqsdf.loc[seqsdf['collection_date'].str.contains('/'), 'collection_date'] = seqsdf['collection_date'].apply(lambda x: x.split('/')[0])
        seqsdf['date'] = pd.to_datetime(seqsdf['collection_date'])
        seqsdf['ins_positions'] = seqsdf['sequence'].apply(find_insertions, args=(insert_positions,))
        return seqsdf


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


def find_deletions(x):
    del_positions = [m.start() for m in re.finditer('-', x)]
    deletions = [list(deletion) for deletion in mit.consecutive_groups(del_positions)]
    return deletions


def find_insertions(x, insert_positions: list):
    ins_positions = [m for m in insert_positions if x[m]!='-' and x[m]!='n']
    insertions = [list(insert) for insert in mit.consecutive_groups(ins_positions)]
    return insertions


def get_seq(all_seqs: Align.MultipleSeqAlignment, sample_name: str) -> str:
    """Fetches the aligned sequence of a specific sample name"""
    for rec in all_seqs:
        if sample_name in rec.id:
            seq = rec.seq
            break
    return str(seq)


def identify_insertion_positions(ref_seq: str) -> list:
    """helper function to identify positions where '-' was found in a sequence"""
    return [m.start() for m in re.finditer('-', str(ref_seq))]


def remove_insertions(seq: str, positions: list) -> str:
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
    return [m.start() for m in re.finditer('-', x)]
    
    
def find_deletions_old(x):
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
    return x['del_positions_x']==x['del_positions_y']