#!/usr/bin/env python
import os
import sys
import gc
import argparse
import time
import json
import numpy as np
import pandas as pd

import bjorn_support as bs
import mutations as bm
from data.mappings import GENE2POS

if __name__=="__main__":
  # COLLECTING USER PARAMETERS
  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--input",
                          type=str,
                          required=True,
                          help="FASTA filepath containing aligned sequences")
  # parser.add_argument("-m", "--meta",
  #                         type=str,
  #                         required=True,
  #                         help="Gzipped TSV filepath containing sequence metadata")
  parser.add_argument("-r", "--patient-zero",
                          type=str,
                          default="NC_045512.2",
                          help="Sample name of reference sequence")
  parser.add_argument("-d", "--data-src",
                        type=str,
                        default="gisaid_feed",
                        help="Data source")
  parser.add_argument("-o", "--outfp",
                          type=str,
                          required=True,
                          help="Output filepath storing mutation information")
  args = parser.parse_args()
  alignment_filepath = args.input
  # gisaid_meta = args.meta
  patient_zero = args.patient_zero
  data_src = args.data_src
  out_fp = args.outfp

  t0 = time.time()
  msa_data = bs.load_fasta(alignment_filepath, is_aligned=True, is_gzip=False)
  msa_load_time = time.time() - t0
  t0 = time.time()
  subs, _ = bm.identify_replacements_per_sample(msa_data,
                                                # gisaid_meta,
                                                gene2pos=GENE2POS,
                                                data_src=data_src,
                                                min_seq_len=20000,
                                                patient_zero=patient_zero
                                              #   test=is_test
                                                )
  subs_time = time.time() - t0
  t0 = time.time()
  dels, _ = bm.identify_deletions_per_sample(msa_data,
                                            #  gisaid_meta,
                                            gene2pos=GENE2POS,
                                            data_src=data_src,
                                            min_del_len=1,
                                            max_del_len=500,
                                            min_seq_len=20000,
                                            patient_zero=patient_zero
                                          #    test=is_test
                                            )
  dels_time = time.time() - t0
  # QC FILTER: remove seqs with >500 nt deletions
  # dels = dels.loc[dels['del_positions'].str.len()<500]
  muts = pd.concat([subs, dels])
  muts['is_synonymous'] = False
  muts.loc[muts['ref_aa']==muts['alt_aa'], 'is_synonymous'] = True
  # muts = muts.astype(str) TAKES FOREVER
  # muts_filename = alignment_filepath.replace('.aligned.fasta', f'_{date}.mutations.csv')

  muts = muts[~(muts['gene'].isin(['5UTR', '3UTR']))]
  # ignore mutations found in non-coding regions
  muts = muts.loc[~(muts['gene']=='Non-coding region')]
  # fuse with metadata
  muts_info = [
      'type', 'mutation', 'gene',
      'ref_codon', 'pos', 'alt_codon',
      'is_synonymous',
      'ref_aa', 'codon_num', 'alt_aa',
      'absolute_coords',
      'change_length_nt', 'is_frameshift',
      'deletion_codon_coords'
  ]
  # If deletions not in chunk add columns
  del_columns = ['is_frameshift', 'change_length_nt', 'deletion_codon_coords', 'absolute_coords']
  muts_columns = muts.columns.tolist()
  for i in del_columns:
      if i not in muts_columns:
          muts[i] = np.nan
  muts = muts.groupby('idx').apply(lambda x: x[muts_info].to_dict('records')).reset_index().rename(columns={0:'mutations'})
  muts['mutations'] = muts['mutations'].map(lambda x: [{k:v for k,v in y.items() if pd.notnull(v)} for y in x])

  muts.to_csv(out_fp, index=False)
  sys.exit(0)
