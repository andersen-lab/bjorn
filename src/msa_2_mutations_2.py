#!/usr/bin/env python
import os
import sys
import gc
import argparse
import time
import json
import pandas as pd
from path import Path

import bjorn_support as bs
import mutations as bm
import data as bd


# def msa_2_mutations(alignment_filepath, patient_zero, out_filepath, config):
#   # date = config['date']
#   patient_zero = config['patient_zero']
#   data_src = config['data_source']
#   # min_date = config['min_date']
#   # unknown_val = config['unknown_value']
#   # countries_fp = config['countries_fp']
#   # divisions_fp = config['divisions_fp']
#   # locations_fp = config['locations_fp']
#   print(f"Loading alignment file at {alignment_filepath}")
#   t0 = time.time()
#   msa_data = bs.load_fasta(alignment_filepath, is_aligned=True, is_gzip=False)
#   msa_load_time = time.time() - t0
#   print(f"Identifying substitution-based mutations...")
#   t0 = time.time()
#   subs, _ = bm.identify_replacements_per_sample(msa_data,
#                                                 # gisaid_meta,
#                                                 gene2pos=bd.GENE2POS,
#                                                 data_src=data_src,
#                                                 min_seq_len=20000,
#                                                 patient_zero=patient_zero
#                                               #   test=is_test
#                                                 )
#   subs_time = time.time() - t0
#   print(f"Identifying deletion-based mutations...")
#   t0 = time.time()
#   dels, _ = bm.identify_deletions_per_sample(msa_data,
#                                             #  gisaid_meta,
#                                             gene2pos=bd.GENE2POS,
#                                             data_src=data_src,
#                                             min_del_len=1,
#                                             max_del_len=500,
#                                             min_seq_len=20000,
#                                             patient_zero=patient_zero
#                                           #    test=is_test
#                                             )
#   dels_time = time.time() - t0
#   # QC FILTER: remove seqs with >500 nt deletions
#   # dels = dels.loc[dels['del_positions'].str.len()<500]
#   print(subs.shape)
#   print(dels.shape)
#   muts = pd.concat([subs, dels])
#   muts['is_synonymous'] = False
#   muts.loc[muts['ref_aa']==muts['alt_aa'], 'is_synonymous'] = True
#   print(muts.shape)
#   # muts = muts.astype(str) TAKES FOREVER
#   # muts_filename = alignment_filepath.replace('.aligned.fasta', f'_{date}.mutations.csv')
#   muts.to_csv(out_filepath, index=False)
#   del muts, subs, dels
#   gc.collect();
#   print(f"Mutations extracted from {alignment_filepath} and saved in {out_filepath}\n")
#   return 0


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
  out_fp = Path(args.outfp)

  print(f"Loading alignment file at {alignment_filepath}")
  t0 = time.time()
  msa_data = bs.load_fasta(alignment_filepath, is_aligned=True, is_gzip=False)
  msa_load_time = time.time() - t0
  print(f"Identifying substitution-based mutations...")
  t0 = time.time()
  subs, _ = bm.identify_replacements_per_sample(msa_data,
                                                # gisaid_meta,
                                                gene2pos=bd.GENE2POS,
                                                data_src=data_src,
                                                min_seq_len=20000,
                                                patient_zero=patient_zero
                                              #   test=is_test
                                                )
  subs_time = time.time() - t0
  print(f"Identifying deletion-based mutations...")
  t0 = time.time()
  dels, _ = bm.identify_deletions_per_sample(msa_data,
                                            #  gisaid_meta,
                                            gene2pos=bd.GENE2POS,
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
  print(subs.shape)
  print(dels.shape)
  muts = pd.concat([subs, dels])
  muts['is_synonymous'] = False
  muts.loc[muts['ref_aa']==muts['alt_aa'], 'is_synonymous'] = True
  print(muts.shape)
  # muts = muts.astype(str) TAKES FOREVER
  # muts_filename = alignment_filepath.replace('.aligned.fasta', f'_{date}.mutations.csv')
  muts.to_csv(out_fp, index=False)
  print(f"Mutations extracted from {alignment_filepath} and saved in {out_fp}\n")
