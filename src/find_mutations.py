import argparse
import pandas as pd
from path import Path
from shutil import copy
from Bio import Seq, SeqIO, AlignIO, Phylo, Align
import bjorn_support as bs
import mutations as bm




parser = argparse.ArgumentParser()

parser.add_argument("-o", "--out-dir",
                        type=str,
                        default="./",
                        help="Output directory")

parser.add_argument("-f", "--fasta-dir",
                    type=str,
                    default="./fa/",
                    help="Directory containing FASTA files")

parser.add_argument("-n", "--cpus",
                    type=int,
                    default=4,
                    help="Number of cpus to use")

parser.add_argument("-r", "--reference",
                    type=str,
                    default="/home/gk/code/hCoV19/db/NC045512.fasta",
                    help="Reference to use")

parser.add_argument("-rn", "--reference-name",
                    type=str,
                    default="NC_045512.2",
                    help="FASTA header name for reference")

args = parser.parse_args()
out_dir = Path(args.out_dir)
seqs_dir = Path(args.fasta_dir)
num_cpus = args.cpus 
ref_path = args.reference
patient_zero = args.reference_name

# fix fasta header names
fa_fps = bs.get_filepaths(seqs_dir, data_fmt='fasta', data_type='', tech='')
for fa_fp in fa_fps['IonXpress']:
    seq = SeqIO.read(fa_fp, 'fasta')
    seq.id = fa_fp.split('/')[-1].split('.')[0]
    SeqIO.write(seq, fa_fp, 'fasta')
msa_dir = out_dir/'msa'
if not Path.isdir(msa_dir):
    Path.mkdir(msa_dir)
# copy reference sequence into folder containing fasta files
copy(ref_path, seqs_dir);
seqs_fp = bs.concat_fasta(seqs_dir, msa_dir/out_dir.basename())
msa_fn = seqs_fp.split('/')[-1].split('.')[0] + '_aligned.fa'
msa_fp = Path(seqs_fp).parent/msa_fn
if not Path.isfile(Path(msa_fp)):
    msa_fp = bs.align_fasta(seqs_fp, msa_fp, num_cpus=num_cpus);
# load multiple sequence alignment
msa_data = bs.load_fasta(msa_fp, is_aligned=True)
muts_dir = out_dir/'mutations'
if not Path.isdir(muts_dir):
    Path.mkdir(muts_dir)
# identify insertions
insertions = bm.identify_insertions(msa_data, 
                                    meta_fp=None, 
                                    data_src='jordan',
                                    patient_zero=patient_zero, 
                                    min_ins_len=1)
# save insertion results to file
insertions.to_csv(muts_dir/'insertions.csv', index=False)
# identify substitution mutations
subs = bm.identify_replacements(msa_data,
                            meta_fp=None,
                            data_src='jordan',
                            patient_zero=patient_zero)
# save substitution results to file
subs.to_csv(muts_dir/'replacements.csv', index=False)
# identify deletions
deletions = bm.identify_deletions(msa_data,
                                  meta_fp=None,
                                  data_src='jordan',
                                  patient_zero=patient_zero,
                                  min_del_len=1)
# save deletion results to file
deletions.to_csv(muts_dir/'deletions.csv', index=False)
# Data logging
with open("{}/data_release.log".format(out_dir), 'w') as f:
    f.write(f"Prepared {len(list(msa_data))-1} samples for review\n")
    f.write(f'Multiple sequence alignment saved in {msa_dir}\n')
    f.write(f'Mutation tables saved in {muts_dir}\n')
    f.write(f'Individual fasta files saved in {seqs_dir}\n')
print(f"Transfer Complete. All results saved in {out_dir}")