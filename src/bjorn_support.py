import sys
import shutil
import glob
import subprocess
import math
from collections import defaultdict
import gzip
from path import Path
import numpy as np
import pandas as pd
from Bio import Seq, SeqIO, AlignIO, Phylo, Align



def get_variant_counts(analysis_filepath: str, search_ids: list, 
                       min_freq: float=0.15, max_freq: float=100.,
                       generalised: bool=False):
    """Takes path to results directory and (optional) list of sample IDs
    Returns data frame of the number of variants identified in each sample 
    with frequencies between `min_freq` and `max_freq` (Default are 0. and 100., respectively)"""
    var_fps = get_filepaths(analysis_filepath, data_fmt='tsv', data_type='variants/illumina', 
                           sample_ids=search_ids, generalised=generalised, return_type='list')
    var_df = pd.concat([pd.read_csv(fp, sep='\t').assign(sample_id=fp.split('/')[-1].split('_')[0]) for fp in var_fps])
    var_df = var_df.loc[(var_df['ALT_FREQ']<max_freq) & (var_df['ALT_FREQ']>min_freq)]
    var_df['variant'] = var_df['POS'].astype(str) + ':' + var_df['ALT'].astype(str)
    return var_df.groupby('sample_id').agg(num_nt_mutations=('variant', 'nunique')).sort_index()


def separate_samples(sample_sheet: pd.DataFrame, 
                     private_ids: list=['LASV'], 
                     primer_ids: list=['iTru']) -> dict:
    """RELEASE functionality: Separate A-lab Sample Sheet into constituent sub-sheets for each primer set."""
    results = {}
    samples_to_exclude = []
    # separate out private sample IDs
    results['_private_sheet'] = pd.concat([sample_sheet.loc[sample_sheet['Sample_ID'].str.contains(idx)] for idx in private_ids])
    for sample_idx in results['_private_sheet']['Sample_ID'].unique():
            samples_to_exclude.append(sample_idx)
    for primer_idx in primer_ids:
        results[f'{primer_idx}_primers'] = sample_sheet.loc[sample_sheet['I7_Index_ID'].str.contains(primer_idx)]
        for sample_idx in results[f'{primer_idx}_primers']['Sample_ID'].unique():
            samples_to_exclude.append(sample_idx)
    results['OG_primers'] = sample_sheet.loc[~sample_sheet['Sample_ID'].isin(samples_to_exclude)]
    return results


def copy_files(filenames: list, destination_dir):
    """Takes list of file paths and copies them to the destination folder"""
    try:
        for fn in filenames:
            shutil.copy(fn, destination_dir)
        return 0
    except:
        raise ValueError("Copy operation failed. Check if all filenames and directories exist")


def generate_release_report(out_dir, report_name='release_report.tar'):
    report_fp = out_dir/report_name
    s1 = "find {out_dir}/ -name '*.csv' -type f -exec tar rfP {report_fp} {{}} \\;".format(out_dir=out_dir, report_fp=report_fp)
    print(s1)
    run_command(s1)
    s2 = f"tar rfP {report_fp} {out_dir}/msa/*"
    print(s2)
    run_command(s2)
    return 0


def separate_alignments(msa_data, sus_ids, out_dir, filename, patient_zero='NC_045512.2'):
    good_seqs = []
    poor_seqs = []
    for rec in msa_data:
        if rec.id in sus_ids:
            poor_seqs.append(rec)
        elif rec.id==patient_zero:
            good_seqs.append(rec)
            poor_seqs.append(rec)
        else:
            good_seqs.append(rec)
    good_msa = Align.MultipleSeqAlignment(good_seqs)   
    good_msa_fn = filename + '_aligned_white.fa'
    good_msa_fp = out_dir/good_msa_fn
    AlignIO.write(good_msa, good_msa_fp, 'fasta')
    poor_msa = Align.MultipleSeqAlignment(poor_seqs) 
    poor_msa_fn = filename + '_aligned_inspect.fa'
    poor_msa_fp = out_dir/poor_msa_fn 
    AlignIO.write(poor_msa, poor_msa_fp, 'fasta')
    return 0
    


def create_chunk_names(meta_filepath: str, chunk_size: int) -> pd.DataFrame:
    # with open(fasta_filepath, 'r') as filehandle:
    #     count = 0
    #     for line in filehandle:
    #         if line.startswith(">"):
    #             count += 1
    meta_df = pd.read_csv(meta_filepath, sep='\t', compression='gzip')
    num_sequences = meta_df['strain'].unique().shape[0]
    chunk_names = [f"chunk_{i+1}" for i in range(math.ceil(num_sequences / chunk_size))]
    return pd.DataFrame(data=chunk_names, columns=['chunk_names'])



def batch_iterator(iterator, chunk_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    Citation: https://biopython.org/wiki/Split_large_file
    """
    record = True
    while record:
        chunk = []
        while len(chunk) < chunk_size:
            try:
                record = next(iterator)
            except StopIteration:
                record = None
            if record is None:
                # End of file
                break
            chunk.append(record)
        if chunk:
            yield chunk


def dict2fasta(seqs: dict, fasta_fp: str, wrap=80):
    with open(fasta_fp, 'w') as f:
        for gid, gseq in seqs.items():
            f.write('>{}\n'.format(gid))
            for i in range(0, len(gseq), wrap):
                f.write('{}\n'.format(gseq[i:i + wrap])) 
    return 0


def integrate_gisaid_meta(old_meta_fp, xtra_fp, msa_fp, rename_cols, drop_cols):
    xtra = pd.read_csv(xtra_fp, sep='\t')
    xtra.rename(columns=rename_cols, inplace=True)
    # xtra.columns
    xtra['country'] = xtra['Location'].apply(lambda x: x.split('/')[1].strip())
    xtra['division'] = xtra['Location'].apply(lambda x: x.split('/')[2].strip())
    xtra['location'] = xtra['Location'].apply(clean_locs)
    xtra['strain'] = xtra['strain'].apply(lambda x: '/'.join(x.split('/')[1:]))
    sois = xtra['strain'].unique().tolist()
    sois = dict.fromkeys(sois)
    meta = pd.read_csv(old_meta_fp, sep='\t', compression='gzip')
    # meta.columns
    meta = meta.loc[~((meta['country'].str.contains('USA'))&(meta['pangolin_lineage']=='B.1.1.7'))]
    meta = pd.concat([meta, xtra])
    meta.drop(columns=drop_cols, inplace=True)
    meta.loc[meta['country']=='USA', 'country'] = 'United States of America'
    return meta


def clean_locs(x):
    if len(x.split('/')) > 3:
        return x.split('/')[3].strip()
    return 'unk'


def fetch_seqs(seqs_filepath, out_fp, sample_idxs: list, is_aligned=False, is_gzip=False):
    if is_aligned:
        if is_gzip:
            with gzip.open(seqs_filepath, "rt") as handle:
                cns = AlignIO.read(handle, 'fasta')
        else:
            cns = AlignIO.read(seqs_filepath, 'fasta')
        my_cns = Align.MultipleSeqAlignment([rec for rec in cns if rec.id in sample_idxs])
        return AlignIO.write(my_cns, out_fp, 'fasta')
    else:
        if is_gzip:
            with gzip.open(seqs_filepath, "rt") as handle:
                cns = SeqIO.parse(handle, 'fasta')
        else:
            cns = SeqIO.parse(seqs_filepath, 'fasta')
        my_cns = [rec for rec in cns if rec.id in sample_idxs]
        return SeqIO.write(my_cns, out_fp, 'fasta')


def get_filepaths(analysis_path: str, data_fmt: str, sample_ids: list=[], 
                  data_type: str='', tech: str='', 
                  generalised=False,
                  return_type='dict') -> dict:
    """Take list of sample IDs, the general area where your files are located, and their format.
       Returns a dictionary with sample IDs as keys and the filepaths as values.
    """
    file_paths = defaultdict(list)
    fs = []
    generalised_flag = ''
    if sample_ids:
        for s_id in sample_ids:
            if generalised:
                generalised_flag = '/**'
            if data_type and tech:
                f = glob.glob(f"{analysis_path}{generalised_flag}/*{data_type}*/*{tech}*/*{s_id}*.{data_fmt}")
            elif data_type:
                f = glob.glob(f"{analysis_path}{generalised_flag}/*{data_type}*/*{s_id}*.{data_fmt}")
            else:
                f = glob.glob(f"{analysis_path}{generalised_flag}/*{s_id}*.{data_fmt}")
            try:
                fs.append(f)
                file_paths[s_id].append(f[0])
            except:
                continue
    else:
        if data_type and tech:
            if generalised:
                fs = glob.glob(f"{analysis_path}/**/*{data_type}*/*{tech}*/*.{data_fmt}")
            else:
                fs = glob.glob(f"{analysis_path}/*{data_type}*/*{tech}*/*.{data_fmt}")
        elif data_type:
            if generalised:
                fs = glob.glob(f"{analysis_path}/**/*{data_type}*/*.{data_fmt}")
            else:
                fs = glob.glob(f"{analysis_path}/*{data_type}*/*.{data_fmt}")
        else:
            if generalised:
                fs = glob.glob(f"{analysis_path}/**/*.{data_fmt}")
            else:
                fs = glob.glob(f"{analysis_path}/*.{data_fmt}")
        for f in fs:
            try:
                sample_id = f.split('/')[-1].split('_')[0]
            except:
                sample_id = f.split('/')[-1]
            file_paths[sample_id].append(f)
    if return_type=='dict':
        return file_paths
    elif return_type=='list':
        return flatten_list(fs)


def get_variant_filepaths(sample_ids: list, analysis_path: str='/home/gk/analysis') -> dict:
    """Takes list of sample IDs and returns filepaths of variant data for corresponding samples"""
    variant_paths = {}
    for s_id in sample_ids:
        f = glob.glob(f"{analysis_path}/**/variants/illumina/*{s_id}*.tsv")
        variant_paths[s_id] = f[0]
    return variant_paths


def find_loc(f: str):
    """helper function to fetch location from filepath"""
    return f.split('/')[-1].split('_')[0].split('-')[-1]


def get_variant_data(variant_filepaths: dict):
    """Takes dict of variant filepaths and loads all variant data into dataframe"""
    df = (pd.concat((pd.read_csv(f, sep='\t')
                     .assign(sample=s_id, location=find_loc(f)) for s_id, f in variant_filepaths.items())))
    return df


def load_fasta(fasta_filepath, is_gzip=False, is_aligned=False):
    if is_gzip:
        with gzip.open(fasta_filepath, "rt") as handle:
            if is_aligned:
                cns = AlignIO.read(handle, 'fasta')
            else:
                cns = SeqIO.parse(handle, 'fasta')
    else:
        if is_aligned:
            cns = AlignIO.read(fasta_filepath, 'fasta')
        else:
            cns = SeqIO.parse(fasta_filepath, 'fasta')
    return cns


def run_datafunk(in_filepath, ref_path, out_filepath):
    df_cmd = f"datafunk sam_2_fasta -s {in_filepath} -r {ref_path} -o {out_filepath} --pad --log-inserts"
    run_command(df_cmd)
    return out_filepath


def run_minimap2(in_filepath, out_filepath, ref_path, num_cpus=25):
    map_cmd = f"minimap2 -a -x asm5 -t {num_cpus} {ref_path} {in_filepath} -o {out_filepath}"
    run_command(map_cmd)
    return out_filepath


def concat_fasta_2(in_filepaths: list, out_filepath):
    """Concatenate fasta sequences into single fasta file.
    Takes a list of fasta filepaths and an output filename for saving"""
    cat_cmd = f"echo {' '.join(in_filepaths)} | xargs cat >> {out_filepath}"
    # cat_cmd = f"cat {' '.join(in_filepaths)} > {out_filepath}"
    # print(cat_cmd)
    run_command(cat_cmd)
    return out_filepath


def sample_fasta(fasta_filepath, out_filepath, sample_size=100):
    "Sample the first n sequences from the input fasta file"
    sampled_seqs = []
    seqs = load_fasta(fasta_filepath)
    for i, rec in enumerate(seqs):
        if i < sample_size:
            sampled_seqs.append(rec)
        else:
            SeqIO.write(sampled_seqs, out_filepath, 'fasta');
            return out_filepath


def concat_fasta(in_dir, out_dir):
    """Concatenate fasta sequences into single fasta file"""
    cat_cmd = f"cat {in_dir}/*.fa* > {out_dir}.fa"
    subprocess.check_call(cat_cmd, shell=True)
    return f"{out_dir}.fa"


def run_command(cmd):
    "helper function used to run bash commands"
    # run it
    p = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
    # But do not wait till netstat finish, start displaying output immediately
    while True:
        out = p.stderr.read(1)
        if out == b'' and p.poll() != None:
            break
        if out != b'':
            sys.stdout.write(out.decode('utf-8'))
            sys.stdout.flush()


def run_command_log(cmd):
    "helper function used to run bash commands"
    # run it
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # But do not wait till netstat finish, start displaying output immediately
    out = p.stdout.read()
    return out.decode('utf-8')


def align_fasta(fasta_filepath, out_filepath, num_cpus=8):
    """Generate Multiple Sequence Alignment of concatenated sequences in input fasta file using mafft.
    TODO: ALLOW USER TO INPUT CUSTOM COMMAND"""
    msa_cmd = f"mafft --auto --thread {num_cpus} {fasta_filepath} > {out_filepath}"
    run_command(msa_cmd)
    return out_filepath


def align_fasta_reference(fasta_filepath, out_filepath, ref_fp: str, num_cpus=8):
    """Generate Multiple Sequence Alignment of concatenated sequences in input fasta file using mafft"""
    msa_cmd = f"mafft --auto --thread {num_cpus} --keeplength --addfragments {fasta_filepath} {ref_fp} > {out_filepath}"
    run_command(msa_cmd)
    return out_filepath


def align_fasta_viralMSA(fasta_filepath, out_filepath, ref_fp: str, 
                         num_cpus=8, email='notmyemail@noway.com',
                         viralmsa_fp = '/home/al/code/ViralMSA.py'):
    """Generate Multiple Sequence Alignment of concatenated sequences in input fasta file using mafft"""
    msa_cmd = f"{viralmsa_fp} -e {email} -t {num_cpus} -s {fasta_filepath} -o {out_filepath} -r {ref_fp}"
    run_command(msa_cmd)
    out_filepath = out_filepath + '/' + fasta_filepath.split('/')[-1] + '.aln'
    return out_filepath


def compute_tree(msa_filepath, num_cpus=8, redo=False):
    """Compute ML tree of aligned sequences in input fasta using iqtree"""
    out_filepath = msa_filepath + '.treefile'
    if redo:
        tree_cmd = f"iqtree -redo -s {msa_filepath} -nt {num_cpus} -m HKY -czb -fast"
    else:
        tree_cmd = f"iqtree -s {msa_filepath} -nt {num_cpus} -m HKY -czb -fast"
    run_command(tree_cmd)
    return out_filepath


def compute_time_tree(msa_filepath, tree_filepath, num_cpus=8):
    """Compute time tree (???)"""
    out_path = '/'.join(msa_filepath.split('/')[:-1]) + '/timetree'
    tree_cmd = f"treetime ancestral --aln {msa_filepath} --tree {tree_filepath} --outdir {out_path}"
    subprocess.check_call(tree_cmd, shell=True)
    return out_path


def add_gene_column(df: pd.DataFrame) -> pd.DataFrame:
    """Takes dataframe containing intra-host variant information and adds a column indicating the gene for each mutation
     Uses the GFF_FEATURE to infer gene, and when GFF_FEATURE is missing (i.e. for indels), it uses nucleotide position 
     to infer gene. The `gff2gene` mapping was obtained from the file /home/al/data/hcov19/gff/NC_045512.2.gff3 and using 
     the create_gff2gene_mapping() function defined below."""
    gff2gene = {
        'cds-YP_009724389.1': 'ORF1ab',
        'cds-YP_009725295.1': 'ORF1ab',
        'cds-YP_009724390.1': 'S',
        'cds-YP_009724396.1': 'ORF8',
        'cds-YP_009724397.2': 'N',
        'cds-YP_009724395.1': 'ORF7a',
        'cds-YP_009724391.1': 'ORF3a',
        'cds-YP_009724393.1': 'M',
        'cds-YP_009724394.1': 'ORF6',
        'cds-YP_009725255.1': 'ORF10',
        'cds-YP_009725318.1': 'ORF7b',
        'cds-YP_009724392.1': 'E'
    }
    # infer the gene from GFF_FEATURE
    try:
        df['gene'] = df['GFF_FEATURE'].apply(lambda x: gff2gene.get(x, 'nan'))
    except:
        raise KeyError('GFF_FEATURE column not found in the input dataframe.')
    # infer the gene from position when GFF_FEATURE is missing
    try:
        df.loc[df['gene']=='nan', 'gene'] = df.loc[df['gene']=='nan', 'POS'].apply(map_gene_to_pos)
    except:
        raise KeyError('POS column not found in the input dataframe.')
    return df
    

def map_gene_to_pos(x):
    """helper function to infer the gene based on nucleotide position of SARS-CoV-2 genome"""
    pos = x
    if pos <= 265:
        return '5UTR'
    elif pos > 265 and pos <= 13466:
        return 'ORF1a'
    elif pos > 13466 and pos <= 21555:
        return 'ORF1b'
    elif pos > 21562 and pos <= 25384:
        return 'S'
    elif pos > 25392 and pos <= 26220:
        return 'ORF3a'
    elif pos > 26244 and pos <= 26472:
        return 'E'
    elif pos > 26522 and pos <= 27191:
        return 'M'
    elif pos > 27201 and pos <= 27387:
        return 'ORF6'
    elif pos > 27393 and pos <= 27759:
        return 'ORF7a'
    elif pos > 27755 and pos <= 27887:
        return 'ORF7b'
    elif pos > 27893 and pos <= 28259:
        return 'ORF8'
    elif pos > 28273 and pos <= 29533:
        return 'N'
    elif pos > 29557 and pos <= 29674:
        return 'ORF10'
    elif pos > 29674:
        return '3UTR'
    return 'nan'


# def create_gff2gene_mapping(variant_data: pd.DataFrame, gff_filepath: str) -> dict:
#     fn = gffutils.example_filename(gff_filepath)
#     db = gffutils.create_db(fn, dbfn='test.db', merge_strategy='merge', force=True)
#     all_gffs = variant_data['GFF_FEATURE'].dropna().unique()
#     gff2gene = {}
#     for gff in all_gffs:
#         if gff:
#             gff2gene[gff] = db[gff]['gene'][0]
#     return gff2gene


def compute_acc_nt_pos(x, gene2pos):
    s = gene2pos.get(x['gene'], 0)
    return s + x['pos']


def compute_acc_aa_pos(x, gene2pos):
    s = gene2pos.get(x['gene'], 0)
    return s + x['codon_num']


def check_state(x, gaps=False):
    if x[-2:].isupper():
        if gaps:
            x = x[:-3]
        else:
            x = x[:-2]
    return x


def flatten_list(_2d_list):
    """flatten any 2D Python list
    Citation: https://stackabuse.com/python-how-to-flatten-list-of-lists/"""
    flat_list = []
    # Iterate through the outer list
    for element in _2d_list:
        if type(element) is list:
            # If the element is of type list, iterate through the sublist
            for item in element:
                flat_list.append(item)
        else:
            flat_list.append(element)
    return flat_list