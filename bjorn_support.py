import sys
import glob
import subprocess
import gzip
from path import Path
import numpy as np
import pandas as pd
from Bio import Seq, SeqIO, AlignIO, Phylo, Align


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
                  data_type: str='consensus', tech: str='illumina') -> dict:
    """Take list of sample IDs, the general area where your files are located, and their format.
       Returns a dictionary with sample IDs as keys and the filepaths as values.
    """
    file_paths = {}
    if sample_ids:
        for s_id in sample_ids:
            f = glob.glob(f"{analysis_path}/*{data_type}*/*{tech}*/*{s_id}*.{data_fmt}")
            try:
                file_paths[s_id] = f[0]
            except:
                continue
    else:
        fs = glob.glob(f"{analysis_path}/*{data_type}*/*{tech}*/*.{data_fmt}")
        for f in fs:
            try:
                sample_id = f.split('/')[-1].split('_')[0].split('-')[1]
            except:
                continue
            file_paths[sample_id] = f
    return file_paths


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
    cat_cmd = f"cat {' '.join(in_filepaths)} > {out_filepath}"
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


def create_gff2gene_mapping(variant_data: pd.DataFrame, gff_filepath: str) -> dict:
    fn = gffutils.example_filename(gff_filepath)
    db = gffutils.create_db(fn, dbfn='test.db', merge_strategy='merge', force=True)
    all_gffs = variant_data['GFF_FEATURE'].dropna().unique()
    gff2gene = {}
    for gff in all_gffs:
        if gff:
            gff2gene[gff] = db[gff]['gene'][0]
    return gff2gene


def compute_acc_nt_pos(x, gene2pos):
    s = gene2pos.get(x['gene'], 0)
    return s + x['pos']


def compute_acc_aa_pos(x, gene2pos):
    s = gene2pos.get(x['gene'], 0)
    return s + x['codon_num']