import pandas as pd
import matplotlib.pylab as plt
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from path import Path
from shutil import copy
from Bio import SeqIO, AlignIO, Phylo
import argparse
import glob
import subprocess
from multiprocessing import Pool
from itertools import repeat
import os
from datetime import datetime as dt
from bjorn_support import concat_fasta, align_fasta, compute_tree, map_gene_to_pos, load_fasta
from mutations import identify_replacements, identify_deletions, identify_insertions
from onion_trees import load_tree, visualize_tree, get_indel2color, get_sample2color
import data as bd


## FUNCTION DEFINTIONS
def create_sra_meta(df: pd.DataFrame, sra_dir: Path):
    biosample_paths = glob.glob(f"{sra_dir}/*.txt")
    if biosample_paths:
        biosample_df = pd.concat((pd.read_csv(f, sep='\t') for f in biosample_paths))
        biosample_df["sample_id"] = biosample_df["Sample Name"].apply(lambda x: "".join(x.split("-")[:2]))
        bam_files = ans.loc[~ans['PATH_y'].isna()][['sample_id', 'PATH_y']].rename(columns={'PATH_y': 'file_name'})
        sra_merged = pd.merge(biosample_df, bam_files, on="sample_id")
        sra_merged[["Accession", "Sample Name", "file_name"]].to_csv(sra_dir/"sra_metadata.csv", index=False)
    return f"SRA metadata saved in {sra_dir/'sra_metadata.csv'}"


def assemble_genbank_release(cns_seqs: list, df: pd.DataFrame, meta_cols: list, genbank_dir: Path):
    # create directory for genbank release
    if not Path.isdir(genbank_dir):
        Path.mkdir(genbank_dir);
    authors = {}
    # group samples by author
    for ctr, (n, grp) in enumerate(df.groupby('authors')):
        authors[ctr+1] = n
        # generate sample metadata
        genbank_meta = create_genbank_meta(grp, meta_cols)
        genbank_meta.to_csv(genbank_dir/f'genbank_metadata_{ctr+1}.tsv', sep='\t', index=False)
        # fetch consensus sequences of those samples
        recs = [i for i in cns_seqs if i.name in genbank_meta['Sequence_ID'].tolist()]
        SeqIO.write(recs, genbank_dir/f'genbank_release_{ctr+1}.fa', 'fasta')
    # write mapping of index to author for later reference
    (pd.DataFrame.from_dict(authors, orient='index')
       .rename(columns={0: 'authors'})
       .to_csv(genbank_dir/'authors.tsv', sep='\t'))
    return f"Genbank data release saved in {genbank_dir}"


def create_genbank_meta(df: pd.DataFrame, meta_cols: list) -> pd.DataFrame:
    genbank_meta = df[meta_cols].copy()
    genbank_meta['country'] = genbank_meta['location'].apply(lambda x: x.split('/')[0])
    genbank_meta['isolate'] = genbank_meta['Virus name'].str.replace('hCoV-19', 'SARS-CoV-2/human')
    genbank_meta['host'] = genbank_meta['Host'].str.replace('Human', 'Homo Sapiens')
    genbank_meta.rename(columns={'Virus name': 'Sequence_ID', 'collection_date': 'collection-date',
                                 'Specimen source': 'isolation-source'}, inplace=True)
    genbank_meta.loc[genbank_meta['country']=='MEX', 'country'] = 'Mexico'
    return genbank_meta[['Sequence_ID', 'isolate', 'country',
                         'collection-date', 'host', 'isolation-source']]


def create_github_meta(new_meta_df: pd.DataFrame, old_meta_filepath: str, meta_cols: list):
    """Generate Github metadata with updated information about newly released samples"""
    old_metadata = pd.read_csv(old_meta_filepath)
    new_metadata = pd.concat([old_metadata, new_meta_df.loc[:, meta_cols]])
    new_metadata.to_csv(out_dir/'metadata.csv', index=False)
    return f"Github metadata saved in {out_dir/'metadata.csv'}"

def create_gisaid_meta(new_meta_df: pd.DataFrame, meta_cols: list):
    """Generate GISAID metadata for newly released samples"""
    new_meta_df[meta_cols].to_csv(out_dir/'gisaid_metadata.csv', index=False)
    return f"GISAID metadata saved in {out_dir/'gisaid_metadata.csv'}"

def get_ids(filepaths: list) -> list:
    "Utility function to get a unified sample ID format from the analysis file paths"
    ids = []
    for fp in filepaths:
        # query = fp.basename().split('-')
        n = fp.basename()
        if n[:6] != "SEARCH":
            n = n[n.find("SEARCH"):]
        query = n.split('_')[0].split('-')
        if len(query) > 1:
            ids.append(''.join(query[:2]))
        else:
            start_idx = fp.find('SEARCH')
            ids.append(fp[start_idx:start_idx+10])
    return ids

def process_coverage_sample_ids(x):
    "Utility function to get a unified sample ID format from the coverage reports"
    if x[:6] != "SEARCH":
        x = x[x.find("SEARCH"):]
    query = x.split('/')
    if len(query) == 1:
        # query = fp.basename().split('_')[0].split('-')
        return ''.join(x.split("_")[0].split('-')[:2]) # Format could be SEARCH-xxxx-LOC or SEARCH-xxxx
    else:
        start_idx = x.find('SEARCH')
        return x[start_idx:start_idx+10] # SEARCHxxxx

def compress_files(filepaths: list, destination='/home/al/tmp2/fa/samples.tar.gz'):
    "Utility function to compress list of files into a single .tar.gz file"
    with tarfile.open(destination, "w:gz") as tar:
        for f in filepaths:
            tar.add(f)
    return 0

def transfer_files(filepaths: pd.DataFrame, destination: str, include_bams=False, ncpus = 1):
    "Utility function to copy consensus and BAM files of given samples from source to destination"
    filepaths = filepaths[['PATH_x', 'PATH_y', "Virus name"]]
    destination = Path(destination)
    if not Path.isdir(destination):
        Path.mkdir(destination)
    if not Path.isdir(destination/'fa/'):
        Path.mkdir(destination/'fa/')
    if not Path.isdir(destination/'bam/'):
        Path.mkdir(destination/'bam/')
    with Pool(ncpus) as pool:
        res = pool.starmap(transfer_cns_sequence,
                           zip(repeat(destination/"fa/"),
                               filepaths["PATH_x"].tolist(), filepaths["Virus name"].tolist()))
        if include_bams:
            res = pool.starmap(transfer_bam, zip(repeat(destination/"bam/"),
                                                 filepaths["PATH_y"].tolist()))
        pool.close()
        pool.join()
    return 0


def transfer_bam(out, bam_fp):
    n = os.path.join(out, os.path.basename(bam_fp))
    if Path.isfile(Path(n)): return 0
    print("Transferring {} to {}".format(bam_fp, n))
    out = subprocess.Popen(["samtools", "view", "-b", "-F", "4", "-o", n, bam_fp],
                           stdout = subprocess.PIPE,
                           stderr = subprocess.PIPE, shell=False)
    stdout, stderr = out.communicate()
    if len(stderr) !=0:
        print(stderr)
    return 0

def transfer_cns_sequence(out, cons_fp, virus_name):
    n = os.path.join(out, os.path.basename(cons_fp))
    if Path.isfile(Path(n)): return 0
    print("Transferring {} to {}".format(cons_fp, n))
    rec = SeqIO.read(cons_fp, format="fasta")
    new_cns_name = virus_name
    rec.name = new_cns_name
    rec.id = new_cns_name
    rec.description = new_cns_name
    with open(n, "w") as fout:
        SeqIO.write(rec, fout, format="fasta")
        fout.close()
    return 0

def process_id(x):
    "Utility function to process sample IDs to fix inconsistencies in the format"
    return ''.join(x.split('-')[:2])

## MAIN

if __name__=="__main__":
    # Input Parameters
    # COLUMNS TO INCLUDE IN GITHUB METADATA
    git_meta_cols = ["ID", "collection_date", "location", "percent_coverage_cds", "avg_depth", "authors", "originating_lab", "fasta_hdr"]
    # COLUMNS TO INCLUDE IN GISAID METADATA
    gisaid_meta_cols = ['Submitter',
                   'FASTA filename', 'Virus name', 'Type', 'Passage details/history',
                   'Collection date', 'location', 'Additional location information',
                   'Host', 'Additional host information', 'Gender', 'Patient age',
                   'Patient status', 'Specimen source', 'Outbreak', 'Last vaccinated',
                   'Treatment', 'Sequencing technology', 'Assembly method', 'Coverage',
                   'originating_lab', 'Address', 'Sample ID given by the sample provider',
                   'Submitting lab', 'Address.1',
                   'Sample ID given by the submitting laboratory', 'authors', 'avg_depth']
    # COLUMNS TO INCLUDE IN GITHUB METADATA
    genbank_meta_cols = ['Sample ID', 'ID', 'Virus name', 'location', 'Specimen source', 'collection_date', 'Host']
    parser = argparse.ArgumentParser()

    parser.add_argument('-d', "--not-dry-run", action='store_false', help="Dry run. Default: True")

    parser.add_argument("-b", "--include-bams", action='store_true', help="Whether or not to include BAM files in the release")

    parser.add_argument("-r", "--reference",
                        type=str,
                        default="/home/gk/code/hCoV19/db/NC045512.fasta",
                        help="Reference to use")

    parser.add_argument("-rn", "--reference-name",
                        type=str,
                        default="NC_045512.2",
                        help="FASTA header name for reference")

    parser.add_argument("-o", "--out-dir",
                        type=str,
                        default="./",
                        help="Output directory")

    parser.add_argument("-n", "--cpus",
                        type=int,
                        default=1,
                        help="Number of cpus to use")

    parser.add_argument("-s", "--sample-sheet",
                        type=str,
                        required=True,
                        help="Sample sheet to use. Please make sure it is the most recent")

    parser.add_argument("-a", "--analysis-folder",
                        type=str,
                        required=True,
                        help="Path to the analysis folder that contains all the files")

    parser.add_argument("-m", "--output-metadata",
                        type=str,
                        required=True,
                        help="Output csv for metadata of the released samples")
    args = parser.parse_args()

    # whether or not to include bam files in the release
    include_bams = args.include_bams
    # these are the columns to include in the metadata.csv output

    # path to reference sequence (used later for MSA and tree construction)
    ref_path = Path(args.reference) if args.reference is not None else None
    patient_zero = args.reference_name
    # this is the directory where results get saved
    out_dir = Path(args.out_dir)
    # number of cores to use
    num_cpus = args.cpus
    # file path to samples sheet (make sure it is the most recent)
    sample_sheet_fpath = args.sample_sheet
    # path to analysis results
    analysis_fpath = args.analysis_folder
    # file path to metadata of samples that have already been released
    released_samples_fpath = args.output_metadata
    # Whether run is dry
    dry_run = args.not_dry_run

    # # Test
    # out_dir = "/home/gk/southpark/2020-11-21_release"
    # sample_sheet_fpath = "/home/gk/code/hCoV19/release_summary_csv/2020-11-20_seq_summary.csv"
    # analysis_fpath = "/home/gk/analysis/"
    # released_samples_fpath = "/home/gk/analysis/hcov-19-genomics/metadata.csv"
    # dry_run = True

    print(f"""User Specified Parameters:
    Dry run: {dry_run}.
    Include BAMS: {include_bams}.
    Reading release summary file from {sample_sheet_fpath}.
    Reading repository metadata from {released_samples_fpath}.
    Searching analysis folder {analysis_fpath}.
    """
    )

    # Collecting Sequence Data

    # grab all filepaths for bam data
    bam_filepaths = glob.glob(f"{analysis_fpath}/**/merged_aligned_bams/illumina/*.bam")
    bam_filepaths = [Path(fp) for fp in bam_filepaths]
    # consolidate sample ID format
    bam_ids = get_ids(bam_filepaths)
    # Turn into dataframe
    bam_data = list(zip(*[bam_ids, bam_filepaths]))
    bam_df = pd.DataFrame(data=bam_data, columns=['sample_id', 'PATH'])
    # grab all paths to consensus sequences
    consensus_filepaths = glob.glob(f"{analysis_fpath}/**/consensus_sequences/illumina/*.fa")
    consensus_filepaths = [Path(fp) for fp in consensus_filepaths]
    # consolidate sample ID format
    consensus_ids = get_ids(consensus_filepaths)
    # Turn into dataframe
    consensus_data = list(zip(*[consensus_ids, consensus_filepaths]))
    consensus_df = pd.DataFrame(data=consensus_data, columns=['sample_id', 'PATH'])
    # clean up cns and bam (remove duplicate IDs)
    bam_df.drop_duplicates(subset=['sample_id'], keep='last', inplace=True)
    consensus_df.drop_duplicates(subset=['sample_id'], keep='last', inplace=True)
    # include only SEARCH samples
    consensus_df = consensus_df[(consensus_df['sample_id'].str.contains('SEARCH'))]
    # merge consensus and bam filepaths for each sample ID
    analysis_df = pd.merge(consensus_df, bam_df, on='sample_id', how='left')
    # load sample sheet data (GISAID) - make sure to download most recent one
    seqsum = pd.read_csv(sample_sheet_fpath)
    # clean up
    seqsum = seqsum[(~seqsum['SEARCH SampleID'].isna()) & (seqsum['SEARCH SampleID']!='#REF!')]
    # consolidate sample ID format
    seqsum.loc[:, 'sample_id'] = seqsum['SEARCH SampleID'].apply(process_id)
    seqsum.drop_duplicates(subset=['sample_id'], keep='last', inplace=True)
    seqsum = seqsum[seqsum['New sequences ready for release'] == 'Yes']
    num_seqs_to_release = seqsum['sample_id'].unique().shape[0]
    # JOIN summary sheet with analysis meta data
    sequence_results = pd.merge(seqsum, analysis_df, on='sample_id', how='inner')
    # compute number of samples with missing consensus and/or bam files
    num_seqs_found = sequence_results['sample_id'].unique().shape[0]
    num_samples_missing_cons = num_seqs_to_release - num_seqs_found
    num_samples_missing_bams = 'NA'
    if include_bams:
        # exclude any samples that do not have BAM data
        num_samples_missing_bams = sequence_results[sequence_results['PATH_y'].isna()].shape[0]
        sequence_results = sequence_results[~sequence_results['PATH_y'].isna()]
    # samples missing consensus or BAM sequence files
    num_samples_missing_bams = sequence_results[sequence_results['PATH_y'].isna()].shape[0]
    num_samples_missing_cons = sequence_results[sequence_results['PATH_x'].isna()].shape[0]
    # ## Make sure to remove any samples that have already been uploaded to github (just an extra safety step)
    # load metadata.csv from github repo, then clean up
    meta_df = pd.read_csv(released_samples_fpath)
    meta_df = meta_df[meta_df['ID'].str.contains('SEARCH')]
    # consolidate sample ID format
    meta_df.loc[:, 'sample_id'] = meta_df['ID'].apply(process_id)
    # meta_df['sample_id']
    # get IDs of samples that have already been released
    released_seqs = meta_df['sample_id'].unique()
    # filter out released samples from all the samples we got
    final_result = sequence_results[~sequence_results['sample_id'].isin(released_seqs)]
    print(f"Preparing {final_result.shape[0]} samples for release")
    # ## Getting coverage information
    cov_filepaths = glob.glob("{}/**/trimmed_bams/illumina/reports/*.tsv".format(analysis_fpath))
    # get_ipython().getoutput("find {analysis_fpath} -type f -path '*trimmed_bams/illumina/reports*' -name '*.tsv'")
    cov_filepaths = [Path(fp) for fp in cov_filepaths]
    # read coverage data and clean it up
    cov_df = pd.concat((pd.read_csv(f, sep='\t').assign(path=f) for f in cov_filepaths))
    cov_df.loc[:,'sample_id'] = cov_df['SAMPLE'].apply(process_coverage_sample_ids)
    cov_df.loc[:,'date'] = cov_df['path'].apply(lambda x: ''.join(x.split('/')[4].split('.')[:3]))
    cov_df = (cov_df.sort_values('date')
              .drop_duplicates(subset=['sample_id'], keep='last'))
    # JOIN results with coverage info
    ans = (
    pd.merge(final_result, cov_df,
             on='sample_id', how='left')
      .assign(
        collection_date = lambda x: pd.to_datetime(x["Collection date"]).dt.strftime("%Y-%m-%d")
    )
      .rename(columns={
        "SEARCH SampleID": "ID",
        "Location": "location",
        "COVERAGE": "percent_coverage_cds",
        "AVG_DEPTH": "avg_depth",
        "Authors": "authors",
        "Originating lab": "originating_lab"
    })
    )
    ans['fasta_hdr'] = ans['Virus name']
    num_samples_missing_coverage = ans[ans['percent_coverage_cds'].isna()].shape[0]
    # compute number of samples below 90% coverage
    low_coverage_samples = ans[ans["percent_coverage_cds"] < 90]
    # ignore samples below 90% coverage
    ans = ans[ans["percent_coverage_cds"] >= 90]
    # generate concatenated consensus sequences
    if not dry_run:
        # Transfer files
        transfer_files(ans, out_dir, include_bams=include_bams, ncpus=num_cpus)
        msa_dir = out_dir/'msa'
        if not Path.isdir(msa_dir):
            Path.mkdir(msa_dir);
        seqs_dir = Path(out_dir/'fa')
        copy(ref_path, seqs_dir);
        seqs_fp = concat_fasta(seqs_dir, msa_dir/out_dir.basename());
        # load concatenated sequences
        cns_seqs = SeqIO.parse(msa_dir/out_dir.basename()+'.fa', 'fasta')
        cns_seqs = list(cns_seqs)
        # generate files containing metadata for Github, GISAID, GenBank
        create_github_meta(ans, released_samples_fpath, git_meta_cols)
        create_gisaid_meta(ans, gisaid_meta_cols)
        assemble_genbank_release(cns_seqs, ans, genbank_meta_cols, out_dir/'genbank')
        sra_dir = out_dir/'sra'
        if not Path.isdir(sra_dir):
            Path.mkdir(sra_dir);
        input(f"\n Have you received the BioSample.txt files and placed them inside {sra_dir}? \n Press Enter to continue...")
        create_sra_meta(ans, sra_dir)
        # generate file containing deletions found
        # generate multiple sequence alignment
        msa_fp = seqs_fp.split('.')[0] + '_aligned.fa'
        if not Path.isfile(Path(msa_fp)):
            msa_fp = align_fasta(seqs_fp, msa_fp, num_cpus=num_cpus);
        # compute ML tree
        tree_dir = out_dir/'trees'
        if not Path.isdir(tree_dir):
            Path.mkdir(tree_dir);
        tree_fp = msa_fp + '.treefile'
        if not Path.isfile(Path(tree_fp)):
            tree_fp = compute_tree(msa_fp, num_cpus=num_cpus)
        tree = load_tree(tree_fp, patient_zero)
        # Plot and save basic tree
        fig1 = visualize_tree(tree)
        fig1.savefig(tree_dir/'basic_tree.pdf')
        # PLOT AND SAVE INDEL TREES
        colors = list(mcolors.TABLEAU_COLORS.keys())
        # path to new github metadata
        meta_fp = out_dir/'metadata.csv'
        # load multiple sequence alignment
        msa_data = load_fasta(msa_fp, is_aligned=True)
        # identify insertions
        insertions = identify_insertions(msa_data, 
                                         meta_fp=meta_fp, 
                                         patient_zero=patient_zero, 
                                         min_ins_len=1,
                                         data_src='alab')
        # save insertion results to file
        insertions.to_csv(out_dir/'insertions.csv', index=False)
        # identify substitution mutations
        subs = identify_replacements(msa_data,
                                    meta_fp=meta_fp,
                                    data_src='alab',
                                    patient_zero=patient_zero)
        # save substitution results to file
        subs.to_csv(out_dir/'replacements.csv', index=False)
        # identify deletions
        deletions = identify_deletions(msa_data,
                                        meta_fp=meta_fp,
                                        data_src='alab',
                                        patient_zero=patient_zero,
                                        min_del_len=1)
        # save deletion results to file
        deletions.to_csv(out_dir/'deletions.csv', index=False)
        # plot Phylogenetic tree with top consensus deletions annotated
        deletions = deletions.nlargest(len(colors), 'num_samples')
        # del2color = get_indel2color(deletions, colors)
        # sample_colors = get_sample2color(deletions, colors)
        # fig2 = visualize_tree(tree, sample_colors,
        #            indels=deletions, colors=colors);
        # fig2.savefig(tree_dir/'deletion_cns_tree.pdf', dpi=300)
        # fig3 = visualize_tree(tree, sample_colors,
        #                   indels=deletions, colors=colors,
        #                   isnv_info=True);
        # fig3.savefig(tree_dir/'deletion_isnv_tree.pdf', dpi=300)
        # plot Phylogenetic tree with top consensus deletions annotated
        insertions = insertions.nlargest(len(colors), 'num_samples')
        # del2color = get_indel2color(insertions, colors)
        # sample_colors = get_sample2color(insertions, colors)
        # fig4 = visualize_tree(tree, sample_colors,
        #            indels=insertions, colors=colors);
        # fig4.savefig(tree_dir/'insertion_cns_tree.pdf', dpi=300)
        # fig5 = visualize_tree(tree, sample_colors,
        #                   indels=insertions, colors=colors,
        #                   isnv_info=True);
        # fig5.savefig(tree_dir/'insertion_isnv_tree.pdf', dpi=300)
    if not Path.isdir(out_dir):
            Path.mkdir(out_dir);
    # Data logging
    with open("{}/data_release.log".format(out_dir), 'w') as f:
        f.write(f"Prepared {final_result.shape[0]} samples for release\n")
        f.write(f'{num_samples_missing_coverage} samples are missing coverage information\n')
        f.write(f'{low_coverage_samples.shape[0]} samples were found to have coverage below 90%\n')
        f.write(f'{num_samples_missing_cons} samples were ignored because they were missing consensus sequence files\n')
        f.write(f'{num_samples_missing_bams} samples were ignored because they were missing BAM sequence files\n')
    print(f"Transfer Complete. All results saved in {out_dir}")
