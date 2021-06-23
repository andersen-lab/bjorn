import sys
sys.path.append('src/')
from path import Path
from datetime import datetime
import pandas as pd
import json
import argparse
import bjorn_support as bs
import json2fasta as bj
import chunk_fasta as bf
import msa_2_mutations as bm
import os

# load user parameters
config_file = "config_test.json"

configfile: config_file

username = config['username']
password = config['password']
out_dir = config['out_dir']
# current_datetime = datetime.now().strftime("%Y-%m-%d-%H-%M") # Set to provision JSON while testing
current_datetime = "2021-06-22-17-25"
gisaid_sequences_filepath = out_dir + '/' + config['gisaid_fasta'] + '_' + current_datetime + '.fasta'
meta_filepath = out_dir + '/' + config['gisaid_meta'] + '_' + current_datetime + '.tsv.gz'
info_filepath = out_dir + '/' + config['chunk_info']
# sam_dir = chunks_dir + '/sam/' + current_date
# msa_dir = chunks_dir + '/msa/' + current_date
# muts_dir = chunks_dir + '/muts/' + current_date
logs_dir = out_dir + '/logs'
chunk_size = int(config['chunk_size'])
num_cpus = int(config['num_cpus'])
reference_filepath = config['ref_fasta']
patient_zero = config['patient_zero']

# Download and pre-process GISAID data
# download_cmd = f"src/json2fasta.py -u {username} -p {password} -s {chunk_size} -t {current_datetime} -c {config_file}"

# bs.run_command(download_cmd)
# info_df = pd.read_csv(info_filepath)
# info_df = bj.download_process_data(username, password, chunk_size)


rule all:
    input:
        "{out_dir}/mutations_{current_datetime}.csv".format(out_dir = out_dir, current_datetime = current_datetime), # output data (signal)

# TODO: create merge_mutations.py
rule merge_results:
    input:
        dynamic("{out_dir}/chunks_muts_{current_datetime}/{sample}.mutations.csv"),
        meta_filepath=meta_filepath
    threads: 1
    params:
        current_datetime=current_datetime,
        configfile = config_file,
        chunk_muts="{out_dir}/chunks_muts_{current_datetime}".format(out_dir = out_dir, current_datetime = current_datetime)
    output:
        "{out_dir}/mutations_{current_datetime}.csv"
    shell:
        """
        src/merge_results.py -i {params.chunk_muts}/ -m {input.meta_filepath} -o {output} -t {params.current_datetime} -c {params.configfile}
        """

# TODO: test msa_2_mutations.py
rule run_bjorn:
    input:
        "{out_dir}/chunks_msa_{current_datetime}/{sample}.aligned.fasta"
    params:
        patient_zero=patient_zero,
    output:
        "{out_dir}/chunks_muts_{current_datetime}/{sample}.mutations.csv"
    shell:
        """
        src/msa_2_mutations.py -i {input} -r {params.patient_zero} -o {output}
        """
        # for i, o in zip(input, output):
        #     _ = bm.msa_2_mutations(i, params.patient_zero, o, config)

rule run_data2funk:
    input:
        "{out_dir}/chunks_sam_{current_datetime}/{sample}.sam"
    params:
        reference_filepath=reference_filepath,
    output:
        "{out_dir}/chunks_msa_{current_datetime}/{sample}.aligned.fasta"
    shell:
        """
        datafunk sam_2_fasta -s {input} -r {params.reference_filepath} -o {output} --pad --log-inserts
        """

rule run_minimap2:
    input:
        "{out_dir}/chunks_fasta_{current_datetime}/fasta/{sample}.fasta"
    params:
        num_cpus=num_cpus,
        reference_filepath=reference_filepath
    output:
        "{out_dir}/chunks_sam_{current_datetime}/{sample}.sam",
    shell:
        """
        minimap2 -a -x asm5 -t {params.num_cpus} {params.reference_filepath} {input} -o {output}
        """

rule convert_to_fasta:
    input:
        "{out_dir}/chunks_json_{current_datetime}/{{sample}}".format(out_dir = out_dir, current_datetime = current_datetime)
    output:
        "{out_dir}/chunks_fasta_{current_datetime}/fasta/{sample}.fasta"
    shell:
        """
        json2fasta.py 
        src/chunk_fasta.py -f {input} -r {params.reference_filepath} -s {params.chunk_size} -o {chunks_dir}/fasta/{current_datetime}
        """

rule chunk_json:
    input:
        "{out_dir}/{gisaid_feed}_{current_datetime}.json".format(gisaid_feed = config['gisaid_feed'], out_dir = out_dir, current_datetime = current_datetime)
    params:
        chunk_dir = "{out_dir}/chunks_json_{current_datetime}/".format(current_datetime = current_datetime, out_dir = out_dir)
    threads: 1
    output:
        "{out_dir}/chunks_json_{current_datetime}/{{sample}}".format(out_dir = out_dir, current_datetime = current_datetime)
    shell:
        """
        split -l10 -d 5 {input} {params.chunk_dir}/{input}.
        """

rule download_sequences:
    output:
        "{out_dir}/{gisaid_feed}_{current_datetime}.json"
    params:
        current_datetime = current_datetime
    shell:
        """
        curl -u {username}:{password} ***REMOVED*** | xz -d -T8 > {output}
        """
