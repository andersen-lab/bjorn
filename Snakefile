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
config_file = "config.json"

configfile: config_file

if "username" and "password" in config:
    username = config['username']
    password = config['password']
out_dir = config['out_dir']
current_datetime = config["current_datetime"] if config["current_datetime"] != False else datetime.now().strftime("%Y-%m-%d-%H-%M")
gisaid_sequences_filepath = out_dir + '/' + config['gisaid_fasta'] + '_' + current_datetime + '.fasta'
meta_filepath = out_dir + '/' + config['gisaid_meta'] + '_' + current_datetime + '.tsv.gz'
info_filepath = out_dir + '/' + config['chunk_info']
logs_dir = out_dir + '/logs'
chunk_size = int(config['chunk_size'])
num_cpus = int(config['num_cpus'])
reference_filepath = config['ref_fasta']
patient_zero = config['patient_zero']
gadm_data = config["gadm_data"]

rule all:
    input:
        "{out_dir}/api_data_{current_datetime}.json.gz".format(current_datetime = current_datetime, out_dir = out_dir)

rule clean_chunks:
    params:
        out_dir = out_dir
    shell:
        """
        rm -r {params.out_dir}/chunks_*
        """

rule merge_json:
    input:
        dynamic("{out_dir}/chunks_apijson_{current_datetime}/{{sample}}.json.gz".format(out_dir = out_dir, current_datetime = current_datetime))
    output:
        "{out_dir}/api_data_{current_datetime}.json.gz"
    shell:
        """
        cat {input} > {output}
        """

rule merge_mutations_metadata:
    input:
        muts="{out_dir}/chunks_muts_{current_datetime}/{sample}.mutations.csv",
        metadata="{out_dir}/chunks_fasta_{current_datetime}/{sample}.tsv.gz"
    threads: 1
    params:
        current_datetime=current_datetime,
        geojson_prefix=config["geojson_prefix"],
        min_date = config["min_date"]
    output:
        temp("{out_dir}/chunks_apijson_{current_datetime}/{sample}.json.gz")
    shell:
        """
        src/merge_results.py -i {input.muts} -m {input.metadata} -o {output} -u None -n {params.min_date}  -g {params.geojson_prefix} -t {params.current_datetime}
        echo "" | gzip - >> {output} # Add new line as delimiter between chunks
        """

# TODO: test msa_2_mutations.py
rule run_bjorn:
    input:
        "{out_dir}/chunks_msa_{current_datetime}/{sample}.aligned.fasta"
    params:
        patient_zero=patient_zero
    output:
        temp("{out_dir}/chunks_muts_{current_datetime}/{sample}.mutations.csv")
    shell:
        """
        src/msa_2_mutations.py -i {input} -r {params.patient_zero} -o {output}
        """

rule run_data2funk:
    input:
        "{out_dir}/chunks_sam_{current_datetime}/{sample}.sam"
    params:
        reference_filepath=reference_filepath,
    output:
        temp("{out_dir}/chunks_msa_{current_datetime}/{sample}.aligned.fasta")
    shell:
        """
        datafunk sam_2_fasta -s {input} -r {params.reference_filepath} -o {output} --pad --log-inserts
        """

rule run_minimap2:
    input:
        "{out_dir}/chunks_fasta_{current_datetime}/{sample}.fasta"
    params:
        num_cpus=num_cpus,
        reference_filepath=reference_filepath
    output:
        temp("{out_dir}/chunks_sam_{current_datetime}/{sample}.sam")
    shell:
        """
        minimap2 -a -x asm5 -t {params.num_cpus} {params.reference_filepath} {input} -o {output}
        """

rule convert_to_fasta:
    input:
        "{out_dir}/chunks_json_{current_datetime}/chunk_json_{{sample}}".format(out_dir = out_dir, current_datetime = current_datetime)
    output:
        fasta=temp("{out_dir}/chunks_fasta_{current_datetime}/{sample}.fasta"),
        metadata=temp("{out_dir}/chunks_fasta_{current_datetime}/{sample}.tsv.gz")
    params:
        tmp="{out_dir}/chunks_fasta_{current_datetime}/{sample}.fasta.tmp",
        reference_filepath=reference_filepath,
        gadm_data = gadm_data,
        output_prefix="{out_dir}/chunks_fasta_{current_datetime}/{{sample}}".format(out_dir = out_dir, current_datetime = current_datetime)
    shell:
        """
        src/json2fasta.py -i {input} -o {params.output_prefix} -g {params.gadm_data}
        cat {output.fasta} {params.reference_filepath} > {params.tmp}
        mv {params.tmp} {output.fasta}
        """

rule chunk_json:
    input:
        "{out_dir}/provision_{current_datetime}.json".format(gisaid_feed = config['gisaid_feed'], out_dir = out_dir, current_datetime = current_datetime)
    params:
        chunk_dir = "{out_dir}/chunks_json_{current_datetime}".format(current_datetime = current_datetime, out_dir = out_dir),
        chunk_prefix = "chunk_json_{current_datetime}.".format(current_datetime = current_datetime),
        chunk_size = chunk_size
    threads: 1
    output:
        temp(dynamic("{out_dir}/chunks_json_{current_datetime}/chunk_json_{{sample}}".format(out_dir = out_dir, current_datetime = current_datetime)))
    shell:
        """
        split --verbose -l{params.chunk_size} -d -a 5 {input} {params.chunk_dir}/{params.chunk_prefix}
        """

rule download_sequences:
    output:
        temp("{out_dir}/provision_{current_datetime}.json")
    params:
        current_datetime = current_datetime
    shell:
        """
        curl -u {username}:{password} https://www.epicov.org/epi3/3p/scripps/export/provision.json.xz | xz -d -T8 > {output}
        """
