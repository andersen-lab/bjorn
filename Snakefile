import sys
sys.path.append('src/')
from path import Path
from datetime import datetime
import pandas as pd
import json
import argparse
import src.bjorn_support as bs
import src.json2fasta as bj
import src.chunk_fasta as bf
import src.msa_2_mutations as bm


# load user parameters
configfile: "config.json"

username = config['username']
password = config['password']
out_dir = config['out_dir']
current_date = config['date']
gisaid_sequences_filepath = out_dir + '/' + config['gisaid_fasta']
meta_filepath = out_dir + '/' + config['gisaid_meta']
chunks_dir = out_dir + '/chunks'
fasta_dir = chunks_dir + '/fasta/' + current_date
# sam_dir = chunks_dir + '/sam/' + current_date
# msa_dir = chunks_dir + '/msa/' + current_date
# muts_dir = chunks_dir + '/muts/' + current_date
logs_dir = out_dir + '/logs'
chunk_size = int(config['chunk_size'])
num_cpus = int(config['num_cpus'])
reference_filepath = config['ref_fasta']
patient_zero = config['patient_zero']

# TODO: create dataframe containing chunked filenames
info_df = bj.download_process_data(username, password, chunk_size)


rule all:
    input:
        "{out_dir}/mutations_{current_date}.csv".format(out_dir = out_dir, current_date = current_date), # output data (signal)
        expand("{chunks_dir}/muts/{current_date}/{sample}.mutations.csv", chunks_dir = chunks_dir, current_date = current_date, sample = info_df['chunk_names']), # bjorn
        expand("{chunks_dir}/msa/{current_date}/{sample}.aligned.fasta", chunks_dir = chunks_dir, current_date = current_date, sample = info_df['chunk_names']), # data2funk -> gofasta
        expand("{chunks_dir}/sam/{current_date}/{sample}.sam", chunks_dir = chunks_dir, current_date = current_date, sample = info_df['chunk_names']), # minimap2 -> mafft
        expand("{chunks_dir}/fasta/{current_date}/{sample}.fasta", chunks_dir = chunks_dir, current_date = current_date, sample = info_df['chunk_names']), # chunk_fasta
        gisaid_sequences_filepath, # input data (signal)
        # reference_filepath, # input data (patient zero)


# TODO: create merge_mutations.py 
rule merge_results:
    input:
        expand("{chunks_dir}/muts/{current_date}/{sample}.mutations.csv", chunks_dir = chunks_dir, sample = info_df['chunk_names'], current_date = current_date),
        meta_filepath=meta_filepath
    threads: 1
    output:
        "{out_dir}/mutations_{current_date}.csv"
    shell:
        """
        python src/merge_results.py -i {chunks_dir}/muts/{current_date}/ -m {input.meta_filepath} -o {output}
        """


# TODO: test msa_2_mutations.py 
rule run_bjorn:
    input:
        "{chunks_dir}/msa/{current_date}/{sample}.aligned.fasta"
    params:
        patient_zero=patient_zero,
    output:
        "{chunks_dir}/muts/{current_date}/{sample}.mutations.csv"
    shell:
        """
        python src/msa_2_mutations.py -i {input} -r {params.patient_zero} -o {output}
        """
        # for i, o in zip(input, output):
        #     _ = bm.msa_2_mutations(i, params.patient_zero, o, config)

rule run_data2funk:
    input:
        "{chunks_dir}/sam/{current_date}/{sample}.sam",
    params:
        reference_filepath=reference_filepath,
    output:
        "{chunks_dir}/msa/{current_date}/{sample}.aligned.fasta",
    shell:
        """
        datafunk sam_2_fasta -s {input} -r {params.reference_filepath} -o {output} --pad --log-inserts
        """

rule run_minimap2:
    input:
        "{chunks_dir}/fasta/{current_date}/{sample}.fasta"
    params:
        num_cpus=num_cpus,
        reference_filepath=reference_filepath
    output:
        "{chunks_dir}/sam/{current_date}/{sample}.sam",
    shell:
        """
        minimap2 -a -x asm5 -t {params.num_cpus} {params.reference_filepath} {input} -o {output}
        """

rule chunk_fasta:
    input:
        gisaid_sequences_filepath,
    params:
        reference_filepath=reference_filepath,
        chunk_size=chunk_size,
        out_dir=fasta_dir
        # out_dir=lambda wildcards, output: Path(output).parent
    threads: 1
    output:
        expand("{chunks_dir}/fasta/{current_date}/{sample}.fasta", chunks_dir = chunks_dir, current_date = current_date, sample = info_df['chunk_names'])
    shell:
        """
        python src/chunk_fasta.py -f {input} -r {params.reference_filepath} -s {params.chunk_size} -o {chunks_dir}/fasta/{current_date}
        """
        # bf.chunk_fasta(gisaid_sequences_filepath, params.reference_filepath, params.chunk_size, Path(params.out_dir))

# rule compute_chunk_size:
#     input:
#         meta_filepath,
#     params:
#         chunk_size=chunk_size
#     threads: 1
#     run:
#         info_df = bs.create_chunk_names(input[0], params.chunk_size)


# rule generate_data:
#     params:
#         username,
#         password
#     threads: 1
#     output:
#         gisaid_sequences_filepath,
#         meta_filepath
#     shell:
#         """
#         python src/json2fasta.py -u {params[0]} -p {params[1]}
#         """