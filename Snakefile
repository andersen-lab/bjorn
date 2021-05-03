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


# load user parameters
configfile: "config.json"

username = config['username']
password = config['password']
out_dir = config['out_dir']
current_datetime = datetime.now().strftime("%Y-%m-%d-%H-%M")
gisaid_sequences_filepath = out_dir + '/' + config['gisaid_fasta'] + '_' + current_datetime + '.fasta'
meta_filepath = out_dir + '/' + config['gisaid_meta'] + '_' + current_datetime + '.tsv.gz'
info_filepath = out_dir + '/' + config['chunk_info']
chunks_dir = out_dir + '/chunks'
fasta_dir = chunks_dir + '/fasta/' + current_datetime
# sam_dir = chunks_dir + '/sam/' + current_date
# msa_dir = chunks_dir + '/msa/' + current_date
# muts_dir = chunks_dir + '/muts/' + current_date
logs_dir = out_dir + '/logs'
chunk_size = int(config['chunk_size'])
num_cpus = int(config['num_cpus'])
reference_filepath = config['ref_fasta']
patient_zero = config['patient_zero']

# Download and pre-process GISAID data
download_cmd = f"src/json2fasta.py -u {username} -p {password} -s {chunk_size} -t {current_datetime}"
bs.run_command(download_cmd)
info_df = pd.read_csv(info_filepath)
# info_df = bj.download_process_data(username, password, chunk_size)


rule all:
    input:
        "{out_dir}/mutations_{current_datetime}.csv".format(out_dir = out_dir, current_datetime = current_datetime), # output data (signal)
        expand("{chunks_dir}/muts/{current_datetime}/{sample}.mutations.csv", chunks_dir = chunks_dir, current_datetime = current_datetime, sample = info_df['chunk_names']), # bjorn
        expand("{chunks_dir}/msa/{current_datetime}/{sample}.aligned.fasta", chunks_dir = chunks_dir, current_datetime = current_datetime, sample = info_df['chunk_names']), # data2funk -> gofasta
        expand("{chunks_dir}/sam/{current_datetime}/{sample}.sam", chunks_dir = chunks_dir, current_datetime = current_datetime, sample = info_df['chunk_names']), # minimap2 -> mafft
        expand("{chunks_dir}/fasta/{current_datetime}/{sample}.fasta", chunks_dir = chunks_dir, current_datetime = current_datetime, sample = info_df['chunk_names']), # chunk_fasta
        gisaid_sequences_filepath, # input data (signal)
        info_filepath
        # reference_filepath, # input data (patient zero)




# TODO: create merge_mutations.py 
rule merge_results:
    input:
        expand("{chunks_dir}/muts/{current_datetime}/{sample}.mutations.csv", chunks_dir = chunks_dir, sample = info_df['chunk_names'], current_datetime = current_datetime),
        meta_filepath=meta_filepath
    threads: 1
    params:
        current_datetime=current_datetime,
    output:
        "{out_dir}/mutations_{current_datetime}.csv"
    shell:
        """
        src/merge_results.py -i {chunks_dir}/muts/{current_datetime}/ -m {input.meta_filepath} -o {output} -t {params.current_datetime}
        """


# TODO: test msa_2_mutations.py 
rule run_bjorn:
    input:
        "{chunks_dir}/msa/{current_datetime}/{sample}.aligned.fasta"
    params:
        patient_zero=patient_zero,
    output:
        "{chunks_dir}/muts/{current_datetime}/{sample}.mutations.csv"
    shell:
        """
        src/msa_2_mutations.py -i {input} -r {params.patient_zero} -o {output}
        """
        # for i, o in zip(input, output):
        #     _ = bm.msa_2_mutations(i, params.patient_zero, o, config)

rule run_data2funk:
    input:
        "{chunks_dir}/sam/{current_datetime}/{sample}.sam",
    params:
        reference_filepath=reference_filepath,
    output:
        "{chunks_dir}/msa/{current_datetime}/{sample}.aligned.fasta",
    shell:
        """
        datafunk sam_2_fasta -s {input} -r {params.reference_filepath} -o {output} --pad --log-inserts
        """

rule run_minimap2:
    input:
        "{chunks_dir}/fasta/{current_datetime}/{sample}.fasta"
    params:
        num_cpus=num_cpus,
        reference_filepath=reference_filepath
    output:
        "{chunks_dir}/sam/{current_datetime}/{sample}.sam",
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
        expand("{chunks_dir}/fasta/{current_datetime}/{sample}.fasta", chunks_dir = chunks_dir, current_datetime = current_datetime, sample = info_df['chunk_names'])
    shell:
        """
        src/chunk_fasta.py -f {input} -r {params.reference_filepath} -s {params.chunk_size} -o {chunks_dir}/fasta/{current_datetime}
        """
        # bf.chunk_fasta(gisaid_sequences_filepath, params.reference_filepath, params.chunk_size, Path(params.out_dir))


# rule load_info_df:
#     input:
#         info_filepath
#     run:
#         info_df=pd.read_csv(input)




# rule generate_data:
#     params:
#         username,
#         password
#     threads: 1
#     output:
#         gisaid_sequences_filepath,
#         meta_filepath,
#         info_filepath
#     shell:
#         """
#         src/json2fasta.py -u {params[0]} -p {params[1]}
#         """


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