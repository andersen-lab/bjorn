from datetime import datetime

username = config['gisaid_username']
password = config['gisaid_password']
out_dir = config['out_dir']
include_hash = config['include_hash']
current_datetime = config["current_datetime"] if config["current_datetime"] != False else datetime.now().strftime("%Y-%m-%d-%H-%M")
chunk_size = int(config['chunk_size'])
num_cpus = int(config['num_cpus'])
reference_filepath = config['ref_fasta']
patient_zero = config['patient_zero']
data_source = config['data_source']
gadm_data = config["gadm_data"]
cleanup = config["cleanup"]
geojson_prefix = config["geojson_prefix"]
min_date = config["min_date"]
test_data = config["test_data"]
is_test = config["is_test"]

rule all:
    input:
        meta=f"{out_dir}/api_metadata_latest.json",
        data=f"{out_dir}/_api_data_{current_datetime}.json.gz"
    shell:
        """
        rm -f {out_dir}/api_metadata_latest.json
        ln -f {input.data} {out_dir}/api_data_latest.json.gz
        """

rule clear:
    shell:
        """
        rm -f {out_dir}/api_data_latest.json.gz
        rm -f {out_dir}/api_metadata_latest.json
        """

rule clean:
    shell:
        """
        rm -rf {out_dir}/chunks_*
        rm -f {out_dir}/*api*.json*
        """

rule build_meta:
    input:
        f"{out_dir}/_api_data_{current_datetime}.json.gz"
    output:
        "{out_dir}/api_metadata_latest.json"
    threads: 1
    shell:
        """
        echo "{{ \\""date_modified\\"": \\""{current_datetime}\\"", \
                 \\""records\\"": $(gunzip -c {input} | wc -l) " \
$( {include_hash} && echo \
               ",\\""hash\\"": \\""$(gunzip -c {input} | jq -cs 'sort|.[]' | md5sum | cut  -d' ' -f1)\\"" " \
) \
             "}}" | jq -c '.' > {out_dir}/_api_metadata_{current_datetime}.json
        cat {out_dir}/_api_metadata_{current_datetime}.json
        ln -f {out_dir}/_api_metadata_{current_datetime}.json {output}
        """

rule merge_json:
    input:
        dynamic(f"{out_dir}/chunks_apijson_{current_datetime}/{{sample}}.json.gz")
    output:
        "{out_dir}/_api_data_{current_datetime}.json.gz"
    threads: 1
    shell:
        """
        cat {input} > {output}
        """

rule merge_mutations_metadata:
    input:
        meta="{out_dir}/chunks_fasta_{current_datetime}/{sample}.tsv.gz",
        data="{out_dir}/chunks_muts_{current_datetime}/{sample}.mutations.csv"
    output:
        temp("{out_dir}/chunks_apijson_{current_datetime}/{sample}.json.gz")
    threads: 1
    shell:
        """
        python/merge_results.py -i {input.data} -m {input.meta} -o {output} -u None -n {min_date}  -g {geojson_prefix} -t {current_datetime}
        echo "" | gzip - >> {output} # Add new line as delimiter between chunks
        """

rule run_bjorn:
    input:
        "{out_dir}/chunks_msa_{current_datetime}/{sample}.aligned.fasta"
    output:
        temp("{out_dir}/chunks_muts_{current_datetime}/{sample}.mutations.csv")
    threads: 1
    shell:
        """
        python/msa_2_mutations.py -i {input} -r {patient_zero} -d {data_source} -o {output}
        """

rule run_data2funk:
    input:
        "{out_dir}/chunks_sam_{current_datetime}/{sample}.sam"
    output:
        temp("{out_dir}/chunks_msa_{current_datetime}/{sample}.aligned.fasta")
    threads: 1
    shell:
        """
        datafunk sam_2_fasta -s {input} -r {reference_filepath} -o {output} --pad --log-inserts
        """

rule run_minimap2:
    input:
        "{out_dir}/chunks_fasta_{current_datetime}/{sample}.fasta"
    output:
        temp("{out_dir}/chunks_sam_{current_datetime}/{sample}.sam")
    threads: num_cpus
    shell:
        """
        minimap2 -a -x asm5 -t {num_cpus} {reference_filepath} {input} -o {output}
        """

fasta_output_prefix = out_dir + "/chunks_fasta_" + current_datetime 
if data_source == "gisaid_feed":
    rule pull_gisaid_sequences:
        output:
            meta=temp(dynamic(f"{fasta_output_prefix}/{{sample}}.tsv.gz")),
            data=temp(dynamic(f"{fasta_output_prefix}/{{sample}}.fasta"))
        threads: 2*num_cpus
        shell:
            """
            mkdir -p /dev/shm/bjorn
            (({is_test} && gunzip -c {test_data}) || (curl -u {username}:{password} ***REMOVED*** | xz -d -T{num_cpus})) |
                    parallel --pipe --tmpdir /dev/shm/bjorn --block 30M -j{num_cpus} \
                        'python/json2fasta.py -o {fasta_output_prefix}/{{#}} -g {gadm_data} -r {reference_filepath}'
            """
elif data_source == "alab_release":
    rule clone_alab_sequences:
        output:
            meta=temp(dynamic(f"{fasta_output_prefix}/{{sample}}.tsv.gz")),
            data=temp(dynamic(f"{fasta_output_prefix}/{{sample}}.fasta"))
        threads: 1
        shell:
            """
            # git clone https://github.com/andersen-lab/HCoV-19-Genomics.git
            # copy to output, filter?
            # append reference sequence
            # line metadata up with the format expected by merge_results.py
            """
else:
    print(f'Error: data_source should be "gisaid_feed" or "alab_release" -- got {data_source}')
    sys.exit()
