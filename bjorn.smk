from datetime import datetime

username = config['gisaid_username']
password = config['gisaid_password']
out_dir = config['out_dir']
include_hash = config['include_hash']
current_datetime = config["current_datetime"] if config["current_datetime"] != False else datetime.now().strftime("%Y-%m-%d-%H-%M")
chunk_size = config['chunk_size']
min_length = int(config['min_length'])
max_unknown_pct = float(config['max_unknown_pct'])
max_task_cpus = int(config['max_task_cpus'])
max_cpus = int(config['max_cpus'])
rem_cpus = max_cpus - max_task_cpus
reference_filepath = "./pangolin/pangolin/data/reference.fasta"
patient_zero = "outgroup_A"
data_source = config['data_source']
gadm_data = config["gadm_data"]
geojson_prefix = config["geojson_prefix"]
min_date = config["min_date"]
test_data = config["test_data"]
is_test = config["is_test"]

rule all:
    input:
        meta=rules.build_meta.output,
        data=rules.merge_json.output
    output:
        meta=config["output_fp"],
        data=config["output_meta"]
    shell:
        """
        ln -f {input.meta} {output.meta}
        ln -f {input.data} {output.data}
        """

rule clear:
    shell:
        """
        rm -f rules.all.output.meta
        rm -f rules.all.output.data
        """

rule build_meta:
    input:
        rules.merge_json.output
    output:
        f"{out_dir}/_api_metadata_{current_datetime}.json"
    threads: max_task_cpus
    shell:
        """
        echo "[{current_datetime}, $(pigz -j{max_task_cpus} -dc {input} | wc -l) ]" |
            jq '{date_modified: .[0], records: .[1]}' > {output}
        cat {output}
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

module pangolin:
    snakefile: "pangolin/pangolin/scripts/pangolearn.smk" 

use rule align_to_reference from pangolin as pangalign
    input:
        fasta = "{out_dir}/chunks_fasta_{current_datetime}/{sample}.fasta"
    output:
        fasta = temp("{out_dir}/chunks_msa_{current_datetime}/{sample}.aligned.fasta")

#rule run_data2funk:
#    input:
#        "{out_dir}/chunks_sam_{current_datetime}/{sample}.sam"
#    output:
#        temp("{out_dir}/chunks_msa_{current_datetime}/{sample}.aligned.fasta")
#    threads: 1
#    shell:
#        """
#        datafunk sam_2_fasta -s {input} -r {reference_filepath} -o {output} --pad --log-inserts
#        """
#
#rule run_minimap2:
#    input:
#        "{out_dir}/chunks_fasta_{current_datetime}/{sample}.fasta"
#    output:
#        temp("{out_dir}/chunks_sam_{current_datetime}/{sample}.sam")
#    threads: max_task_cpus
#    shell:
#        """
#        minimap2 -a -x asm5 -t {max_task_cpus} {reference_filepath} {input} -o {output}
#        """

fasta_output_prefix = out_dir + "/chunks_fasta_" + current_datetime 
if data_source == "gisaid_feed":
    rule pull_gisaid_sequences:
        output:
            meta=temp(dynamic(f"{fasta_output_prefix}/{{sample}}.tsv.gz")),
            data=temp(dynamic(f"{fasta_output_prefix}/{{sample}}.fasta"))
        threads: max_cpus
        shell:
            """
            mkdir -p /dev/shm/bjorn;
            (({is_test} && gunzip -c {test_data}) || (curl -u {username}:{password} https://www.epicov.org/epi3/3p/scripps/export/provision.json.xz | xz -d -T{max_task_cpus})) |
                    parallel --pipe --tmpdir /dev/shm/bjorn --block {chunk_size} -j{rem_cpus} \
                        'jq -cr "select((.covv_host|ascii_downcase == \\"human\\") and (.sequence|length > {min_length}) and ((.sequence|split(\\"N\\")|length) < (.sequence|length * {max_unknown_pct})))" | \
                            python/json2fasta.py -o {fasta_output_prefix}/{{#}} -g {gadm_data} -r {reference_filepath}' || true
            """
elif data_source == "alab_release":
    rule clone_alab_sequences:
        output:
            meta=temp(dynamic(f"{fasta_output_prefix}/{{sample}}.tsv.gz")),
            data=temp(dynamic(f"{fasta_output_prefix}/{{sample}}.fasta"))
        threads: 1
        shell:
            """
            echo {fasta_output_prefix}
            git clone https://github.com/andersen-lab/HCoV-19-Genomics.git
            cp HCoV-19-Genomics/consensus_sequences/*.fasta {fasta_output_prefix} 
            python/manipulate_metadata.py -i HCoV-19-Genomics/metadata.csv -o {fasta_output_prefix}
            for file in {fasta_output_prefix}/*.fasta;do
                cat {reference_filepath} > "$file"
            done
            """
else:
    print(f'Error: data_source should be "gisaid_feed" or "alab_release" -- got {data_source}')
    sys.exit()
