from datetime import datetime
configfile:"example_config.json"
username = config['gisaid_username']
password = config['gisaid_password']
work_dir = config['work_dir']
current_datetime = config["current_datetime"] if config["current_datetime"] != False else datetime.now().strftime("%Y-%m-%d-%H-%M")
chunk_size = config['chunk_size']
min_length = int(config['min_length'])
max_unknown_pct = float(config['max_unknown_pct'])
max_task_cpus = int(config['max_task_cpus'])
max_cpus = workflow.cores
rem_cpus = max_cpus - max_task_cpus
reference_fp = "./pangolin/pangolin/data/reference.fasta"
patient_zero = "outgroup_A"
data_source = config['data_source']
geojson_prefix = config["geojson_prefix"]
min_date = config["min_date"]
gisaid_data = config["gisaid_data"]
is_manual_in = config["is_manual_in"]
gisaid_uri = config["gisaid_uri"]
unknown_value = config["unknown_value"]

fasta_output_prefix = work_dir + "/chunks_fasta_" + current_datetime

if data_source == "gisaid_feed":
    rule pull_gisaid_sequences:
        output:
            meta=temp(dynamic(f"{fasta_output_prefix}/{{sample}}.tsv")),
            data=temp(dynamic(f"{fasta_output_prefix}/{{sample}}.fasta.gz"))
        threads: max_cpus - 1
        shell:
            """
            mkdir -p {work_dir};
            mkdir -p {work_dir}/parallel;
            mkdir -p {fasta_output_prefix};
            ( ( ! ({is_manual_in}) && (curl -u {username}:{password} {gisaid_uri} | xz -d -T4) ) ||
              (gunzip -c {gisaid_data} ) ) |
                    parallel --pipe --tmpdir {work_dir}/parallel --block {chunk_size} -j4 \
                        'jq -cr " \
                            select( ( .covv_host|ascii_downcase == \\"human\\" ) \
                                and ( .sequence|length > {min_length} ) \
                                and ( ( .sequence|length * {max_unknown_pct} ) > (.sequence|split(\\"N\\")|length ) ) ) | \
                            {{  strain: .covv_virus_name, \
                                loc: .covv_location|split(\\"/\\"), \
                                date: .covv_collection_date, \
                                odate: .covv_submission_date,
                                lin: .covv_lineage, \
                                id: .covv_accession_id, \
                                seq: .sequence|split(\\"\\n\\")|join(\\"\\") }} | \
                            select( .loc|length <= 6 ) | \
                            {{  c: (.loc[1] // \\"{unknown_value}\\") | ascii_downcase, \
                                d: (.loc[2] // \\"{unknown_value}\\") | ascii_downcase, \
                                l: (.loc[3] // \\"{unknown_value}\\") | ascii_downcase }} + (.) | \
                            \\">\(.strain)\n\(.seq)\t\([.strain, .id, .date, .odate, .lin, .c, .d, .l, (.loc|join(\\"/\\"))]|join(\\"\\t\\"))\\"" | \
                            tee >(cut -f1 | gzip -c > {fasta_output_prefix}/{{#}}.fasta.gz) | \
                            cut -sf2- | sed "1s/^/strain\\taccession_id\\tdate_collected\\tdate_submitted\\tpangolin_lineage\\tcountry\\tdivision\\tlocation\\tlocstring\\n/" > {fasta_output_prefix}/{{#}}.tsv' || true
            for file in {fasta_output_prefix}/*.fasta.gz; do
                cat {reference_fp} | gzip -c >> "$file"
            done || true
            """

elif data_source == "alab_release":
    rule clone_alab_sequences:
        output:
            meta=temp(dynamic(f"{fasta_output_prefix}/{{sample}}.tsv.gz")),
            data=temp(dynamic(f"{fasta_output_prefix}/{{sample}}.fasta.gz"))
        threads: 1
        shell:
            """
            echo {fasta_output_prefix}
            mkdir -p {work_dir}
            mkdir -p {work_dir}/parallel
            mkdir -p {fasta_output_prefix}
            git clone https://github.com/andersen-lab/HCoV-19-Genomics.git
            gzip -rk HCoV-19-Genomics/consensus_sequences/*.fasta
            
            #select out the files we want
            #INPUT=HCoV-19-Genomics/metadata.csv
            #OLDIFS=$IFS
            #IFS=','
            #while read flname dob ssn tel status
            #do
            #    echo "Name : $flname"
            #    echo "DOB : $dob"
            #    echo "SSN : $ssn"
            #    echo "Telephone : $tel"
            #    echo "Status : $status"
            #done < $INPUT
            #IFS=$OLDIFS
            
            mv HCoV-19-Genomics/consensus_sequences/*.fasta.gz {fasta_output_prefix} 
            
            #parallel process from here on out in chunks
            find {fasta_output_prefix} -type f -name "*.fasta.gz" | \
            parallel --tmpdir {work_dir}/parallel -l {chunk_size} -j4 python/manipulate_metadata.py -i HCoV-19-Genomics/metadata.csv -o {fasta_output_prefix} -f {{}}
            
            for file in {fasta_output_prefix}/*.fasta.gz;do
                cat {reference_fp} | gzip -c >> "$file"
            done  
            """
else:
    print(f'Error: data_source should be "gisaid_feed" or "alab_release" -- got {data_source}')
    sys.exit()

rule align_to_reference:
    input:
        f"{fasta_output_prefix}/{{sample}}.fasta.gz"
    output:
        f"{fasta_output_prefix}/{{sample}}.mutations.csv"
    threads: max_task_cpus
    shell:
        """
        echo "Here"
        minimap2 -a -x asm5 --sam-hit-only --secondary=no -t{max_task_cpus} {reference_fp} <(gzip -dc {input}) |
            gofasta sam toMultiAlign -t{max_task_cpus} --reference {reference_fp} --trimstart 265 --trimend 29674 --trim --pad |
            python/msa_2_mutations.py -r '{patient_zero}' -d '{data_source}' -i /dev/stdin -o {output}
        """

rule merge_mutations_metadata:
    input:
        meta=f"{fasta_output_prefix}/{{sample}}.tsv.gz",
        data=f"{fasta_output_prefix}/{{sample}}.mutations.csv"
    output:
        temp(f"{fasta_output_prefix}/{{sample}}.jsonl.gz")
    threads: 1
    shell:
        """
        python/normalize_and_merge.py -i {input.data} -m {input.meta} -o {output} -u None -n {min_date}  -g {geojson_prefix} -t {current_datetime}
        echo "" | gzip - >> {output} # Add new line as delimiter between chunks
        """

rule merge_json:
    input:
        dynamic(f"{fasta_output_prefix}/{{sample}}.jsonl.gz")
    output:
        temp(f"{work_dir}/_api_data_{current_datetime}.json.gz")
    threads: 1
    shell:
        """
        cat {input} > {output}
        """

rule build_meta:
    input:
        rules.merge_json.output
    output:
        temp(f"{work_dir}/_api_metadata_{current_datetime}.json")
    threads: 1
    shell:
        """
        echo '["{current_datetime}",'"$(gzip -dc {input} | wc -l)]" |
            jq '{{date_modified: .[0], records: .[1]}}' > {output}
        cat {output}
        """

rule clear:
    shell:
        """
        rm -rf {fasta_output_prefix}
        rm -f {rules.build_meta.output}
        rm -f {rules.merge_json.output} 
        rm -f rules.all.output.meta
        rm -f rules.all.output.data
        """
rule all:
    input:
        meta=rules.build_meta.output,
        data=rules.merge_json.output
    output:
        meta=config["output_meta_fp"],
        data=config["output_data_fp"]
    shell:
        """
        cp {input.meta} {output.meta}
        cp {input.data} {output.data}
        """

