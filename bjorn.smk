from datetime import datetime
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
        conda: "envs/env_min.yml"
        shell:
            """
            mkdir -p {work_dir};
            mkdir -p {work_dir}/parallel;
            mkdir -p {work_dir}/pangolin;
            pip install rapidfuzz
            mkdir -p {fasta_output_prefix};
            ( ( ! ({is_manual_in}) && (curl -u {username}:{password} {gisaid_uri} | xz -d -T8) ) ||
              (cat {gisaid_data} ) ) |
                    parallel --pipe --tmpdir {work_dir}/parallel --block {chunk_size} -j8 \
                        'jq -cr " \
                            select( ( .covv_host|ascii_downcase == \\"human\\" ) \
                                and ( .sequence|length > {min_length} ) \
                                and ( ( .sequence|length * {max_unknown_pct} ) > (.sequence|split(\\"N\\")|length ) ) ) | \
                            {{  strain: .covv_virus_name|gsub(\\" \\";\\"\\"), \
                                loc: .covv_location|split(\\"/\\"), \
                                date: .covv_collection_date, \
                                odate: .covv_subm_date,
                                lin: .covv_lineage, \
                                id: .covv_accession_id, \
                                seq: .sequence|split(\\"\\n\\")|join(\\"\\") }} | \
                            select( .loc|length <= 6 ) | \
                            \\">\(.strain)\n\(.seq)\t\([.strain, .id, .date, .odate, .lin, (.loc|join(\\"/\\"))]|join(\\"\\t\\"))\\"" | \
                            tee >(cut -f1 | gzip -c > {fasta_output_prefix}/{{#}}.fasta.gz) | \
                            cut -sf2- | sed "1s/^/strain\\taccession_id\\tdate_collected\\tdate_submitted\\tpangolin_lineage\\tlocstring\\n/" > {fasta_output_prefix}/{{#}}.tsv' || true
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
        conda: "envs/env_min.yml"
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

rule update_pangolin:
    output: "pangover"
    conda: "pangolin/environment.yml"
    shell:
        """
        cd pangolin
        pip install path
        pip install snakemake
        pip install git+https://github.com/cov-lineages/pangoLEARN.git
        pip install git+https://github.com/cov-lineages/scorpio.git
        pip install git+https://github.com/cov-lineages/constellations.git
        pip install git+https://github.com/cov-lineages/pango-designation.git
        pip install .
        pangolin --update
        pangolin --decompress-model
        cd ../
        pangolin --version > pangover
        """

rule run_pangolin:
    input:
        "pangover",
        fasta=f"{fasta_output_prefix}/{{sample}}.fasta.gz"
    output:
        analysis=temp(f"{fasta_output_prefix}/{{sample}}.analysis.csv"),
        fasta=temp(f"{fasta_output_prefix}/{{sample}}.afasta.gz")
    threads: max_task_cpus
    conda: "pangolin/environment.yml"
    shell:
        """
        tmp=`mktemp -d -p {work_dir}/pangolin`
        pangolin -t{max_task_cpus} --skip-designation-hash --tempdir $tmp --alignment --outdir $tmp {input.fasta}
        cp $tmp/lineage_report.csv {output.analysis}
        gzip -c $tmp/sequences.aln.fasta > {output.fasta}
        rm {input.fasta}
        rm -rf $tmp
        """

#rule align_to_reference:
#    input:
#        f"{fasta_output_prefix}/{{sample}}.fasta.gz"
#    output:
#        temp(f"{fasta_output_prefix}/{{sample}}.mutations.csv")
#    threads: max_task_cpus
#    conda: "envs/env_min.yml"
#    shell:
#        """
#        minimap2 -a -x asm20 --score-N=0 --sam-hit-only --secondary=no -t{max_task_cpus} {reference_fp} {input} |
#            gofasta sam toMultiAlign -t{max_task_cpus} --reference {reference_fp} --trimstart 265 --trimend 29674 --trim --pad |
#            python/msa_2_mutations.py -r '{patient_zero}' -d '{data_source}' -i /dev/stdin -o {output}
#        """

rule merge_mutations_metadata:
    input:
        meta=f"{fasta_output_prefix}/{{sample}}.tsv",
        pango=f"{fasta_output_prefix}/{{sample}}.analysis.csv",
        fasta=f"{fasta_output_prefix}/{{sample}}.afasta.gz"
    output:
        temp(f"{fasta_output_prefix}/{{sample}}.jsonl.gz")
    threads: 1
    conda: "envs/env_min.yml"
    shell:
        """
        gunzip -c {input.fasta} | python/msa_2_mutations.py -r '{patient_zero}' -d '{data_source}' -i /dev/stdin -o /dev/stdout |
            python/normalize_and_merge.py -i /dev/stdin -m {input.meta} -p {input.pango} -o {output} -u None -n {min_date}  -g {geojson_prefix} -t {current_datetime}
        echo "" | gzip - >> {output} # Add new line as delimiter between chunks
        """

rule merge_json:
    input:
        dynamic(f"{fasta_output_prefix}/{{sample}}.jsonl.gz")
    output:
        temp(f"{work_dir}/_api_data_{current_datetime}.json.gz")
    threads: 1
    conda: "envs/env_min.yml"
    shell:
        """
        gunzip -c {input} | parallel --pipe --tmpdir {work_dir}/parallel -j {max_task_cpus} --quote jq -cr 'select((.mutations | length > 0) or .pangolin_lineage == "B") | if ((has("mutations")|not) or .mutations == "") then .mutations = [] else . end' | gzip > {output}
        """

rule build_meta:
    input:
        rules.merge_json.output
    output:
        temp(f"{work_dir}/_api_metadata_{current_datetime}.json")
    threads: 1
    conda: "envs/env_min.yml"
    shell:
        """
        echo '["{current_datetime}",'"$(gzip -dc {input} | wc -l)]" |
            jq '{{date_modified: .[0], records: .[1]}}' > {output}
        cat {output}
        """

rule all:
    conda: "envs/env_min.yml"
    input:
        meta=rules.build_meta.output,
        data=rules.merge_json.output
    output:
        meta=config["output_meta_fp"],
        data=config["output_data_fp"]
    shell:
        """
        if jq -cr '.records' {input.meta} > jq -cr '.records' {output.meta}:
            cp {input.meta} {output.meta}
            cp {input.data} {output.data}
        elif
        """

rule clear:
    conda: "envs/env_min.yml"
    shell:
        """
        rm -rf {fasta_output_prefix}
        rm -f {rules.build_meta.output}
        rm -f {rules.merge_json.output}
        rm -f {rules.all.output.meta}
        rm -f pangover
        """
