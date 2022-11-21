# ![LOGO](figs/bjorn_logo_yellow_grey-1.png)
This is the code repository for `bjorn` - a suite of tools for processing SARS-CoV-2 sequences to support large-scale genomic surveillance. This functionality relies on external tools such as pangolin, UsHER, and GNU parallel. 

## Installation
```bash
cd bjorn
docker build -t'bjorn_container' .
```

## Usage
Launching the container
```bash
docker run --v {data_dir}:/data -v {temp_dir}:/temp -it bjorn_container
```

Running bjorn on a provision of new sequences. See example config file
```bash
./bjorn.sh {config.json} {provision.xz} [/data] [/temp]
```
(Existing sequence db in datadir will be auto-detected according to config.)

Processing a data provision from GISAID's jsonl format to tsv
```bash
cat {provision.xz} | ./readseqs.sh {provision_decoder} {provision_parser} {treeinfo_dir} {tempdir} {work_groups} {workers_per_group} > {provision.tsv}
```

Identifing changed records
```bash
./fastdiff.sh {old_records.tsv} {new_records.tsv} {deletes_out.tsv} {insertions_out.tsv} {tempdir}
```

Analyzing sequences (alignment and mutation- and lineage-calling)
```bash
./analysis.sh {provision.tsv} {workers} {subworkers} {blocksize} {treeinfo_dir} {geoinfo_dir} > {analysed_sequences.tsv}
```

Exporting to outbreak.info's jsonl format
```bash
parallel -j{workers} --block {blocksize} --pipepart "./norm_jsonl_output.py -i /dev/stdin -o /dev/stdout -u {unknown_value} -g {geoinfo_dir}" :::: {analysed_sequences.tsv} | gzip -c > {out.jsonl.gz}
```
