# ![LOGO](figs/bjorn_logo_yellow_grey-1.png)
This is the code repository for `bjorn` - a suite of tools for processing SARS-CoV-2 sequences to support large-scale genomic surveillance. This functionality relies on external tools such as pangolin, and GNU parallel. 

## Installation
```bash
cd bjorn
docker build -t'bjorn_container' .
```

## Usage
Launching the container
```bash
docker run -it --volume [your_data_dir]:/data bjorn_container
```
Processing a data provision from GISAID's jsonl format to tsv
```bash
cat /data/gisaid.xz | ./readseqs.sh "xz -dc" ./gisaid.jq ./data /data [work_groups] [workers_per_group] [sequence_count] > /data/gisaid.tsv
```

Identifing changed records
```bash
./fastdiff.sh old_records.tsv new_records.tsv [deletes_out] [insertions_out] temp
```

Analyzing sequences (alignment and mutation- and lineage-calling)
```bash
./analysis.sh sequences.tsv [workers] [subworkers] ./data [temp] > [analysed_sequences.sh]
```

Exporting to outbreak.info's jsonl format
```bash
parallel -j[workers] --block 10M --pipepart "./norm_jsonl_output.py -i /dev/stdin -o /dev/stdout -u 'None' -g data/geo.jsonl" :::: [analysed_sequences.tsv] | gzip -c > [out.jsonl.gz]
```
