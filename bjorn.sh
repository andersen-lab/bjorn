#!/usr/bin/env bash

configfile=$1

function conf {
  (( $# > 0 )) && printf '"%q" ' "$(jq -j ".$1" $configfile)" && conf "${@:2}"
}

data=$(eval echo `[ ! -z $3 ] || conf 'datadir'`)
temp=$(eval echo `[ ! -z $4 ] || conf 'tempdir'`)
delete_mode=`$5 || conf 'delete_mode'`

cat "$2" | eval "./readseqs.sh $(conf 'provision_decoder' 'provision_parser' 'zinfo_dir' 'datadir' 'reader_workgroups' 'reader_workersper')" > "$temp/provision.tsv" && \
./fastdiff.sh <(cut -f-5 "$data/bjorn.tsv") "$temp/provision.tsv" "$temp/deleted_seqs.tsv" "$temp/new_seqs.tsv" "$temp/difftemp" && \
([ ! $delete_mode ] || (echo "" > "$temp/deleted_seqs.tsv")) && \
eval "./analysis.sh $temp/new_seqs.tsv $(conf 'analysis_workgroups' 'analysis_workersper' 'analysis_blocksize' 'treeinfo_dir') $temp/anatemp" > "$data/bjorn.new.tsv" && \
sort -mu "$data/bjorn.new.tsv" <(grep -vFf "$temp/deleted_seqs.tsv" "$data/bjorn.tsv") > "$data/bjorn.swap.tsv" && \
mv "$data/bjorn.swap.tsv" "$data/bjorn.tsv" && \
rm -rf "$temp/provision.tsv" "$temp/deleted_seqs.tsv" "$temp/new_seqs.tsv" "$temp/difftemp" "$temp/anatemp"


#parallel -j100 --block 10M --pipepart "./norm_jsonl_output.py -i /dev/stdin -o /dev/stdout -u 'None' -g data/geo.jsonl" :::: /data/gisaid_1114.ana.tsv | pv --line-mode | gz
