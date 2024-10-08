#!/usr/bin/env bash

source ./init.sh

parallel --pipepart --blocksize $BJORN_OUT_BLOCKSIZE "./norm_jsonl_output.py -i /dev/stdin -o /dev/stdout -g '$BJORN_GEOINFO_DIR/geo.jsonl' -u '$BJORN_UNKNOWN_VAL' | gzip" :::: "$BJORN_DATADIR/new.$BJORN_PERSISTFILE" | pv > "$BJORN_DATADIR/new.$BJORN_OUTFILE"
parallel --pipepart --regexp --recstart $'\x1f\x8b\x08\\x00' $'gunzip -c | jq -cr \'del(.pangolin_lineage_crumbs)\' | gzip -c' :::: /data/new.bjorn.jsonl.gz | pv > /data/new.bjorn_nc.jsonl.gz
parallel --pipepart --regexp --recstart $'\x1f\x8b\x08\\x00' $'gunzip -c | jq -cr \'(.mutations|=map(.mutation))|(.accession_id|=tonumber)|del(.strain, .locstring)\' | gzip' :::: /data/new.bjorn.jsonl.gz | pv > /data/new.bjorn.mutless.jsonl.gz
parallel --pipepart --regexp --recstart $'\x1f\x8b\x08\\x00' $'gunzip -c | jq -cr \'.pangolin_lineage_crumbs|=(split("|")|join(";")|";\(.)")\' | gzip' :::: /data/new.bjorn.mutless.jsonl.gz | pv > /data/new.bjorn.mutless2.jsonl.gz
parallel --pipepart --regexp --recstart $'\x1f\x8b\x08\\x00' $'gunzip -c | jq -cr \'.mutations|map(.id=(.mutation+"_"+.ref_codon+"_"+.alt_codon+"_"+.pos)|(.codon_num|=tonumber)|(.pos|=tonumber)|(.is_synonymous=(.is_synonymous=="True")))|.[]|"\(.id)\t\(.)"\'' :::: /data/new.bjorn.jsonl.gz | pv --line-mode | parsort -u > /data/new.bjorn.muts.jsonl
parallel --pipepart --regexp --recstart $'\x1f\x8b\x08\\x00' $'gunzip -c | jq -cr \'.mutations|map(.id=(.mutation+"_"+.ref_codon+"_"+.alt_codon+"_"+.pos)|(.codon_num|=tonumber)|(.pos|=tonumber)|(.is_synonymous=(.is_synonymous=="True")))|.[]|"\(.id)\"\'' :::: /data/new.bjorn.jsonl.gz | pv --line-mode | parsort | uniq -c | sed -E 's/\s+/\t/g' | cut -f2,3 | sort -k2,2 > /data/new.bjorn.mut_counts.jsonl
join -t$'\t' -2 2 /data/new.bjorn.muts.jsonl /data/new.bjorn.mut_counts.jsonl | pv --line-mode | cut -f2,3 | sed -E 's/"\}\s/", "count": /g; s/$/\}/g' | cut -f2 > /data/new.bjorn.muts_w_counts.jsonl
parallel --pipepart $'./hash.py | gzip' :::: /data/new.bjorn.muts_w_counts.jsonl | pv > /data/new.bjorn.muts_w_counts.hashid.jsonl.gz
parallel --pipepart --regexp --recstart $'\x1f\x8b\x08\\x00' $'gunzip -c | jq -cr \'"\(.accession_id), \(.date_collected), \(.country_id), \(.division_id), \(.pangolin_lineage)"\' | ./copy_unalias.py' :::: /data/new.bjorn.jsonl.gz | pv --line-mode > /data/seqs22.tsv
parallel --pipepart --regexp --recstart $'\x1f\x8b\x08\\x00' $'gunzip -c' :::: /data/new.bjorn.mutless.jsonl.gz | pv --line-mode | wc -l
echo '{"date_sequences": "'`cat /data/start.txt`'", "date_tree": "'`date -r data/tree.prune.latest.pb '+%Y-%m-%d'`'", "date_modified":"'`date '+%Y-%m-%d'`'", "records":"'`cat /data/count.txt`'"}' > /data/metadata.txt
