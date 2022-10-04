#!/usr/bin/env bash
# ./tsv_commdiff.sh /old.tsv /new.tsv deletes insertions temp
comm -3 $1 $2 | tee >(grep -v $'^\t' | cut -f1 > $5) | grep $'^\t' | cut -f2- > $4 # fastest diff
grep -vFf <(cut -f1 $4) $5 > $3 # remove altered / reinserted records from deletions list
rm $4
