#!/usr/bin/env bash
# ./parallelreadtotsv.sh /workdir/ /gisaid.xz /dictionary /prefix compression streams substreams buffsize sample
set -- `realpath $1` `realpath $2` `realpath $3` `realpath $4` $5 $6 $7 $8 $9
cd `dirname $0`
mkdir -p $1
cat $2 | parallel -j$6 --tmpdir $1 --pipe --tee --linebuffered \
  " xz -dc | sed -n "'"'"{}~$6p"'"'" | buffer -s$8K -m$8M | \
    parallel --pipe --roundrobin -j$7 --linebuffered      `# <^ parallel division of labor` \
     ' jq --unbuffered -crf ./filterextract_gisaid.jq ' | `# process gisaid records` \
    python3 ./compresslines.py --level $5 --dict $3 |     `# <v compress sequence part` \
    awk '"'BEGIN{FS="\t"; OFS="\t"}{ $1 = sprintf("%09i", substr($1, 9)); $5 = substr($5, '`cat $4 | wc -c`')}1'"'" \
  ::: $( seq -s' ' 0 `echo "$6 - 1" | bc -l`) |            # shard sequences \
head -n $9 | pv --line-mode | parsort -T $1                # count and sort
