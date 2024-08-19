#!/usr/bin/env bash
# ./readseqs.sh "xz -dc" ./filter.jq /compdata/ /workdir/ shards shardworkers
expand=$1 fltr=`realpath $2`
compdata=`realpath $3` lvl=`cat $compdata/level`
pfx=`cat $compdata/prefix` dct="$compdata/dictionary"
TMPDIR=`realpath $4` shards=$5 wrks=$6
pfx="s~\\\t$pfx~\\\t~g"
mkdir -p $TMPDIR && cd `dirname $0`
parallel -j$shards --pipe --tee --lb          `# branch into decompression shards` \
  " eval $expand | sed -n {}~$shards'p' |            `# each filtering for only its seqs` \
    parallel -j$wrks --pipe --round --lb          `# branch into workers for each shard` \
      ' jq --unbuffered -crf $fltr |                `# process gisaid records` \
        ./zipl.py -l $lvl -d $dct | sed $pfx ' "      `# compress sequence part` \
  ::: $( seq -s' ' 0 `echo "$shards - 1" | bc -l`) |    `# assign shard ids` \
pv --line-mode | sort --parallel $wrks `# count and sort `
