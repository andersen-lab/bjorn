# generate_dict.sh /workdir/ /gisaid.xz skip_sample total_sample compression /dictionary > prefix
set -- `realpath $1` `realpath $2` $3 $4 $5 `realpath $6`
cd `dirname $0`
mkdir -p $1
xz -dc $2 | sed -n "0~$3p" | pv -ls $4 | head -n $4 | jq -r '.sequence|split("\n")|join("")' | parallel -j1 "echo {} > $1/{#}.seq"
zstd -l $5 --train $1/*.seq -o $6
cat $1/*seq | python ./compresslines.py --level $5 --dict $6 | sed -e '1{h;d;}' -e 'G;s,\(.*\).*\n\1.*,\1,;h;$!d'
rm $1/*.seq
