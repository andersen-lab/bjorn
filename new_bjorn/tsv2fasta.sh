# ./tsv2fasta.sh dictionary prefix <(input.tsv) > output.fasta
set -- `realpath $1` `realpath $2`
cd `dirname $0`
awk 'BEGIN{FS="\t"; OFS="\t"}{ print(">"$1, "'`cat $2`'"$5) }' | python ./decompresslines.py --dict $1 | tr '\t' '\n'
