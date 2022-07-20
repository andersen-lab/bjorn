set -- `realpath $1`
cd `dirname $0`
awk 'BEGIN{FS="\t"; OFS="\t"}{ print(">"$1, "'`cat $1`'"$5) }' | python ./decompresslines.py --dict dictionary | tr '\t' '\n'
