# ./tsv2fasta.sh /compdata/ <input.tsv >output.fasta
compdata=`realpath $1` lvl=`cat $compdata/level`
pfx=`cat $compdata/prefix` dct="$compdata/dictionary"
cd `dirname $0`
awk 'BEGIN{FS="\t"; OFS="\t"}{ print(">"$1, "'$pfx'"$5) }' | ./uzipl.py -d $dct | tr '\t' '\n'
