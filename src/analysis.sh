#!/usr/bin/env bash

TMPDIR=`realpath $6`
stamp="bj_ana_temp_"`date '+%Y%m%d%H%M'`"_"
mkdir -p $TMPDIR
cd `dirname $0`
parallel -j$2 --block $4 --delay 1s --pipepart \
 " TMPSUB=$TMPDIR/$stamp{#} && mkdir \$TMPSUB; \
    ./getfasta.sh $5 | pangolin - --threads $3 \
     --skip-scorpio --expanded-lineage -o \$TMPSUB \
     --usher-tree $5/tree.prune.latest.pb --no-temp --tempdir \$TMPSUB > /dev/stderr; \
   gofasta sam toma -r $5/reference.fasta -t$3 <\$TMPSUB/mapped.sam | \
   parallel -j$3 --block 1M --recstart $'>' --pipe --lb --quote ./msa_2_mutations.py -i /dev/stdin -r $5/reference.fasta -o /dev/stdout | tr $'\t' ' ' | sort -k1,1 | \
   join <(tail -n +2 \$TMPSUB/lineage_report.csv | cut -d',' -f1,2,5,17) - -t',' | \
   tr $'\t,' $',\t' > $TMPDIR/$stamp{#}.tsv && rm -rf \$TMPSUB" :::: $1
sort --parallel $2 -m -k1,1 $TMPDIR/$stamp*.tsv | join $1 - -t$'\t' | sed -E $'s/^([0-9]+)\t/\\1\t'`date '+%Y-%m-%d'`$'\t/g' | cut -f-7,10-
rm -f $TMPDIR/$stamp*.tsv
