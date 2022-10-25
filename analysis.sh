#!/usr/bin/env bash

TMPDIR=`realpath $5`
stamp="bj_ana_temp_"`date +"%Y%m%d%H%M"`"_"
mkdir $TMPDIR
cd `dirname $0`
parallel -j$2 --block 1M --pipepart \
 " TMPSUB=$TMPDIR/$stamp{#} && mkdir \$TMPSUB; \
    ./getfasta.sh $4 | pangolin - --threads $3 \
     --analysis-mode fast --expanded-lineage -o \$TMPSUB \
     --no-temp --tempdir \$TMPSUB > /dev/stderr; \
   gofasta sam toma -r $4/reference.fasta -t$3 <\$TMPSUB/mapped.sam | \
   ./msa_2_mutations.py -i /dev/stdin -r $4/reference.fasta -o /dev/stdout | tr $'\t' ' ' | sort -k1,1 | \
   join <(tail -n +2 \$TMPSUB/lineage_report.csv | cut -d',' -f1,2,5,17) - -t',' | \
   tr $'\t,' $',\t' > $TMPDIR/$stamp{#}.tsv && rm -rf \$TMPSUB" :::: $1
sort --parallel $2 -m -k1,1 $TMPDIR/$stamp*.tsv | join $1 - -t$'\t'
rm -f $TMPDIR/$stamp*.tsv
