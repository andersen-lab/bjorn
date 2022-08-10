#!/usr/bin/env bash
TMPDIR=`realpath $4`
stamp="bj_ana_temp_"`date +"%Y%m%d%H%M"`"_"
mkdir -p $TMPDIR
cd `dirname $0`
parallel -j$2 --lb --pipepart --delay 0.1 --block 1M \
 " TMPSUB=$TMPDIR/$stamp{#} && mkdir -p \$TMPSUB ; \
   ./getfasta.sh $3 | pangolin - --threads 1 \
     --analysis-mode fast --expanded-lineage -o \$TMPSUB \
     --no-temp --tempdir \$TMPSUB > /dev/stderr; \
   gofasta sam variants -g $3/reference.gb -t1 <\$TMPSUB/mapped.sam | tail -n +2 | \
   join <(tail -n +2 \$TMPSUB/lineage_report.csv | cut -d',' -f1,2,5,17) - -t',' | \
   tr ',' $'\t' > $TMPDIR/$stamp{#}.tsv" :::: $1
sort --parallel $2 -m -k1,1 $TMPDIR/$stamp*.tsv | join $1 - -t$'\t'