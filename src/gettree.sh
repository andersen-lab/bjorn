#!/bin/bash

wget https://hgwdev.gi.ucsc.edu/~angie/b7f72f3423a52f4dfdc9fc37a09a9347/gisaidAndPublic.latest.masked.pb.gz -O gisaidAndPublic.latest.masked.pb.gz
rm gisaidAndPublic.latest.masked.pb
gunzip gisaidAndPublic.latest.masked.pb.gz
matUtils summary -i gisaidAndPublic.latest.masked.pb -C sample_clades
matUtils extract -i gisaidAndPublic.latest.masked.pb -p <(grep -iE $'.+\t.+\t.*(\_|prop|misc).*' sample_clades) -o tree.prune.pb
mv tree.prune.pb ../data/tree.prune.latest.pb
#wget https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json
#cp alias_key.json ../data/
wget https://raw.githubusercontent.com/outbreak-info/outbreak.info/master/curated_reports_prep/lineages.yml
mv lineages.yml ../data
