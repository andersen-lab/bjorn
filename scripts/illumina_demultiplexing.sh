#!/bin/bash

d=$1
c=$2

cd $d/
ulimit -n 5000
bcl2fastq --create-fastq-for-index-reads -r $c -w 4 -p $c --barcode-mismatches 0 --output-dir output --sample-sheet SampleSheet.csv

