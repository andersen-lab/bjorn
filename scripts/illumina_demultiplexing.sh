#!/bin/bash

d=$1

cd $d/
ulimit -n 5000
bcl2fastq --create-fastq-for-index-reads -r 4 -w 4 -p 4 --barcode-mismatches 0 --output-dir output --sample-sheet SampleSheet.csv

