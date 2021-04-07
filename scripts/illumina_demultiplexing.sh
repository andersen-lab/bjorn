#!/bin/bash

in_dir=$1
cpus=$2
sheet=$3
out_dir=$4

cd $in_dir/
ulimit -n 5000
bcl2fastq --create-fastq-for-index-reads -r $cpus -w 4 -p $cpus --barcode-mismatches 0 --output-dir out_dir --sample-sheet $sheet

