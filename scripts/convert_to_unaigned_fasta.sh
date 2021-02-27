#!/bin/bash

fasta=$1
awk '{if($0 ~ "^>"){print $0;}else{gsub("-", "", $0);print $0;}}' $fasta


