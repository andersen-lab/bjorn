#!/bin/bash

fasta=$1
# Convert to single line
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $fasta > ${fasta}.tmp
# Split into individual fasta files
awk '{if( $0 ~ /^>/ ){if(s!=""){close(s);}split($0, n, "/");s=n[3]".fasta";h=$0;}else {if(s!=""){print h > s;print $0 > s}}}' ${fasta}.tmp

