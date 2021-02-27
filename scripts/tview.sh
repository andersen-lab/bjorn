#!/bin/bash
# echo $(find ./ -name *$1*.trimmed.sorted.bam)
samtools tview --reference /home/gk/code/hCoV19/db/NC045512.fasta -p NC_045512.2:$2 $(find ./ -name *$1*.trimmed.sorted.bam)