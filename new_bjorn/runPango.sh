#!/usr/bin/env bash
parallel --linebuffered --pipepart --delay 0.1 --block 1M -j48 "/bjorn/tsv2fasta.sh /bjorn/prefix | pangolin --alignment -o /out/{#} -" :::: $1
