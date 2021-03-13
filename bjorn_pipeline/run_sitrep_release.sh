#!/bin/bash
cd /home/al/code/bjorn/bjorn_pipeline
source ~/anaconda3/etc/profile.d/conda.sh
conda activate pipeline
snakemake --cores 15 --config username=***REMOVED*** password=***REMOVED*** --use-conda
