#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh
conda activate bjorn
snakemake --configfile $1 clear
snakemake --cores 20 --configfile $1 all upload
