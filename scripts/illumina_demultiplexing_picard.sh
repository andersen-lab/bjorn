#!/bin/bash

indir=$1
outdir=$2
quality=$3
cpus=$4

cd $outdir/
mkdir all_fastqs
mkdir lane_1_tmp
mkdir lane_2_tmp
mkdir barcodes_1
mkdir barcodes_2
python /home/al/code/bwa_pipeline/scripts/create_picard_sheet.py --sample-sheet $indir/SampleSheet.csv --out-dir $outdir --test true
# lane 1 demultiplex
# cd $outdir/lane_1_tmp
# java -jar /home/al/code/picard/picard.jar ExtractIlluminaBarcodes \
#               BASECALLS_DIR=$indir/Data/Intensities/BaseCalls \
#               LANE=1 \
#               READ_STRUCTURE=151T8B8B151T \
#               MINIMUM_BASE_QUALITY=$quality \
#               BARCODE_FILE=$outdir/picard_experiment_barcodes.tsv \
#               NUM_PROCESSORS=$cpus \
#               OUTPUT_DIR=$outdir/barcodes_1 \
#               METRICS_FILE=$outdir/metrics_output_1.txt > $outdir/barcode_1.log
# java -Xms20G -Xmx60G -jar /home/al/code/picard/picard.jar ExtractIlluminaBarcodes \
#                      BASECALLS_DIR=/valhalla/fastq/2021.04.12/210409_A01255_0036_BH5T25DRXY/Data/Intensities/BaseCalls \
#                      LANE=1 \
#                      READ_STRUCTURE=151T8B8B151T \
#                      MINIMUM_BASE_QUALITY=20 \
#                      MAX_RECORDS_IN_RAM=100000 \
#                      BARCODE_FILE=/valhalla/ancestral_anomaly/picard_test/picard_experiment_barcodes.tsv \
#                      NUM_PROCESSORS=25 \
#                      OUTPUT_DIR=/valhalla/ancestral_anomaly/picard_test/barcodes_1 \
#                      TMP_DIR=/valhalla/al_tmp/barcodes_1 \
#                      METRICS_FILE=/valhalla/ancestral_anomaly/picard_test/barcodes_1/metrics_output_1.txt
# java -Xms20G -Xmx60G -jar /home/al/code/picard/picard.jar IlluminaBasecallsToFastq \
#       READ_STRUCTURE=151T8B8B151T \
#       BASECALLS_DIR=$indir/Data/Intensities/BaseCalls \
#       BARCODES_DIR=$outdir/barcodes_1 \
#       TMP_DIR=/valhalla/al_tmp \
#       MAX_RECORDS_IN_RAM=100000 \
#       MAX_READS_IN_RAM_PER_TILE=100000 \
#       LANE=1 \
#       NUM_PROCESSORS=$cpus \
#       IGNORE_UNEXPECTED_BARCODES=true \
#       COMPRESS_OUTPUTS=true \
#       MULTIPLEX_PARAMS=$outdir/picard_experiment_sheet.tsv \
#       RUN_BARCODE=picardrun \
#       MACHINE_NAME=NovaSeq \
#       FLOWCELL_BARCODE=abcdeACXX
# # lane 2 demultiplex
# cd $outdir/lane_2_tmp
# java -jar /home/al/code/picard/picard.jar ExtractIlluminaBarcodes \
#               BASECALLS_DIR=$indir/Data/Intensities/BaseCalls \
#               LANE=2 \
#               READ_STRUCTURE=151T8B8B151T \
#               MINIMUM_BASE_QUALITY=$quality \
#               BARCODE_FILE=$outdir/picard_experiment_barcodes.tsv \
#               NUM_PROCESSORS=$cpus \
#               OUTPUT_DIR=$outdir/barcodes_2 \
#               METRICS_FILE=$outdir/metrics_output_2.txt > $outdir/barcode_2.log
# java -Xms20G -Xmx60G -jar /home/al/code/picard/picard.jar IlluminaBasecallsToFastq \
#       READ_STRUCTURE=151T8B8B151T \
#       BASECALLS_DIR=$indir/Data/Intensities/BaseCalls \
#       BARCODES_DIR=$outdir/barcodes_2 \
#       TMP_DIR=/valhalla/al_tmp \
#       MAX_RECORDS_IN_RAM=100000 \
#       MAX_READS_IN_RAM_PER_TILE=100000 \
#       LANE=2 \
#       NUM_PROCESSORS=$cpus \
#       IGNORE_UNEXPECTED_BARCODES=true \
#       COMPRESS_OUTPUTS=true \
#       MULTIPLEX_PARAMS=$outdir/picard_experiment_sheet.tsv \
#       RUN_BARCODE=picardrun \
#       MACHINE_NAME=NovaSeq \
#       FLOWCELL_BARCODE=abcdeACXX
# rename fastq files and consolidate 
python /home/al/code/bwa_pipeline/scripts/rename_picard_fastq.py --out-dir $outdir