#!/bin/bash
# convert BCL to Fastq
#PBS -l walltime=48:00:00 -l mem=46gb -l nodes=1:ppn=4 -q new -o /gpfs/home/gkarthik/logs/new.txt -j oe
directory=/gpfs/group/andersen/gkarthik/data/160812_M01244_0071_000000000-ARNGJ
out=/gpfs/home/gkarthik/analysis/2016.08.17.1471459055
date=$(date +%s)
lane=1
read_structure=251T8B251T
max_mismatches=0
min_quality_score=20
temp=/scratch/gkarthik
bar_codes=$out/_params/barcodes.txt
library_parameters=$out/_params/library_parameters.txt
#start pipeline
echo Start
date
# mkdir $temp/$date.$barcode && mkdir $temp/$date.$barcode/Data && mkdir $temp/$date.$barcode/Data/Intensities && mkdir $temp/$date.$barcode/Data/Intensities/BaseCalls && ln -s $bustard_dir/Data/Intensities/L00* $temp/$date.$barcode/Data/Intensities && ln -s $bustard_dir/Data/Intensities/BaseCalls/L00* $temp/$date.$barcode/Data/Intensities/BaseCalls && ln -s $bustard_dir/Data/Intensities/BaseCalls/barcodeData* $temp/$date.$barcode/Data/Intensities/BaseCalls
#extract barcodes
echo ExtractIlluminaBarcodes
date
java -Xmx42g -jar /gpfs/home/gkarthik/bin/picard/picard.jar ExtractIlluminaBarcodes BASECALLS_DIR=$directory/Data/Intensities/BaseCalls/ LANE=$lane READ_STRUCTURE=$read_structure BARCODE_FILE=$bar_codes METRICS_FILE=$out/_logs/barcode.metrics.$date.$lane.txt MAX_MISMATCHES=$max_mismatches MINIMUM_BASE_QUALITY=$min_quality_score NUM_PROCESSORS=4 TMP_DIR=$temp
#make bams
echo IlluminaBasecallsToSam
date
java -Xmx42g -jar /gpfs/home/andersen/bin/picard/picard.jar IlluminaBasecallsToSam BASECALLS_DIR=$directory/Data/Intensities/BaseCalls/ LANE=$lane READ_STRUCTURE=$read_structure LIBRARY_PARAMS=$library_parameters SEQUENCING_CENTER=STSI RUN_BARCODE=$date NUM_PROCESSORS=4 ADAPTERS_TO_CHECK=PAIRED_END MAX_READS_IN_RAM_PER_TILE=100000 MAX_RECORDS_IN_RAM=100000 FORCE_GC=false TMP_DIR=$temp
exit
