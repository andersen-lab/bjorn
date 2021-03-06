#!/bin/bash
# TEST RUN
# python3 src/alab_release.py --out-dir /valhalla/2021-03-06_release --sample-sheet /home/al/analysis/alab_release/COVID_sequencing_summary-GISAID.csv --cpus 25 --analysis-folder /valhalla/analysis --output-metadata /home/al/code/HCoV-19-Genomics/metadata.csv
# REAL RUN
python3 src/alab_release.py --not-dry-run --include-bams --out-dir /valhalla/2021-03-06_release --sample-sheet /home/al/analysis/alab_release/COVID_sequencing_summary-GISAID.csv --cpus 25 --analysis-folder /valhalla/analysis --output-metadata /home/al/code/HCoV-19-Genomics/metadata.csv