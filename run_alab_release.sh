#!/bin/bash
# TEST RUN
# python3 src/alab_release.py --min-coverage 95 --out-dir /valhalla/2021-03-17_release --sample-sheet /home/al/analysis/alab_release/SARS-CoV-2_sequence_tracker-GISAID.csv --cpus 25 --analysis-folder /valhalla/analysis --output-metadata /home/al/code/HCoV-19-Genomics/metadata.csv
# DEV RUN
python3 src/alab_release.py --min-coverage 95 --not-dry-run --include-bams --out-dir /valhalla/tests/test_release --sample-sheet /valhalla/tests/test_analysis/SARS-CoV-2_sequence_tracker-GISAID.csv --cpus 25 --analysis-folder /valhalla/tests/test_analysis --output-metadata /valhalla/tests/test_analysis/metadata.csv
# REAL RUN
# python3 src/alab_release.py --min-coverage 95 --not-dry-run --include-bams --out-dir /valhalla/2021-03-17_release --sample-sheet /home/al/analysis/alab_release/SARS-CoV-2_sequence_tracker-GISAID.csv --cpus 25 --analysis-folder /valhalla/analysis --output-metadata /home/al/code/HCoV-19-Genomics/metadata.csv
