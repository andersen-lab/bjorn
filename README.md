# `björn`
This is the code repository for `bjorn` - a suite of miscellaneous tools that can be used to:

* generate information for large-scale genomic surveillance of SARS-CoV-2 sequences. This functionality relies on external tools such as `datafunk`, `minimap2`, and `pandas`. 

* prepare results and data files from SARS-CoV-2 sequencing analysis for release to public databases such as GISAID, Google Cloud, and GitHub

## Installation
* Install Anaconda: [instructions can be found here](https://docs.anaconda.com/anaconda/install/)
* Create the `bjorn` environment
```bash
conda env create -f env/linux.yml -n bjorn
```
* Activate environment
```bash
conda activate bjorn
```
* Install datafunk (inside the activated environment): [instructions (ensure environment is activated during installation)](https://github.com/cov-ert/datafunk)

## Usage
### Information for Surveillance of SARS-CoV-2 Genomic Mutations
* Activate `bjorn` environment
```bash
conda activate bjorn
```
* Open `config.json` to specify your parameters such as
    * current date
    * date appended to each filepath 
    * output directory where results are saved
    * number of CPU cores available for use 
* Run the `run_sitrep.sh` script to initiate the Snakemake pipeline
```bash
bash run_sitrep.sh
```

### Post-processing of SARS-CoV-2 Sequencing Results for Release to public databases
* Activate `bjorn` environment
```bash
conda activate bjorn
```
* Open `run_alab_release.sh` to specify your parameters such as
    * filepath to sample sheet containing sample metadata (input)
    * filepath to updated metadata of samples that have already been uploaded
    * output directory where results are saved
    * number of CPU cores available for use
    * DEFAULT: test parameters
* Open `config.json` to specify your parameters such as
    * list of SARS-CoV-2 genes that are considered non-concerning
        * i.e. the occurrence of open-read frame (ORF) altering mutations can be accepted
        * e.g. ['ORF8', 'ORF10']
    * list of SARS-CoV-2 mutations that are considered non-concerning
        * i.e. the occurrence of `ORF8:Q27_` can be accepted (B117 exists)
        * e.g. ['ORF8:Q27_']
* Run the `run_alab_release.sh` script to initiate the data release pipeline
```bash
bash run_alab_release.sh
```
