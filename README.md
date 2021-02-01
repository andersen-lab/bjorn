This is the code repository for `bjorn` - a suite of tools that can be used to generate information for large-scale genomic surveillance of SARS-CoV-2 sequences. `bjorn` heavily relies on external tools such as `datafunk`, `minimap2`, and `pandas`. 

## Installation
* Install Anaconda: [instructions can be found here](https://docs.anaconda.com/anaconda/install/)
* Create the `bjorn` environment
```bash
conda env create -f envs/macos.yml -n bjorn
```
* Activate environment
```bash
conda activate bjorn
```
* Install datafunk: [instructions (ensure environment is activated during installation)](https://github.com/cov-ert/datafunk)

## Usage
* Activate `bjorn` environment
```bash
conda activate bjorn
```
* Open `config.json` to specify your parameters, then save file
* NB: the config will run a test by default. Once its tested, make sure to change the `is_test` value to `false` in `config.json`
* Run `count_variants` to generate mutation information
```bash
python count_variants.py
```