"""
Replaces 
Combine msa/*_white.fa and *_inspect.fa
# Convert aligned to unaligned fasta
bash ~/code/bjorn/scripts/convert_to_unaigned_fasta.sh 2021-08-25_release_aligned_combined.fa > 2021-08-25_release_unaligned_combined.fa
# Convert multifastas to sindividual fasta preferably in a subfolder
bash ~/code/bjorn/scripts/multifasta_to_fasta.sh ../2021-08-25_release_unaligned_combined.fa
# move files from repo to HCoV-19-Genomics and update readme 
with a single step process
"""

import pandas as pd

def generate_location_summary_table(metadata: str, out_file: str, location_combination_config: str) -> None:
    """
    Get original ReadME 
    Generate a file based on metadata to get the summary
    Map similar locations
    Generate a new readme for the HCoV-19-Genomics Repo
    Return that new ReadME
    """
    # read file
    df = pd.read_csv(metadata)
    # generate locations df
    locations = pd.DataFrame(df.location.value_counts())
    locations.index.rename("Location", inplace=True)
    locations.rename(columns={"location": "Number of Sequences"}, inplace=True)
    # 
    # write dataframe to a str
    locations.to_markdown(out_file)
    return

def update_and_combine_readme(metadata_md: str, readme_md: str):
    """
    Takes the normal readme, updates it with new table info
    """