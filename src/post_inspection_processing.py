"""
#TODO: Expand on function header
Replaces 
Combine msa/*_white.fa and *_inspect.fa
# Convert aligned to unaligned fasta
bash ~/code/bjorn/scripts/convert_to_unaigned_fasta.sh 2021-08-25_release_aligned_combined.fa > 2021-08-25_release_unaligned_combined.fa
# Convert multifastas to sindividual fasta preferably in a subfolder
bash ~/code/bjorn/scripts/multifasta_to_fasta.sh ../2021-08-25_release_unaligned_combined.fa
# move files from repo to HCoV-19-Genomics and update readme 
with a single step process
"""

import os
import sys
import subprocess
from typing import List

def concat_fastas(file_1: str, file_2: str, combined_aligned_fasta: str) -> None:
    """
    Concatenates a list of files into a combined fasta format
    #TODO: Expand on function header
    """
    with open(combined_aligned_fasta, 'w') as outfile:
        for fname in [file_1, file_2]:
            with open(fname, 'r') as infile:
                for line in infile:
                    outfile.write(line)
    return

def unalign_fasta(combined_aligned_fasta: str, combined_unaligned_fasta: str) -> None:
    """
    Takes an aligned fasta and returns an unaligned fasta at the destination required
    #TODO: Expand on function header
    """
    # if the line starts with a > (fasta header) print that line
    # else substitute "-" with nothing
    with open(combined_unaligned_fasta, 'w') as outfile:
        with open(combined_aligned_fasta, 'r') as infile:
            for line in infile:
                if line[0] == ">":
                    outfile.write(line)
                else:
                    outfile.write(line.replace("-", ""))
    return

def multifasta_to_fasta(combined_unaligned_fasta: str) -> None:
    """
    Takes a combined fasta and splits it into separate files, then copying
    back the fasta file to the main repo
    #TODO: Expand on function header
    """
    # create new directory for consensus sequences
    base_dir = os.path.dirname(combined_unaligned_fasta)
    os.mkdir(os.path.join(base_dir, "consensus_sequences"))
    # generate a separated fasta file for each sequence
    header = ""
    lines = []
    with open(combined_unaligned_fasta, 'r') as infile:
        for line in infile:
            if line[0] == ">" and header == "":
                header = get_fasta_true_name(line)
            elif line[0] != ">" and header != "":
                lines.append(line)
            elif line[0] == ">" and header != "":
                # write out the previous sequence to disk
                file_path = os.path.join(base_dir, "consensus_sequences", header[1:] + ".fa")
                with open(file_path, 'w') as outfile:
                    outfile.write(header)
                    outfile.write("".join(lines))
                header = get_fasta_true_name(line)
                lines = []
    return

def get_fasta_true_name(header: str) -> str:
    """
    Take a header line which is supposed to be a fasta name and then get the true name
    without any '/'
    """
    if '/' in header:
        return header.split('/')[2]
    else:
        return header

if __name__=="__main__":
    # generate filenames
    combined_aligned_fasta = os.path.join(
        os.path.dirname(sys.argv[1]),
        "_".join(
        os.path.basename(
            os.path.normpath(sys.argv[1])
            ).split("_")[:2] + ["combined.fa"]
        )
    )
    combined_unaligned_fasta = os.path.join(
        os.path.dirname(sys.argv[1]),
        "_".join(
        os.path.basename(
            os.path.normpath(sys.argv[1])
            ).split("_")[:2] + ["unaligned_combined.fa"]
        )
    )
    # concat, unalign, and multifasta to fasta consensus sequences
    concat_fastas(sys.argv[1], sys.argv[2], combined_aligned_fasta)
    unalign_fasta(combined_aligned_fasta, combined_unaligned_fasta)
    multifasta_to_fasta(combined_unaligned_fasta)  