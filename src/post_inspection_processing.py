"""
Replaces separate bash concatenation, unalignment, and multifasta to fasta
commands with a single pythonic interface

Can either be called from command line directly as a cohesive package or as
separate package elements
"""

import os
import sys
import subprocess
from typing import List


def concat_fastas(file_1: str, file_2: str, combined_aligned_fasta: str) -> None:
    """
    Concatenates a two files into a combined fasta format
    """
    with open(combined_aligned_fasta, "w") as outfile:
        for fname in [file_1, file_2]:
            with open(fname, "r") as infile:
                for line in infile:
                    outfile.write(line)
    return


def unalign_fasta(combined_aligned_fasta: str, combined_unaligned_fasta: str) -> None:
    """
    Takes an aligned fasta and writes an unaligned fasta to the same directory
    """
    # if the line starts with a > (fasta header) print that line
    # else substitute "-" with nothing
    with open(combined_unaligned_fasta, "w") as outfile:
        with open(combined_aligned_fasta, "r") as infile:
            for line in infile:
                if line[0] == ">":
                    outfile.write(line)
                else:
                    outfile.write(line.replace("-", ""))
    return


def multifasta_to_fasta(combined_unaligned_fasta: str) -> None:
    """
    Takes a combined fasta and splits it into separate files in a "consensus
    sequences" folder within the main directory
    """
    # create new directory for consensus sequences
    base_dir = os.path.dirname(combined_unaligned_fasta)
    cons_dir = os.path.join(base_dir, "consensus_sequences")
    if not os.path.exists(cons_dir):
        os.mkdir(cons_dir)
    # generate a separated fasta file for each sequence
    file_name = ""
    header = ""
    lines = []
    with open(combined_unaligned_fasta, "r") as infile:
        for line in infile:
            if line[0] == ">" and header == "":
                file_name = _get_fasta_true_name(line)
                header = line
            elif line[0] != ">" and header != "":
                lines.append(line)
            elif line[0] == ">" and header != "":
                # write out the previous sequence to disk
                file_path = os.path.join(
                    base_dir, "consensus_sequences", file_name + ".fasta"
                )
                with open(file_path, "w") as outfile:
                    outfile.writelines([header, "".join(lines)])
                file_name = _get_fasta_true_name(line)
                header = ""
                lines = []
    return


def _get_fasta_true_name(header: str) -> str:
    """
    Take a header line which is supposed to be a fasta name and then get the
    true name without any '/' or '>'
    """
    if "/" in header:
        return header.split("/")[2]
    else:
        return header[1:]


if __name__ == "__main__":
    # generate filenames
    combined_aligned_fasta = os.path.join(
        os.path.dirname(sys.argv[1]),
        "_".join(
            os.path.basename(os.path.normpath(sys.argv[1])).split("_")[:2]
            + ["aligned_combined.fa"]
        ),
    )
    combined_unaligned_fasta = os.path.join(
        os.path.dirname(sys.argv[1]),
        "_".join(
            os.path.basename(os.path.normpath(sys.argv[1])).split("_")[:2]
            + ["unaligned_combined.fa"]
        ),
    )
    # concat, unalign, and multifasta to fasta consensus sequences
    concat_fastas(sys.argv[1], sys.argv[2], combined_aligned_fasta)
    unalign_fasta(combined_aligned_fasta, combined_unaligned_fasta)
    multifasta_to_fasta(combined_unaligned_fasta)
