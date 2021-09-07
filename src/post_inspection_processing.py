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

import os
import sys
import subprocess
from typing import List

def concat_fastas(file_1: str, file_2: str, combined_aligned_fasta: str) -> None:
    """
    Concatenates a list of files into a combined fasta format
    """
    subprocess.check_call(["cat", file_1, file_2, ">", combined_aligned_fasta])
    return

def unalign_fasta(combined_aligned_fasta: str, combined_unaligned_fasta: str) -> None:
    """
    Takes an aligned fasta and returns an unaligned fasta at the destination required
    """
    subprocess.check_call(["awk", '{if($0 ~ "^>"){print $0;}else{gsub("-", "", $0);print $0;}}', combined_aligned_fasta, ">", combined_unaligned_fasta])
    return

def multifasta_to_fasta(combined_unaligned_fasta: str) -> None:
    """
    Takes a combined fasta and splits it into separate files, then copying
    back the fasta file to the main repo
    """
    os.mkdir("consensus_sequences")
    subprocess.check_call(
        ["awk",
        '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}',
        combined_unaligned_fasta,
        ">",
        combined_unaligned_fasta+".tmp"
        ]
    )
    os.chdir("consensus_sequences")
    subprocess.check_call(
        ["awk",
        '{if( $0 ~ /^>/ ){if(s!=""){close(s);}split($0, n, "/");s=n[3]".fasta";h=$0;}else {if(s!=""){print h > s;print $0 > s}}}', 
        "../" + combined_unaligned_fasta + ".tmp"
        ]
    )
    return

if __name__=="__main__":
    # generate filenames
    combined_aligned_fasta = "_".join(sys.argv[1].split("_")[:2] + ["combined.fa"])
    combined_unaligned_fasta = "_".join(sys.argv[1].split("_")[:1] + ["unaligned_combined.fa"])
    concat_fastas(sys.argv[1], sys.argv[2], combined_aligned_fasta)
    unalign_fasta(combined_aligned_fasta, combined_unaligned_fasta)
    multifasta_to_fasta(combined_unaligned_fasta)