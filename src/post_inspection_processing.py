"""
Replaces 
Combine msa/*_white.fa and *_inspect.fa
# Convert aligned to unaligned fasta
bash ~/code/bjorn/scripts/convert_to_unaigned_fasta.sh 2021-08-25_release_aligned_combined.fa > 2021-08-25_release_unaligned_combined.fa
# Convert multifastas to sindividual fasta preferably in a subfolder
bash ~/code/bjorn/scripts/multifasta_to_fasta.sh ../2021-08-25_release_unaligned_combined.fa

with a single step process
"""

