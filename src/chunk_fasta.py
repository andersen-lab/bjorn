#!/usr/bin/env python
import argparse
from Bio import SeqIO
from path import Path
import bjorn_support as bs



# if not Path.isdir(out_dir):
#     Path.mkdir(out_dir)

def chunk_fasta(fasta_filepath, ref_path, chunk_size, out_dir):
    fasta_data = SeqIO.parse(open(fasta_filepath), "fasta")
    reference_data = SeqIO.read(ref_path, 'fasta')
    for i, batch in enumerate(bs.batch_iterator(fasta_data, chunk_size)):
        filepath = out_dir/f"chunk_{i+1}.fasta"
        batch.append(reference_data)
        with open(filepath, "w") as handle:
            num_sequences = SeqIO.write(batch, handle, "fasta")
        print(f"Wrote {num_sequences} records to {filepath}")
    print(f"FASTA chunking complete. All fasta files saved in {out_dir}")
    return 0

if __name__=="__main__":
    # COLLECTING USER PARAMETERS
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta",
                            type=str,
                            required=True,
                            help="FASTA filepath containing all sequences")
    parser.add_argument("-r", "--reference",
                            type=str,
                            required=True,
                            help="FASTA filepath containing reference sequence")
    parser.add_argument("-s", "--size",
                            type=str,
                            required=True,
                            help="Maximum number of sequences in each chunked output FASTA file")
    parser.add_argument("-o", "--outdir",
                            type=str,
                            required=True,
                            help="Output directory")
    args = parser.parse_args()
    fasta_filepath = args.fasta
    ref_path = args.reference
    chunk_size = int(args.size)
    out_dir = Path(args.outdir)
    fasta_data = SeqIO.parse(open(fasta_filepath), "fasta")
    reference_data = SeqIO.read(ref_path, 'fasta')
    for i, batch in enumerate(bs.batch_iterator(fasta_data, chunk_size)):
        filepath = out_dir/f"chunk_{i+1}.fasta"
        batch.append(reference_data)
        with open(filepath, "w") as handle:
            num_sequences = SeqIO.write(batch, handle, "fasta")
        print(f"Wrote {num_sequences} records to {filepath}")
    print(f"FASTA chunking complete. All fasta files saved in {out_dir}")
