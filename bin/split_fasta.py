#!/usr/bin/env python3

import argparse
import os
from Bio import SeqIO

def split_fasta(input_file, output_dir):
    """Split a multi-FASTA file into individual files, one per sequence."""
    os.makedirs(output_dir, exist_ok=True)

    for record in SeqIO.parse(input_file, "fasta"):
        output_file = os.path.join(output_dir, f"{record.id}.fna")
        with open(output_file, "w") as out:
            SeqIO.write(record, out, "fasta")

    return True

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split a multi-FASTA file into individual files")
    parser.add_argument("input", help="Input multi-FASTA file")
    parser.add_argument("output_dir", help="Output directory for individual sequence files")

    args = parser.parse_args()
    split_fasta(args.input, args.output_dir)
    print(f"Split completed. Files saved to {args.output_dir}")