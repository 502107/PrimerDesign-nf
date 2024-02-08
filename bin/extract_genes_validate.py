#!/pythonloc
import sys
from Bio import SeqIO
import os
import re

def extract_genes(ref_fasta, ref_annotation, gene_list, output_dir, padding=500, max_padding=1000):
    with open(gene_list) as f:
        # We need to consider each gene in separate lines
        gene_list = f.read().splitlines()
        gene_names = [gene.strip() for gene in gene_list]
        
    # Before attempting to obtain the gene sequence, we need to sort the reference fasta file first
    sorted_fasta_file = sort_fasta(ref_fasta)
    sequences = {record.id: str(record.seq) for record in SeqIO.parse(sorted_fasta_file, "fasta")}
    out_sequences = {}
    out_padding = {}
    
    for gene_name in gene_names:
        gene_info = find_gene(ref_annotation, gene_name)
        if gene_info:
            gene_id, start, end = gene_info
            unique = False
            
            # Change amount of padding depending on gene length, so that at most we have 2500 bp's per gene
            # (assuming the gene selection remains the same, where it only considers genes with length < 2000)
            # We're now also checking that the sequence is unique, and will probably need to create multiple primer
            # pairs for the same gene, if the sequence is huge ...
            
            # Resetting the padding for each gene
            if 1500 <= end - start:
                current_padding = 250
            else:
                current_padding = padding
                
            while not unique and current_padding <= max_padding:
                padded_seq = pad_sequence(sequences[gene_id], start, end, current_padding)
                unique = check_unique(padded_seq, sequences, gene_id)
                
                if not unique:
                    current_padding += 500
                
                if unique or current_padding > max_padding:
                    out_sequences[gene_name] = padded_seq
                    out_padding[gene_name] = current_padding
                    break
        
    return out_sequences, out_padding
        
def check_unique(padded_seq, sequences, gene_id):
    for seqid, seq in sequences.items():
        if seqid != gene_id and padded_seq in seq:
            return False
    return True
        
def sort_fasta(ref_fasta):
    
    # Sorting the reference fasta file before gene sequence extraction
    with open(ref_fasta) as f:
        records = SeqIO.parse(f, "fasta")
        sorted_records = sorted(records, key=lambda x: int(re.findall(r'\d+', x.id)[0]))
        
    base = os.path.basename(ref_fasta)
    sorted_fasta_file = os.path.join(os.path.dirname(ref_fasta), "sorted_" + base)
    
    with open(sorted_fasta_file, 'w') as output:
        SeqIO.write(sorted_records, output, 'fasta')
    return sorted_fasta_file

def find_gene(ref_annotation, gene_name):
    with open(ref_annotation) as f:
        for line in f:
            if gene_name in line:
                parts = line.split()
                # Getting the gene ID, start and end positions from the line containing the gene name
                return parts[0], int(parts[3]), int(parts[4])
    return None

# Add padding to the gene (bp's defined in extract_genes())
def pad_sequence(sequence, start, end, padding):
    new_start = max(1, start - padding)
    new_end = min(end + padding, len(sequence))
    return sequence[new_start-1:new_end]

def main():
    if len(sys.argv) not in range(5,6):
        print("Usage: extract_genes.py <FASTA_FILE> <REF_FILE> <GENE_NAMES_FILE> <OUTPUT_FILE> <PADDING (default=500)>")
        sys.exit(1)
    if len(sys.argv) == 6:
        fasta_file, ref_file, gene_names_file, output_file, padding = sys.argv[1:6]
        padding=int(padding)
    else:
        fasta_file, ref_file, gene_names_file, output_file = sys.argv[1:5]
        padding=500
    
    out_sequences, out_padding = extract_genes(fasta_file, ref_file, gene_names_file, output_file, padding=padding)

    with open(output_file, 'w') as output_file:
        for gene_name, seq in out_sequences.items():
            for geneID, padding in out_padding.items():
                if geneID == gene_name:
                    print(f">{gene_name} (padding={padding})\n{seq}", file=output_file)

if __name__ == "__main__":
    main()
