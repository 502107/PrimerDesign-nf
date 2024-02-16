#!/pythonloc
import os
import re
import sys
import primer3
from Bio import SeqIO
from Bio.Seq import Seq
from extract_genes_validate import find_gene
from conserved_sites import read_blast


def primer3_design(gene_id, sequence):
    primer3_result = primer3.design_primers(
        seq_args={
            'SEQUENCE_ID': gene_id,
            'SEQUENCE_TEMPLATE': sequence,
            'SEQUENCE_INCLUDED_REGION': [0, len(sequence)]
        },
        global_args={
            'PRIMER_NUM_RETURN': 10000,
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_PICK_INTERNAL_OLIGO': 1,
            'PRIMER_INTERNAL_MAX_SELF_END': 8,
            'PRIMER_MIN_SIZE': 15,
            'PRIMER_MAX_SIZE': 30,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 57.0,
            'PRIMER_MAX_TM': 63.0,
            'PRIMER_MIN_GC': 30.0,
            'PRIMER_MAX_GC': 70.0,
            'PRIMER_MAX_POLY_X': 100,
            'PRIMER_INTERNAL_MAX_POLY_X': 4,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 50.0,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_MAX_SELF_ANY': 8,
            'PRIMER_MAX_SELF_END': 3,
            'PRIMER_PAIR_MAX_COMPL_ANY': 8,
            'PRIMER_PAIR_MAX_COMPL_END': 3,
            'PRIMER_PRODUCT_SIZE_RANGE': [1000, 3000]
        }
    )

    return primer3_result

def search_mismatch(primer, sequence, max_mismatch=4, forward=True):
    """
    Return True if the primer is diverse enough ( at least 1 snp in the last 3bp )

    """
    primer_len, seq_len = len(primer), len(sequence)
    assert seq_len >= primer_len, 'primer is larger than the sequence'

    if primer in sequence: # exact match bad primer
        return False

    for i in range(seq_len - primer_len + 1):
        segment = sequence[i:i+primer_len]
        mm = sum(1 for a, b in zip(primer, segment) if a != b)
        if mm < max_mismatch:
            if forward:
                last_mm = sum(1 for a, b in zip(primer[-3:], segment[-3:]) if a!=b) # found 1 snp window
            else:
                last_mm = sum(1 for a, b in zip(primer[:3], segment[:3]) if a!=b)
            if last_mm > 0:
                return True
            else:
                return False
            
    return True

def check_host(primer, host_genome, checked_primers, forward=True):
    host_genome = str(next(SeqIO.parse(host_genome, "fasta")).seq)
    primer_rc = str(Seq(primer).reverse_complement())
    
    if primer_rc in checked_primers:
        return True
    if forward:
        if not search_mismatch(primer, host_genome): # or not search_mismatch(primer_rc, host_genome):
            checked_primers.add(primer_rc)
            return True
    else:
        if not search_mismatch(primer_rc, host_genome, forward=False):
            checked_primers.add(primer_rc)
            return True
    
    return False

def check_ref_isols(primer, checked_primers, modified_sequence, forward=True):
    # Use the modified sequence from design_primers
    primer_rc = str(Seq(primer).reverse_complement())

    if primer in checked_primers:
        return True
    if forward:
        # Check for matches in the modified sequence
        if not search_mismatch(primer, modified_sequence): # or not search_mismatch(primer_rc, modified_sequence):
            checked_primers.add(primer)
            return True
    else:
        if not search_mismatch(primer_rc, modified_sequence, forward=False): # or not search_mismatch(primer_rc, modified_sequence):
            checked_primers.add(primer)
            return True

    return False

def extract_gene_hits(all_blast_results, gene_id):
    gene_hits = {}
    
    for qseqid, hit_lists in all_blast_results.items():
        if qseqid == gene_id:
            for hits in hit_lists:
                for hit in hits:
                    sseqid = hit["sseqid"]
                    if sseqid not in gene_hits:
                        gene_hits[sseqid] = []
                    gene_hits[sseqid].append((int(hit["sstart"]), int(hit["send"])))
            
    return gene_hits

def modify_sequence(ref_fnas, gene_hits):
    modified_sequence = ""
    
    with open(ref_fnas, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            if record.id in gene_hits:
                sequence = str(record.seq)
                for start, end in gene_hits[record.id]:
                    sequence = sequence[:start-1] + "N"*(end-start+1) + sequence[end:]
                modified_sequence += sequence
            else:
                modified_sequence += str(record.seq)
        
    return modified_sequence

def design_primers(rust, conserved_sites, ref_annotation, host_genome, blast_out, ref_fnas, forward=True):
    valid_primers = []
    sequences = {}
    primer_pairs = []
    all_blast_results = {}
    checked_primers = set()
    padding_info = {}
    failed_genes = []
    
    # Here we're using the read_blast function to save the blast results
    # for each gene into a dictionary, where the key is the gene ID and
    # the value is a list of hits
    for blast in os.listdir(blast_out):
        if blast.startswith(rust):
            blast_path = os.path.join(blast_out, blast)
            gene_blast_results = read_blast(blast_path)
            
            for qseqid, hit_info in gene_blast_results.items():
                if qseqid not in all_blast_results:
                    all_blast_results[qseqid] = [hit_info]
                else:
                    all_blast_results[qseqid].append(hit_info)  
    
    with open(conserved_sites) as f:
        for record in SeqIO.parse(f, "fasta"):
            sequences[record.id] = str(record.seq).replace('-', 'N')
           
    for gene_id, sequence in sequences.items():
        gene_hits = extract_gene_hits(all_blast_results, gene_id)
        modified_sequence = modify_sequence(ref_fnas, gene_hits)

        valid_gene = False
        gene_info = find_gene(ref_annotation, gene_id)
        if not gene_info:
            continue
            
        sseqid, sstart, send = gene_info
        
        padding = 500
        
        # Replacing the gene sequence with Ns to avoid creating of primer pairs for the internal regions (+/- 100bp from the ends)
        new_sequence = ""
        sequence = sequence[:(padding - 100)] + "N"*((len(sequence) - padding + 100) - (padding - 100)) + sequence[(len(sequence) - (padding - 100)):]
        new_sequence += sequence
        
        primer3_result = primer3_design(gene_id, new_sequence)
        
        primer_pairs.clear()
        # To find the best_hit (i.e. the primer pair that covers the region as close to the
        # actual gene as possible, without the padding, we're sorting by closest distance)
        for i in range(primer3_result['PRIMER_PAIR_NUM_RETURNED']):
            primer_left_pos = primer3_result[f'PRIMER_LEFT_{i}'][0]
            primer_right_pos = primer3_result[f'PRIMER_RIGHT_{i}'][0]
            seqlen_amplifying = primer_right_pos - primer_left_pos
            primer_pairs.append((seqlen_amplifying, i))
        primer_pairs.sort(key=lambda x: x[0])
        
        for dist, i in primer_pairs:
            primer_left_seq = primer3_result[f'PRIMER_LEFT_{i}_SEQUENCE']
            primer_right_seq = primer3_result[f'PRIMER_RIGHT_{i}_SEQUENCE']
            primer_left_pos = primer3_result[f'PRIMER_LEFT_{i}'][0]
            primer_right_pos = primer3_result[f'PRIMER_RIGHT_{i}'][0]
            
            # I think the best way to do this is to check which part of the blast results
            # the primer hits amplify, then get that sseqid, sstart and send to check against
            # the reference isolates, but skip those regions
            # I think it's faster to check the isolates first and then the host (takes ~5 mins/seq)
            if forward:
                if not check_ref_isols(primer_left_seq, checked_primers, modified_sequence):
                    # Checking if the primer pair amplifies on the host genome
                    if not check_host(primer_left_seq, host_genome, checked_primers):
                        valid_primers.append({
                            'GeneID': gene_id,
                            'PrimerPair': i+1,
                            'PrimerSeq_F': primer_left_seq,
                            'PrimerLoc_F': primer_left_pos,
                            'PrimerSeq_R': primer_right_seq,
                            'PrimerLoc_R': primer_right_pos,
                            'GeneLength_no_padding': len(sequence) - 2*padding,
                            'GeneLength_with_padding': len(sequence),
                            'Primer_coverage': primer_right_pos - primer_left_pos
                        })
                        valid_gene = True
                        break

            else:
                if not check_ref_isols(primer_right_seq, checked_primers, modified_sequence, forward=False):
                    # Checking if the primer pair amplifies on the host genome
                    if not check_host(primer_right_seq, host_genome, checked_primers, forward=False):
                        valid_primers.append({
                            'GeneID': gene_id,
                            'PrimerPair': i+1,
                            'PrimerSeq_F': primer_left_seq,
                            'PrimerLoc_F': primer_left_pos,
                            'PrimerSeq_R': primer_right_seq,
                            'PrimerLoc_R': primer_right_pos,
                            'GeneLength_no_padding': len(sequence) - 2*padding,
                            'GeneLength_with_padding': len(sequence),
                            'Primer_coverage': primer_right_pos - primer_left_pos
                        })
                        valid_gene = True
                        break
                    
        if not valid_gene:
            failed_genes.append(gene_id)

    # If no valid primers are identified then just fail the gene
    for gene_id in failed_genes:
        valid_primers.append({
            'GeneID': gene_id,
            'PrimerPair': None,
            'PrimerSeq_F': None,
            'PrimerLoc_F': None,
            'PrimerSeq_R': None,
            'PrimerLoc_R': None,
            'GeneLength_no_padding': None,
            'GeneLength_with_padding': None,
            'Primer_coverage': None
        })
    return valid_primers

def main():
    
    if len(sys.argv) != 9:
        print("Usage: design_primers.py <rust> <ref_dir> <conserved_sites.fasta> <host_genome.fasta> <blast_output_dir> <rust_all.fna> <output_file.csv> <forward/reverse>")
        sys.exit(1)

    rust = sys.argv[1]
    ref_dir = sys.argv[2]
    conserved_sites = sys.argv[3]
    host_genome = sys.argv[4]
    blast_out = sys.argv[5]
    ref_fnas = sys.argv[6]
    output_file = sys.argv[7]
    direction = sys.argv[8]
    
    if direction == "forward":
        direction = True
    elif direction == "reverse":
        direction = False
    
    # This gets all the gff files in the directory (although there should only be one)
    gff_files = [os.path.join(ref_dir, file) for file in os.listdir(ref_dir) if file.endswith(".gff") and os.path.isfile(os.path.join(ref_dir, file))]
    # Select the first (and hopefully only) one
    if gff_files:
        ref_annotation = gff_files[0]
    else:
        raise FileNotFoundError("No .gff file found in the specified directory")

    valid_primers = design_primers(rust, conserved_sites, ref_annotation, host_genome, blast_out, ref_fnas, forward=direction)
        
    
    with open(output_file, 'w') as csvfile:
        csvfile.write("GeneID,PrimerPair,PrimerSeq_F,PrimerLoc_F,PrimerSeq_R,PrimerLoc_R,GeneLength_no_padding,GeneLength_with_padding,Primer_coverage\n")
        for pair in valid_primers:
            csvfile.write(f"{pair['GeneID']},{pair['PrimerPair']},{pair['PrimerSeq_F']},{pair['PrimerLoc_F']},{pair['PrimerSeq_R']},{pair['PrimerLoc_R']}, {pair['GeneLength_no_padding']}, {pair['GeneLength_with_padding']}, {pair['Primer_coverage']}\n")

if __name__ == "__main__":
    main()
