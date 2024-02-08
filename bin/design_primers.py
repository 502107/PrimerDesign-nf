#!/pythonloc
import primer3
import os
import sys
from extract_genes_validate import find_gene
from conserved_sites import read_blast
from Bio import SeqIO
from Bio.Seq import Seq
import re


def primer3_design(gene_id, sequence):
    primer3_result = primer3.design_primers(
        seq_args={
            'SEQUENCE_ID': gene_id,
            'SEQUENCE_TEMPLATE': sequence,
            'SEQUENCE_INCLUDED_REGION': [0, len(sequence)]
        },
        global_args={
            # From primer3 manual:
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
            # From the p3_settings Tony and Dominik were using:
            # 'PRIMER_NUM_RETURN': 10000,
            # 'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT': 1,
            # 'PRIMER_MAX_SELF_ANY_TH': 45.00,
            # 'PRIMER_PAIR_MAX_COMPL_ANY_TH': 45.00,
            # 'PRIMER_MAX_SELF_END_TH': 35.00,
            # 'PRIMER_PAIR_MAX_COMPL_END_TH': 35.00,
            # 'PRIMER_MAX_HAIRPIN_TH': 24.0,
            # 'PRIMER_MAX_END_STABILITY': 9.0,
            # 'PRIMER_MAX_TEMPLATE_MISPRIMING': 12.00,
            # 'PRIMER_MAX_TEMPLATE_MISPRIMING_TH': 40.00,
            # 'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING': 24.00,
            # 'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH': 70.00,
            # 'PRIMER_MAX_POLY_X': 4,
            # 'PRIMER_MIN_SIZE': 15,
            # 'PRIMER_MAX_SIZE': 30,
            # 'PRIMER_MIN_GC': 30.0,
            # 'PRIMER_INTERNAL_OPT_GC_PERCENT': 50,
            # 'PRIMER_MAX_GC': 70,
            # 'PRIMER_GC_CLAMP': 1,
            # 'PRIMER_MIN_TM': 58,
            # 'PRIMER_MAX_TM': 62,
            # 'PRIMER_PAIR_MAX_DIFF_TM': 2,
            # 'PRIMER_SALT_DIVALENT': 2.0,
            # 'PRIMER_DNTP_CONC': 0.8
        }
    )

    return primer3_result

# Used to check if the primer pair amplifies on the host genome by BLASTing against
# the host genome -- VERY computationally expensive
# I'll try and implement the ipcress methods Tony was using to find mismatches instead

def search_mismatch(primer, sequence, max_mismatch=4):
    """
    Return True if the primer is diverse enough

    """
    primer_len, seq_len = len(primer), len(sequence)

    assert seq_len >= primer_len, 'primer is larger than the sequence'

    if primer in sequence: # exact match bad primer
        return False

    for i in range(seq_len - primer_len + 1):
        mm = sum(1 for a, b in zip(primer, sequence[i:i+primer_len]) if a != b)
        if mm < max_mismatch: # found a 4-snp window
            return False           
    return True


def check_host(primer, host_genome, checked_primers):
    host_genome = str(next(SeqIO.parse(host_genome, "fasta")).seq)
    primer_rc = str(Seq(primer).reverse_complement())
    primer_c = str(Seq(primer).complement())
    
    if primer_rc in checked_primers:
        return True
    
    if not search_mismatch(primer_rc, host_genome) or not search_mismatch(primer, host_genome):
        checked_primers.add(primer_rc)
        return True
    
    return False

def check_ref_isols(primer, checked_primers, modified_sequence):
    # Use the modified sequence from design_primers

    if primer in checked_primers:
        return True

    primer_c = str(Seq(primer).complement())
    primer_rc = str(Seq(primer).reverse_complement())

    # Check for matches in the modified sequence
    if not search_mismatch(primer, modified_sequence): #or not search_mismatch(primer_rc, modified_sequence):
        checked_primers.add(primer)
        return True

    return False

def extract_gene_hits(all_blast_results, gene_id):
    gene_hits = {}
    
    for qseqid, hits in all_blast_results.items():
        if qseqid == gene_id:
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

# Now I need to separate the blast results based on the sseqid
# and then figure out the position of the blast hits based on
# sstart and send.

def find_sseqid_hits(all_blast_results, primer_seqs, gene_id):
    # Will need an array here to save all the hits from the different sseqids
    hits_info = []
    # for primer_seq in primer_seqs:
    #     primer_rc = str(Seq(primer_seq).reverse_complement())
    #     sstart, send = None, None
    
        # for qseqid, hits in all_blast_results.items():
        #     if qseqid == gene_id:
        #         for hit in hits:
        #             if primer_seq in hit["sseq"]:
        #                 sstart, send = int(hit["sstart"]), int(hit["send"])
        #             elif primer_rc in hit["sseq"]:
        #                 sstart, send = int(hit["send"]), int(hit["sstart"])
        #         hits_info.append((hit["sseqid"], sstart, send))
        
    for qseqid, hits in all_blast_results.items():
        if qseqid == gene_id:
            for hit in hits:
                for primer_seq in primer_seqs:
                    primer_rc = str(Seq(primer_seq).reverse_complement())
                    if primer_seq in hit["sseq"] or primer_rc in hit["sseq"]:
                        sstart, send = int(hit["sstart"]), int(hit["send"])
                        hits_info.append((hit["sseqid"], sstart, send))
        
    return hits_info

def design_primers(rust, conserved_sites, ref_annotation, host_genome, blast_out, ref_fnas, genes_file):
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
    
    # Get padding from gene extraction
    with open(genes_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                gene_id = line.split(' ')[0][1:]
                
                padding_match = re.search(r'\(padding=(\d+)\)', line)
                if padding_match:
                    padding = int(padding_match.group(1))
                    padding_info[gene_id] = padding 
           
    for gene_id, sequence in sequences.items():
        gene_hits = extract_gene_hits(all_blast_results, gene_id)
        modified_sequence = modify_sequence(ref_fnas, gene_hits)

        valid_gene = False
        gene_info = find_gene(ref_annotation, gene_id)
        if not gene_info:
            continue
            
        sseqid, sstart, send = gene_info
        
        padding = padding_info[gene_id]
        
        primer3_result = primer3_design(gene_id, sequence)
        # To find the best_hit (i.e. the primer pair that covers the region as close to the
        # actual gene as possible, without the padding, I'll be sorting by closest distance)
        
        primer_pairs.clear()
        for i in range(primer3_result['PRIMER_PAIR_NUM_RETURNED']):
            # primer_left_seq = primer3_result[f'PRIMER_LEFT_{i}_SEQUENCE']
            # primer_right_seq = primer3_result[f'PRIMER_RIGHT_{i}_SEQUENCE']
            primer_left_pos = primer3_result[f'PRIMER_LEFT_{i}'][0]
            primer_right_pos = primer3_result[f'PRIMER_RIGHT_{i}'][0]
            if primer_left_pos > padding and primer_right_pos < len(sequence) - padding:
                continue
            seqlen_amplifying = primer_right_pos - primer_left_pos
            # distance_from_padding = abs(primer_left_pos - padding) + abs((len(sequence) - padding) - primer_right_pos)
            primer_pairs.append((seqlen_amplifying, i))

        # Ok, now sorting by distance from the padding
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
            # sseqid, sstart, send = find_sseqid_hits(all_blast_results, [primer_left_seq, primer_right_seq], gene_id)
            hits_info = find_sseqid_hits(all_blast_results, [primer_left_seq, primer_right_seq], gene_id)

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
        print("Usage: design_primers.py <rust> <ref_dir> <conserved_sites.fasta> <host_genome.fasta> <blast_output_dir> <rust_all.fna> <genes_fna> <output_file.csv>")
        sys.exit(1)

    rust = sys.argv[1]
    ref_dir = sys.argv[2]
    conserved_sites = sys.argv[3]
    host_genome = sys.argv[4]
    blast_out = sys.argv[5]
    ref_fnas = sys.argv[6]
    genes_file = sys.argv[7]
    output_file = sys.argv[8]
    
    
    # This gets all the gff files in the directory (although there should only be one)
    gff_files = [os.path.join(ref_dir, file) for file in os.listdir(ref_dir) if file.endswith(".gff") and os.path.isfile(os.path.join(ref_dir, file))]
    # Select the first (and hopefully only) one
    if gff_files:
        ref_annotation = gff_files[0]
    else:
        raise FileNotFoundError("No .gff file found in the specified directory")

    valid_primers = design_primers(rust, conserved_sites, ref_annotation, host_genome, blast_out, ref_fnas, genes_file)
        
    
    with open(output_file, 'w') as csvfile:
        csvfile.write("GeneID,PrimerPair,PrimerSeq_F,PrimerLoc_F,PrimerSeq_R,PrimerLoc_R,GeneLength_no_padding,GeneLength_with_padding,Primer_coverage\n")
        for pair in valid_primers:
            csvfile.write(f"{pair['GeneID']},{pair['PrimerPair']},{pair['PrimerSeq_F']},{pair['PrimerLoc_F']},{pair['PrimerSeq_R']},{pair['PrimerLoc_R']}, {pair['GeneLength_no_padding']}, {pair['GeneLength_with_padding']}, {pair['Primer_coverage']}\n")

if __name__ == "__main__":
    main()
