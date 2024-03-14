#!/pythonloc
import os
import re
import sys
import primer3
from Bio import SeqIO
from Bio.Seq import Seq
from extract_genes_validate_splt import find_gene
from conserved_sites_splt import read_blast


def primer3_design(gene_id, sequence):
    primer3_result = primer3.design_primers(
        seq_args={
            'SEQUENCE_ID': gene_id,
            'SEQUENCE_TEMPLATE': sequence,
            'SEQUENCE_INCLUDED_REGION': [0, len(sequence)]
        },
        global_args={
            'PRIMER_NUM_RETURN': 20000,
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
            'PRIMER_PRODUCT_SIZE_RANGE': [1500, 3000]
        }
    )

    return primer3_result

def search_mismatch(primer, sequence, max_mismatch=3, forward=True):
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
        
    return modified_sequence.upper()

def design_splt_primers(rust, conserved_sites, ref_annotation, host_genome, blast_out, ref_fnas):
    valid_primers = []
    sequences = {}
    primer_pairs = []
    all_blast_results = {}
    checked_primers = set()
    padding_info = {}
    failed_genes = []
    
    for blast in os.listdir(blast_out):
        if blast.startswith(rust):
            blast_path = os.path.join(blast_out, blast)
            gene_blast_results = read_blast(blast_path)
            
            for qseqid, hit_info in gene_blast_results.items():
                if qseqid not in all_blast_results:
                    all_blast_results[qseqid] = [hit_info]
                else:
                    all_blast_results[qseqid].append(hit_info)

    # need to do something about BLAST results that are identical for different genes..this will fix the issue for now, but will probably also need to check the start/end positions of the hits
    # because we might have some genes that exist in the same scaffold...
    # After the all_blast_results dictionary is fully populated
    sseqid_dict = {}
    to_remove = set()

    # First pass to find sseqid that exist in more than one qseqid
    for qseqid, hit_infos in all_blast_results.items():
        for hit_info_list in hit_infos:
            for hit_info in hit_info_list:
                sseqid = hit_info['sseqid']
                if sseqid in sseqid_dict and sseqid_dict[sseqid] != qseqid:
                    to_remove.add(sseqid)
                else:
                    sseqid_dict[sseqid] = qseqid

    # Second pass to create a new dictionary excluding the sseqid in to_remove
    new_all_blast_results = {}
    for qseqid, hit_infos in all_blast_results.items():
        new_hit_infos = []
        for hit_info_list in hit_infos:
            if all(hit_info['sseqid'] not in to_remove for hit_info in hit_info_list):
                new_hit_infos.append(hit_info_list)
        if new_hit_infos:
            new_all_blast_results[qseqid] = new_hit_infos

    with open(conserved_sites) as f:
        for record in SeqIO.parse(f, "fasta"):
            sequences[record.id] = str(record.seq).replace('-', 'N')
           
    for gene_id, sequence in sequences.items():
        gene_hits = extract_gene_hits(new_all_blast_results, gene_id)
        modified_sequence = modify_sequence(ref_fnas, gene_hits)

        valid_gene = False
        gene_info = find_gene(ref_annotation, gene_id)
        if not gene_info:
            continue
            
        sseqid, sstart, send = gene_info
        
        padding = 3000
        
        primer3_result = primer3_design(gene_id, sequence)
        
        primer_pairs.clear()

        for i in range(primer3_result['PRIMER_PAIR_NUM_RETURNED']):
            primer_left_pos = primer3_result[f'PRIMER_LEFT_{i}'][0]
            primer_right_pos = primer3_result[f'PRIMER_RIGHT_{i}'][0]
            seqlen_amplifying = primer_right_pos - primer_left_pos
            primer_pairs.append((primer_left_pos, primer_right_pos, i))
        
        primer_pairs.sort(key=lambda x: x[0])
                
        for forward,reverse, i in primer_pairs:
            forpos = primer3_result[f'PRIMER_LEFT_{i}'][0]
            revpos = primer3_result[f'PRIMER_RIGHT_{i}'][0]
            forseq = primer3_result[f'PRIMER_LEFT_{i}_SEQUENCE']
            revseq = primer3_result[f'PRIMER_RIGHT_{i}_SEQUENCE']
            
            if not check_ref_isols(forseq, checked_primers, modified_sequence) and not check_host(forseq, host_genome, checked_primers):
                break        
            
        print('primer1: ',forseq, revseq, forpos, revpos)        
        resorted_pp=sorted(primer_pairs, key=lambda x: -x[1] if (x[0] > forward+1000 and x[0] < reverse and x[1] > reverse) else float('inf'))
            
        for forward2,reverse2, j in resorted_pp:
            if i == j or (forward2 > revpos) or (abs(forpos - forward2) <= 500) or (reverse2 <= revpos+500):
                continue
            if forward2 < revpos:
                forpos2 = primer3_result[f'PRIMER_LEFT_{j}'][0]
                revpos2 = primer3_result[f'PRIMER_RIGHT_{j}'][0]
                forseq2 = primer3_result[f'PRIMER_LEFT_{j}_SEQUENCE']
                revseq2 = primer3_result[f'PRIMER_RIGHT_{j}_SEQUENCE']

                if not check_ref_isols(revseq2, checked_primers, modified_sequence, forward=False) and not check_host(revseq2, host_genome, checked_primers, forward=False):
                    print('primer2: ',forseq2, revseq2, forpos2, revpos2)
                    valid_primers.append({
                        'GeneID': gene_id,
                        'PrimerPair1': i+1,
                        'PrimerSeq_F1': forseq,
                        'PrimerLoc_F1': forpos,
                        'PrimerSeq_R1': revseq,
                        'PrimerLoc_R1': revpos,
                        'PrimerPair2': j+1,
                        'PrimerSeq_F2': forseq2,
                        'PrimerLoc_F2': forpos2,
                        'PrimerSeq_R2': revseq2,
                        'PrimerLoc_R2': revpos2,
                        'GeneLength_no_padding': len(sequence) - 2*padding,
                        'GeneLength_with_padding': len(sequence),
                        'Primer_coverage': revpos2 - forpos
                    })
                    valid_gene = True
                    break
            
        if not valid_gene:
            failed_genes.append(gene_id)

    # If no valid primers are identified then just fail the gene
    for gene_id in failed_genes:
        valid_primers.append({
        'GeneID': gene_id,
        'PrimerPair1': None,
        'PrimerSeq_F1': None,
        'PrimerLoc_F1': None,
        'PrimerSeq_R1': None,
        'PrimerLoc_R1': None,
        'PrimerPair2': None,
        'PrimerSeq_F2': None,
        'PrimerLoc_F2': None,
        'PrimerSeq_R2': None,
        'PrimerLoc_R2': None,
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

    valid_primers = design_splt_primers(rust, conserved_sites, ref_annotation, host_genome, blast_out, ref_fnas)
        
    with open(output_file, 'w') as csvfile:
        csvfile.write("GeneID,PrimerPair1,PrimerSeq_F1,PrimerLoc_F1,PrimerSeq_R1,PrimerLoc_R1,PrimerPair2,PrimerSeq_F2,PrimerLoc_F2,PrimerSeq_R2,PrimerLoc_R1GeneLength_no_padding,GeneLength_with_padding,Primer_coverage\n")
        for pair in valid_primers:
            csvfile.write(f"{pair['GeneID']},{pair['PrimerPair1']},{pair['PrimerSeq_F1']},{pair['PrimerLoc_F1']},{pair['PrimerSeq_R1']},{pair['PrimerLoc_R1']},{pair['PrimerPair2']},{pair['PrimerSeq_F2']},{pair['PrimerLoc_F2']},{pair['PrimerSeq_R2']},{pair['PrimerLoc_R2']},{pair['GeneLength_no_padding']},{pair['GeneLength_with_padding']},{pair['Primer_coverage']}\n")

if __name__ == "__main__":
    main()

