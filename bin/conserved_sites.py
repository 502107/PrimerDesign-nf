#!/pythonloc
import os
import sys
import subprocess
import tempfile
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import AlignInfo
from Bio import Seq

# Now that we have the BLAST output files, each for a different isolate, and containing all the genes
# we need to extract the gene sequences, based on the locations specified in the BLAST output files
# and identify the conserved regions between all the isolates.

def read_blast(blast_output, pident_thresh=94, sseq_len_thresh=1000):
    results = {}
    # I will be using two scoring functions from blastn here, just to make sure that nothing funny is happening
    # - both of them should be the highest amongst the hits on the same isolate, for the same gene
    with open(blast_output, 'r') as f:
        for line in f:
            parts = line.strip().split("\t")
            qseqid, sseqid, pident, bitscore, score, sstart, send, sseq = parts
            pident, bitscore, score, sstart, send = map(float, [pident, bitscore, score, sstart, send])
            sseq_len = len(sseq)

            if (pident == 100 or (pident >= pident_thresh and sseq_len >= sseq_len_thresh)):
                if qseqid not in results or bitscore > results[qseqid]['bitscore'] or (bitscore == results[qseqid]['bitscore'] and score > results[qseqid]['score']):
                    results[qseqid] = {'sseqid': sseqid, 'pident': pident, 'bitscore': bitscore, 'score': score, 'sstart': sstart, 'send': send, 'sseq': sseq}
    return results

# The output is an interleaved fasta file by default, so defining the output as fasta to have a sequential output
def align_seqs(sequences):
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_input:
        for i, seq in enumerate(sequences):
            temp_input.write(f">seq{i}\n{seq}\n")
        temp_input.flush()

        temp_output_name = temp_input.name + ".aln"
        subprocess.run(["clustalw2", "-infile=" + temp_input.name, "-outfile=" + temp_output_name, "-output=fasta"], check=True)

    return temp_output_name

def find_conserved_seq(sequences):
    alignment_file_path = align_seqs(sequences)
    with open(alignment_file_path, 'r') as f:
        alignment = AlignIO.read(f, "fasta")
    os.remove(alignment_file_path)
    
    conserved_seq = ""
    
    for i in range(len(alignment[0])):
        col = alignment[:,i]
        if col.count(col[0]) == len(col):
            conserved_seq += col[0]
        else:
            conserved_seq += "-"
            
    return conserved_seq


def find_conserved_sites(blast_results):
    conserved_sites = {}
    # The dictionary looks like that: qseqid: {sseqid: sseq}
    for qseqid, hit_lists in blast_results.items():
        sequences = [hit['sseq'] for hit in hit_lists]
        
        # Now that we have the sequences in a dictionary, we need to find the conserved sites between them, but first need to convert the sequences to a list
        if sequences:
            conserved_seq = find_conserved_seq(sequences)
            if conserved_seq:
                conserved_sites[qseqid] = conserved_seq
    return conserved_sites

def main():
    if len(sys.argv) < 3:
        raise ValueError("Usage: python conserved_sites.py <rust_type> <blast_process_output_dir> <output_file>")
    
    rust_type = sys.argv[1]
    blast_process_output_dir = sys.argv[2]
    output_file = sys.argv[3]
    
    all_blast_results = {}
    for file in os.listdir(blast_process_output_dir):
        if file.startswith(rust_type):
            results = read_blast(os.path.join(blast_process_output_dir, file))
            for qseqid,hit_info in results.items():
                if qseqid not in all_blast_results:
                    all_blast_results[qseqid] = [hit_info]
                else:
                    all_blast_results[qseqid].append(hit_info)
                
    conserved_sites = find_conserved_sites(all_blast_results)
    
    with open(output_file, 'w') as f:
        for qseqid, seq in conserved_sites.items():
            f.write(f">{qseqid}\n{seq}\n")

if __name__ == "__main__":
    main()

