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

# First try at aligning, but we need it to handle MSA better
# def align_seqs(sequences):
#         aligner = Align.PairwiseAligner()
#         aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
#         aligned_seqs = []
#         # Aligning sequences to the first sequence
#         sequences_edit = [sub.replace("-", "") for sub in sequences]
#         ref_seq = sequences_edit[0]
#         for seq in sequences_edit[1:]:
#             alignments = aligner.align(ref_seq, seq)
#             aligned_seq = str(alignments[0].aligned[1])
#             aligned_seqs.append(aligned_seq)
#         aligned_seqs.insert(0, ref_seq)
#         return aligned_seqs

# Second try, should've worked but getting clustalw2 to do what i want is a pain
# def align_seqs(sequences):
#     with open("tmp.fna", "w") as f:
#         for i, seq in enumerate(sequences):
#             f.write(f">seq{i}\n{seq}\n")
#     # clustalomega_cline = ClustalwCommandline("clustalw2",infile="tmp.fna", outfile="tmp-aln.fna")
#     subprocess.run(["clustalw2", "-infile=tmp.fna", "-outfile=tmp-aln.fna"], check=True)
#     # clustalomega_cline()
#     align = AlignIO.read("tmp-aln.fna", "fasta")
#     return [str(rec.seq) for rec in align]

# Third try, basically same as second but calling clustalw2 directly
# Also, the output is an interleaved fasta file by default, so defining the output as fasta to have a sequential output
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
        
    # alignment_file_path = align_seqs(sequences)
    # with open(alignment_file_path, 'r') as alignment_file:
    #     alignment = AlignIO.read(alignment_file, "fasta")
    # summary_align = AlignInfo.SummaryInfo(alignment)
    # consensus = summary_align.dumb_consensus(ambiguous='-')
    # os.remove(alignment_file_path)
    # return str(consensus)


def find_conserved_sites(blast_results):
    conserved_sites = {}
    # The dictionary looks like that: qseqid: {sseqid: sseq},
    # so need to iterate over it twice, and use a list to get each of the sub-items
    for qseqid, hit_lists in blast_results.items():
        sequences = [hit['sseq'] for hit in hit_lists]
        
        # Now that we have the sequences in a dictionary, we need to find the conserved sites between them, but first need to convert the sequences to a list
        # for hit in hit_lists: 
        #     if 'sseq' in hit:
        #         sequences.append(hit['sseq'])
        
        if sequences:
            conserved_seq = find_conserved_seq(sequences)
            if conserved_seq:
                conserved_sites[qseqid] = conserved_seq
    return conserved_sites

# OK, let's pull it all together and see what happens
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

