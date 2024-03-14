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
def read_blast(blast_output, pident_thresh=93, sseq_len_thresh=1000):
    valid_sseqids = {}
    checked_qseqids = set()
    all_hits = {}
    
    # First pass to identify valid sseqids and collect all hits
    with open(blast_output, 'r') as f:
        for line in f:
            parts = line.strip().split("\t")
            qseqid, sseqid, pident, bitscore, sstart, send, sseq = parts
            pident, bitscore = map(float, [pident, bitscore])
            sstart, send = map(int, [sstart, send])
            sseq_len = len(sseq)

            if qseqid not in valid_sseqids:
                valid_sseqids[qseqid] = set()

            if ((pident == 100 and sseq_len >= sseq_len_thresh) or (pident >= pident_thresh and sseq_len >= sseq_len_thresh)):
                if qseqid not in checked_qseqids:
                    checked_qseqids.add(qseqid)
                    valid_sseqids[qseqid].add(sseqid)
                    
            if sseqid in valid_sseqids[qseqid] and qseqid in checked_qseqids:          
                if qseqid not in all_hits:
                    all_hits[qseqid] = []
                all_hits[qseqid].append({'qseqid': qseqid, 'sseqid': sseqid, 'pident': pident, 'bitscore': bitscore, 'sstart': sstart, 'send': send, 'sseq': sseq})

    new_hits = {}
    for i in checked_qseqids:
        for qseqid, hits in all_hits.items():
            if i == qseqid:
                combined_hits = adjust_blast_hits(hits)
                if combined_hits:
                    new_hits[qseqid] = combined_hits

    return new_hits

def adjust_blast_hits(hits, gap_between_hits=250):
    """
    Combine hits if their combined range is within len(seq)
    """
    combined_hits = []
    # for sseqid, sseq_hits in hits.items():
    # sorted_hits = sorted(hits, key=lambda x: x['sstart'])
    current_hit = hits[0]

    for next_hit in hits[1:]:
        if ((abs(next_hit['sstart'] - current_hit['send']) < gap_between_hits) or (abs(next_hit['send'] - current_hit['sstart']) < gap_between_hits)):
            combined_hits.append(current_hit)
            current_hit = next_hit
        else:
            continue

    combined_hits.append(current_hit)
    
    return combined_hits

def align_seqs(sequences):
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_input:
        for i, seq in enumerate(sequences):
            temp_input.write(f">seq{i}\n{seq}\n")
        temp_input.flush()

        temp_output_name = temp_input.name + ".aln"
        subprocess.run(["clustalw2", "-infile=" + temp_input.name, "-outfile=" + temp_output_name, "-output=fasta"], check=True)

    return temp_output_name

def find_conserved_sites(blast_results):
    conserved_sites = {}

    for qseqid, hits in blast_results.items():
        sequences = []
        edited_sequences = []
        
        for hit_dict in hits:
            smallest_send = float('inf')
            highest_sstart = 0
            smallest_sstart = float('inf')
            highest_send = 0
            combined_seq = ''
            for hit in hit_dict:
                # Need to account for the direction (i.e. sstart/send) of the hit, when combining sequences
                if hit['send'] < hit['sstart']:
                    sorted_hits = sorted(hit_dict, key=lambda x: -x['send'])
                else:
                    sorted_hits = sorted(hit_dict, key=lambda x: x['sstart'])

            for hit in sorted_hits:
                if hit['send'] < hit['sstart']:
                    smallest_sstart = min(smallest_sstart, hit['sstart'])
                    highest_send = max(highest_send, hit['send'])

                    # Calculate the number of ambiguities to insert
                    amb_count = max(0, hit['sstart'] - highest_send - 1)
                    if combined_seq:
                        combined_seq += '-' * amb_count
                        
                    combined_seq += hit['sseq']
                    smallest_sstart = hit['sstart']
                else:
                    smallest_send = min(smallest_send, hit['send'])
                    highest_sstart = max(highest_sstart, hit['sstart'])

                    # Calculate the number of ambiguities to insert
                    amb_count = max(0, hit['sstart'] - smallest_send - 1)
                    if combined_seq:
                        combined_seq += '-' * amb_count
                        
                    combined_seq += hit['sseq']
                    smallest_send = hit['send']
                    
            sequences.append(combined_seq)

# Here, we're checking that there are enough bases in all blast results, so that we're not removing a lot of bases for hits that are empty [mod:ls120324]
        max_seqlen=max([len(seq) for seq in sequences])
        for seq in sequences:
            if len(seq) < max_seqlen - 500:
                continue
            else:
                edited_sequences.append(seq)

# Checking if we have >1 sequences after editing. Otherwise, just take all of them -- the gene will most probably fail anw.
        if edited_sequences:
            if len(edited_sequences) > 1:
                conserved_seq = find_conserved_seq(edited_sequences)
            else:
                conserved_seq = find_conserved_seq(sequences)
            if conserved_seq:
                conserved_sites[qseqid] = conserved_seq
            
    return conserved_sites

def find_conserved_seq(sequences):
    """
    Finds conserved sequences from the alignment
    """
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

def main():
    if len(sys.argv) < 5:
        raise ValueError("Usage: python conserved_sites.py <rust_type> <blast_process_output_dir> <genes.fna> <output_file>")
    
    rust_type = sys.argv[1]
    blast_process_output_dir = sys.argv[2]
    genes_seq = sys.argv[3] # Don't need it
    output_file = sys.argv[4]
    
    seqs = {record.id: str(record.seq) for record in SeqIO.parse(genes_seq, 'fasta')}
    
    all_blast_results = {}
    for file in os.listdir(blast_process_output_dir):
        if file.startswith(rust_type):
            results = read_blast(os.path.join(blast_process_output_dir, file))
            for qseqid, hit_info in results.items():
                if qseqid not in all_blast_results:
                    all_blast_results[qseqid] = [hit_info]
                else:
                    all_blast_results[qseqid].append(hit_info)
    
# For the fungicide genes, where we have very similar sequences
# I will be checking the blast results, and drop the hits that are
# repeated for another gene, as we're probably unsure of that site.

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

    all_blast_results = new_all_blast_results

# Done with the change
            
    conserved_sites = find_conserved_sites(all_blast_results)
    
    with open(output_file, 'w') as f:
        for qseqid, seq in conserved_sites.items():
            f.write(f">{qseqid}\n{seq}\n")

if __name__ == "__main__":
    main()
