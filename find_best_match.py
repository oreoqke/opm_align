from Bio import SeqIO
from Bio.Seq import Seq
from Bio.pairwise2 import align
import time
'''This file contains functions to align proteins by sequence using BioPython. 
This works a lot slower than using blastp, so I ended up not using this.'''



def format_time(seconds):
    """
    Converts seconds to a human-readable format (hours, minutes, seconds).
    """
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return "{:02}:{:02}:{:02}".format(int(hours), int(minutes), int(seconds))

def find_best_match(target_sequence, sequences_file):
    # Read target sequence
    target_record = SeqIO.read(target_sequence, "fasta")
    target_seq = str(target_record.seq)

    # Read sequences file
    best_match = None
    best_score = float('-inf')

    for record in SeqIO.parse(sequences_file, "fasta"):
        sequence = str(record.seq)
        alignments = align.globalxx(target_seq, sequence)
        top_alignment = alignments[0]
        score = top_alignment[2]
        if score > best_score:
            best_score = score
            best_match = record

    return best_match, best_score

# Usage

target_sequence_file = "seq2.fasta"
sequences_file = "opm_fasta.txt"

start_time = time.time()
best_match, best_score = find_best_match(target_sequence_file, sequences_file)
end_time = time.time()

if best_match:
    print("Best match:")
    print("ID:", best_match.id)
    print("Sequence:", best_match.seq)
    print("Score:", best_score)
    print("Time:", format_time(end_time - start_time))
else:
    print("No match found.")
   
