"""
* Filename: sequence_segments.py
* Part of: 1KGP_sweeps_enformer_project
* Contains: Script to split FASTA files into overalapping subregions.
* Debbie Lilly, July 2023
"""

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import math
import sys

# Replace with location of your FASTA file.
genome_record = SeqIO.parse(
    "fasta_files/complete_genome_records_hc2.fa", "fasta")
# Replace with location for output FASTA file.
file_out = "../data/hg38_fasta_files/segmented_seq_chr2.fa"
# Clearing contents of file_out to ensure it's empty.
with open(file_out, 'w') as f_out:
    SeqIO.write('', f_out, 'fasta')

# Can specify your own seq lengths and increment as command line args.
windowSize = int(sys.argv[1]) if (len(sys.argv) > 1) else 393216
step = int(sys.argv[2]) if (len(sys.argv) > 1) else 57344
overlap = windowSize - step

# Input: Single FASTA entry.
# Returns: The start position of the FASTA entries
#          (should be the same for all FASTA entries).
# Does: Returns the start position of your FASTA range, needed for split() function.
def get_start_pos(sequence):
    id = sequence.id
    start = id.index(":") + 1
    end = id.index("-")
    return int(id[start:end])

# Input: Individual FASTA entry, 393216bp subregion as string,
#        start position of subregion as int, end position of subregion as int.
# Returns: None.
# Does: Creates FASTA entry for subregion and appends it to file_out.
# Note: Sequence range of subregion specified in description field of FASTA entry.
def writeFasta(seqObj, seq, start, end):
    fasta_info = SeqRecord(
        Seq(seq), id=seqObj.id, description=("chr2:" + str(start) + "-" + str(end)))
    with open(file_out, 'a') as f_out:
        SeqIO.write(fasta_info, f_out, 'fasta')

# Input: Single FASTA entry.
# Returns: None.
# Does: Splits given FASTA entry into numRegions number of smaller entries.
# Note: Converts idvSeq's sequence into a list to extract 393216bp sequence. 
#       (Faster than trying to index strings in Python).
def split(idvSeq):
    lst = list()
    lst.extend(str(idvSeq.seq))
    numRegions = math.ceil((len(lst) - overlap)/float(step))
    for i in range(int(numRegions)):
        seqStart = get_start_pos(idvSeq)
        startPos = step * i
        endPos = startPos + windowSize
        segment = lst[startPos:endPos]
        strSeq = ''.join(segment)
        writeFasta(idvSeq, strSeq, seqStart + startPos, seqStart + startPos + len(segment) - 1)

# Applying split() to all FASTA entries.
map(split, genome_record)
