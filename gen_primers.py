from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import primer3

# High-level function per sample.

# Actually generate primers in these steps: a. use primer3 to generate primers
# for each CDS, b. BLAST (or primersearch?) against the other whole genomes.
# a. Use primer3. Need to pipe or put the file. Or should I use python? More
# people write python. Unless I can just pipe? I feel like Python might be
# better and more readable? OK then I should set up virtualenv and use primer3
# and SeqIO (from biopython) to process this.

# a. improvement in data storage in future is probably to put annotations and
# genomes in same folder. for now, just notice whether

# list samples
target_dir = "annotated_nucs" # as files are currently organized,
# you can directly substitute in whole genomes for coding sequences.
# sample_filename = f"{target_dir}/LIB060784_GEN00271243_1_S1.ffn"
def gen_primers(sample_filename, max_num_primers=5):
    seqs = SeqIO.parse(sample_filename, "fasta")
    primers = []
    for seq in seqs:
        # Run primer3; use default PRIMER_PRODUCT_SIZE_RANGE but shorten if
        # sequence is shorter than minimum.
        primer3_out = primer3.bindings.design_primers({
            "SEQUENCE_ID": sample_filename,
            "SEQUENCE_TEMPLATE": str(seq.seq)
        }, {'PRIMER_PRODUCT_SIZE_RANGE': [min(len(str(seq.seq)), 100), 300]})
        # Extract primer pairs from primer3_out in list of (forward, reverse)
        # tuples.
        seq_primers = zip([primer['SEQUENCE'] for primer in
                                primer3_out['PRIMER_LEFT']],
                               [primer['SEQUENCE'] for primer in
                                primer3_out['PRIMER_RIGHT']])
        # Add to primers list and exit loop if we've found enough primers.
        primers += seq_primers
        if len(primers) >= max_num_primers:
            break
    return primers

# Convert primer lists from gen2primers into SeqRecords, for one and both
# orientations.
def primers2records_oriented(primers, sample_filename, is_forward):
    tuple_idx = 0 if is_forward else 1
    orientation = "forward" if is_forward else "reverse"
    return [SeqRecord(Seq(primer[tuple_idx]),
                      id=f"{sample_filename} {orientation} primer {idx}",
                      name="", description="") for
            idx, primer in enumerate(primers)]

def primers2records(primers, sample_filename):
    return (primers2records_oriented(primers, sample_filename, True) +
            primers2records_oriented(primers, sample_filename, False))

# desired data structure: dict of {file/sample name: {PRIMER_LEFT: list,
# PRIMER_RIGHT: list}, ...}.
seq_records = []
for target_file in os.listdir(target_dir):
    print(sample_filename)
    sample_filename = f"{target_dir}/{target_file}"
    primers = gen_primers(sample_filename)
    # Convert primers to SeqRecord objects (may not need to reverse
    # complement for now).
    seq_records.extend(primers2records(primers, sample_filename))
with open("primers.fasta", "w") as primers_file:
    SeqIO.write(seq_records, primers_file, "fasta")