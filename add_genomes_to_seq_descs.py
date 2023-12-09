"""Script to add genome names to FASTA sequence descriptions.

For each file in shared_vars.TARGET_DIR, this script changes sequence
descriptions like ">NZ_BCZQ01000001.1_RagTag" (a contig name without the
genome name) to ">LIB060784_GEN00271243_1_S1:NZ_BCZQ01000001.1_RagTag" (
<genome_name>:<contig_name>). This metadata is useful generally and helps us
identify when primers align to a genome from which they were not generated (
which is a cross-amplification).
"""

from Bio import SeqIO
import argparse
import os

import shared_vars


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_dir', metavar='input-dir',
                        help="Directory with input FASTA genomes.")
    args = parser.parse_args()
    for genome_file in os.listdir(args.input_dir):
        genome = genome_file.split(".")[
            0]  # Genome name is filename pre-extension.
        file_path = os.path.join(args.input_dir, genome_file)
        # Listing out the iterator may not be efficient for larger files.
        seq_records = list(SeqIO.parse(file_path, "fasta"))
        for i, _ in enumerate(seq_records):
            seq_records[i].id = (f"{genome}{shared_vars.GENOME_CONTIG_SEP}"
                                 f"{seq_records[i].id}")
            seq_records[i].name = "";
            seq_records[i].description = ""
        with open(file_path, "w") as file:
            SeqIO.write(seq_records, file, "fasta")


if __name__ == "__main__":
    main()
