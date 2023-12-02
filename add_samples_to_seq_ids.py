from Bio import SeqIO
import os

import shared_vars

genomes_dir = "ragtag_scaffolds"
for sample_file in os.listdir(genomes_dir):
    sample = sample_file.split(".")[0] # Sample name is filename pre-extension.
    file_path = os.path.join(genomes_dir, sample_file)
    # Listing out the iterator may not be efficient for larger files.
    seq_records = list(SeqIO.parse(file_path, "fasta"))
    for i, _ in enumerate(seq_records):
        seq_records[i].id = (f"{sample}{shared_vars.SAMPLE_CONTIG_SEP}"
                             f"{seq_records[i].id}")
        seq_records[i].name = ""; seq_records[i].description = ""
    with open(file_path, "w") as file:
        SeqIO.write(seq_records, file, "fasta")