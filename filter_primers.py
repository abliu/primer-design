import pandas as pd
# Might be a shell file actually. desired output: just the primers that don't
# cross-amplify. Might just be a bowtie and samtools run.

# Run bowtie for each genome to align its primers against other genomes.


# Filter SAM files for unmapped, specific primers.
# We use columns 0 (primer sequence ID), 1 (SAM flag indicating alignment or
# not), 2 (reference sequence ID) and 9 (sequence).
primer_file = "LIB060784_GEN00271243_1_S1.sam"
alignments = pd.read_csv(primer_file, sep="\t", header=0,
                         names=["primer_id","sam_flag","genome_contig_id",
                                "primer_seq"],
                         usecols=[0,1,2,9], comment="@")
# SAM flag of 4 means read is unmapped to off-target genomes; this means primer
# is specific.
specific_primers = alignments[alignments.sam_flag == 4]
filt_primer_file = f"{primer_file.split('.')[0]}_specific.csv"
specific_primers[["primer_id","primer_seq"]].to_csv(filt_primer_file)