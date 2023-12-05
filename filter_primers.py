import os
import pandas as pd
import time

from shared_vars import PRIMERS_DIR, OFF_TARGET_DIR, BOWTIE_IDX_DIR
# Might be a shell file actually. desired output: just the primers that don't
# cross-amplify. Might just be a bowtie and samtools run.

# Filter SAM files for unmapped, specific primers. Take in a primer alignment
# file (SAM) and output a csv of the unmapped primers only.
def filter_primers_sam(primer_sam_file):
    # We use columns 0 (primer sequence ID), 1 (SAM flag indicating alignment or
    # not), 2 (reference sequence ID) and 9 (sequence) of the SAM file.
    alignments = pd.read_csv(primer_sam_file, sep="\t", header=0,
                             names=["primer_id","sam_flag",
                                    "ref_genome_contig_id","primer_seq"],
                             usecols=[0,1,2,9], comment="@")
    alignments = alignments.assign(primer_genome = pd.DataFrame(alignments[
        'primer_id'].str.split(":", expand=True))[0],
                                   ref_genome = pd.DataFrame(alignments[
        'ref_genome_contig_id'].str.split(":", expand=True))[0]) # TODO: NO COPY PASTE!
    # SAM flag of 4 means read is unmapped to off-target genomes; this means primer
    # is specific.
    specific_primers = alignments[alignments.sam_flag == 4]
    # To instead write the above to throw out primers, we might do this. We
    # either throw out any primer for which it maps onto any genome that's
    # not itself, or we only keep primers where all the maps are unmapped or
    # onto itself. Either suggests a groupby approach.
    nonspec_primers =
    filt_primer_file = f"{primer_sam_file.split('.')[0]}_specific.csv"
    specific_primers[["primer_id","primer_seq"]].to_csv(filt_primer_file)

# Run bowtie for each genome to align its primers against other genomes.
# I am writing these together
# for primers_file in os.listdir(PRIMERS_DIR):
#     sample = primers_file.split("primers_")[1].split(".")[0]
#     off_target_genome_files = [os.path.join(OFF_TARGET_DIR, f) for f in
#                                os.listdir(OFF_TARGET_DIR) if not sample in f]
#     bowtie_idx = os.path.join(BOWTIE_IDX_DIR, f"{sample}_excl")
#     start = time.time()
#     os.system(f"bowtie2-build -f {','.join(off_target_genome_files)} "
#               f"{bowtie_idx}")
#     end = time.time()
#     print(end - start)
#     primer_sam_file = os.path.join(PRIMERS_DIR, f"{sample}_primers.sam")
#     os.system(f"bowtie2 -f -k 1 -x {bowtie_idx} -U "
#               f"{os.path.join(PRIMERS_DIR, primers_file)} -S {primer_sam_file}")
#     filter_primers_sam(primer_sam_file)

# 1. Be OK with the state of things, OR
# 2a. Confirm that primers match? Maybe write my script to check primers against
# reference genomes with SeqIO and "in."
# 2b. Check which primers are not matching and run against their individual
# indices?
(set(alignments.primer_id.unique()) -
 set(alignments[(alignments['primer_genome'] == alignments['ref_genome'])
           & (alignments.sam_flag != 4)].primer_id.unique()))
# LIB060784_GEN00271270_28_S28:reverse_primer_4
# CGTTTCCGTTGAAGTGGACG
# 2c. Check the 1k are distinct.