#!/bin/bash
# This script designs qPCR primers for a given set of input genomes, one for
# each genome. Each primer pair is specific to one genome (versus all others).
# The script chains together scripts for generating primers (gen_primers.py),
# checking cross-amplifications (bowtie2 commands), and filtering out
# non-specific primers.
# Usage: $ ./gen_filter_primers.sh

# Create bowtie2 index of all genomes, which will be used to check primer
# specificity.
export BOWTIE_IDX="bowtie_idx/all_genomes"
bowtie2-build -f ragtag_scaffolds/*.fasta $BOWTIE_IDX

# Add genome names to fasta sequence descriptions (e.g. change
# ">NZ_BCZQ01000001.1_RagTag" to ">LIB060784_GEN00271243_1_S1:NZ_BCZQ01000001
# .1_RagTag"). This edits the genome files in ragtag_scaffolds directly and is
# not idempotent (it should only be run once after a clean download of the
# genomes). (There might also be a sed/awk 1-liner that achieves the same goal
# as the python script, but I wrote the python script for more readability.)
python3 add_genomes_to_seq_descs.py $GENOMES_DIR

# Generate primers.
export GENOMES_DIR="ragtag_scaffolds"
export ALL_PRIMERS_FILE="primers/primers.fasta"
export SPEC_PRIMER_FILE="primers/specific_primers.csv"
source venv/bin/activate
python3 gen_primers.py $GENOMES_DIR $ALL_PRIMERS_FILE --ignore-file $SPEC_PRIMER_FILE

# Check for cross-amplifications by running bowtie2.
export PRIMERS_SAM="primers/primers.sam"
bowtie2 -f -a -x $BOWTIE_IDX -U $ALL_PRIMERS_FILE -S $PRIMERS_SAM

# Filter out non-specific primer pairs.
Rscript filter_primers.R --input $PRIMERS_SAM --output $SPEC_PRIMER_FILE