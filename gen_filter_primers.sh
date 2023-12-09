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