#!/bin/bash
# This script contains commands to preprocess and download Xanthobacter and
# related genomes from the Springer lab diretory on O2 to the user's machine.
# Usage: run "On the remote" commands on O2; then run "Locally" commands on
# your machine.

## On the remote

# Copy existing ragtag scaffolded genomes (e.g.
# "LIB060784_GEN00271243_1_S1/ragtag.scaffold.fasta") to a new directory,
# adding the genome name from the directory to the filename. This file
# structure with distinctly named files enables scp. We don't use file links
# because scp doesn't copy links.
export LAB_DIR="/n/groups/springer"
export NEW_DIR="$LAB_DIR/andrew/ragtag_scaffolds"
mkdir -p $NEW_DIR
cd "$LAB_DIR/amogh/analysis/xa45-redone/ragtag_scaffold"
for genome in *; do
  cp "$genome/ragtag.scaffold.fasta" "$NEW_DIR/$genome.ragtag.scaffold.fasta"
done

## Locally

# Scp whole genomes from O2.
export HOST_LAB_DIR="abl19@o2.hms.harvard.edu:/n/groups/springer"
scp -r "$HOST_LAB_DIR/andrew/ragtag_scaffolds" .

# Add genome names to fasta sequence descriptions (e.g. change
# ">NZ_BCZQ01000001.1_RagTag" to ">LIB060784_GEN00271243_1_S1:NZ_BCZQ01000001
# .1_RagTag"). This edits the genome files in ragtag_scaffolds directly and is
# not idempotent (it should only be run once after a clean download of the
# genomes). (There might also be a sed/awk 1-liner that achieves the same goal
# as the python script, but I wrote the python script for more readability.)
python3 add_genomes_to_seq_descs.py

# Scp the genome's coding sequences (prokka annotations) as well.
scp $HOST_LAB_DIR/cornucopia/xa45-wgs-results/annotated-nuc/*.ffn .

#2. gen_primers script into one primers file (possibly doable by command line).
# note that primer3 outputs reverse primers exactly, but bowtie reverses these
# for some reason ().
#3. One bowtie command:
#bowtie2-build -f ragtag_scaffolds/*.fasta bowtie_idx/all_genomes
#bowtie2 -f -a -x bowtie_idx/all_genomes -U primers/primers.fasta -S primers/primers.sam
#4. The final filtering step/script (maybe doable by command line, but if we
# use samtools, we're relying on something anyway, and it feels better to
# just use my own script).
