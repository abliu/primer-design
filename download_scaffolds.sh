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
  # Remove sequencing files with pooled genomes.
  if [$genome != "LIB060784_GEN00271393_151_S151"] && [$genome != "LIB060784_GEN00271394_152_S152"]; then
    cp "$genome/ragtag.scaffold.fasta" "$NEW_DIR/$genome.ragtag.scaffold.fasta"
  fi
done

## Locally

# Scp whole genomes from O2.
export HOST_LAB_DIR="abl19@o2.hms.harvard.edu:/n/groups/springer"
scp -r "$HOST_LAB_DIR/andrew/ragtag_scaffolds" .

# Scp the genome's coding sequences (prokka annotations) as well.
scp $HOST_LAB_DIR/cornucopia/xa45-wgs-results/annotated-nuc/*.ffn .