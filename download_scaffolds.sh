# Some useful commands for moving data:

# Copy data locally.
# 1. On the remote: create folder with copies of ragtag scaffolds so that names
# are distinguished where each file is the sample name (e.g.
# LIB060784_GEN00271243_1_S1.ragtag.scaffold.fasta) rather than a bunch of files
# named ragtag.scaffold.fasta in different folders. This enables scp and might
# be a more usable storage format for certain use cases. Don't use links because
# scp doesn't copy links.
LAB_DIR="/n/groups/springer"
cd "$LAB_DIR/amogh/analysis/xa45-redone/ragtag_scaffold"
for sample in *; do
  cp "$sample/ragtag.scaffold.fasta" "$LAB_DIR/andrew/ragtag_scaffolds/$sample.ragtag.scaffold.fasta"
done

# 2. Locally: scp.
HOST_LAB_DIR="abl19@o2.hms.harvard.edu:/n/groups/springer"
scp -r "$HOST_LAB_DIR/andrew/ragtag_scaffolds" .
# Add sample names to fasta sequence ids. Note that this edits the genome
# files in ragtag_scaffolds directly and is not idempotent (it should only be
# run once after a clean download of the genomes). There might also be a
# sed/awk 1-liner that achieves the same goal as the python script, but I
# wrote the python script for more readability.
python3 add_samples_to_seq_ids.py

# 3. Repeat both steps for annotations.
scp $HOST_LAB_DIR/cornucopia/xa45-wgs-results/annotated-nuc/*.ffn .