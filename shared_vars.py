import os

# Choose a character that is not found in sample names (and preferably contig
# names) because we will split on this character.
SAMPLE_CONTIG_SEP = ":"

# you can directly substitute in whole genomes for coding sequences.
TARGET_DIR = "ragtag_scaffolds" # as files are currently organized,
OFF_TARGET_DIR = "ragtag_scaffolds"
PRIMERS_DIR = "primers"
BOWTIE_IDX_DIR = "bowtie_idx"
os.environ['PATH'] += (os.pathsep +
                       "/Users/abl19/.local/bin/bowtie2-2.5.2-macos-arm64")