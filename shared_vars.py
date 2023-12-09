"""Shared variables between primer generation and filtration scripts.
"""

# Choose a character that is not found in genome names (and preferably contig
# names) because we will split on this character.
GENOME_CONTIG_SEP = ":"