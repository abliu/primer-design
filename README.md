# PCR primer design tool

This pipeline (`gen_filter_primers.sh`) designs qPCR primers for a given set of input genomes, one for each genome. Each primer pair is specific to one genome (versus all others). The script chains together scripts for generating primers (`gen_primers.py`), checking cross-amplifications (`bowtie2` commands), and filtering out non-specific primers (`filter_primers.R`).

This builds upon an existing set of tools for primer design, such as [SpeciesPrimer](https://peerj.com/articles/8544), [fdp](https://github.com/widdowquinn/find_differential_primers), and [RUCS](https://bitbucket.org/genomicepidemiology/rucs/src/master/). Our tool differs from these previous tools:

1. We achieve significant speed gains by replacing `primersearch` with `bowtie2`, repurposing the `bowtie2` aligner to identify cross-amplifications.
2. We focus on the use case of generating one primer per genome rather than per set of genomes; this allows us to simplify our pipeline by avoiding steps like identifying conserved genes.

## Dependencies

- primer3_core, part of [primer3](https://primer3.org/manual.html), for generating primers in `gen_primers.py` (v2.6.1)
- bowtie2 and bowtie2-build, part of [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml), for aligning primers to genomes to check cross-amplification (v2.5.2)
- Python (v3.8.9):
  - [Biopython](https://biopython.org/wiki/Download) for parsing and writing FASTA files (v1.81)
  - [primer3-py](https://libnano.github.io/primer3-py/) as a Python wrapper around primer3_core, used in `gen_primers.py` (v2.0.1)
- R (v4.1.2):
  - [tidyverse](https://www.tidyverse.org/) for processing SAM file of primer-genome alignments into list of specific primer pairs (v1.3.1)
  - [optparse](https://cran.r-project.org/web/packages/optparse/index.html) for enabling R scripts to accept command line arguments (v.1.7.3)

You can probably use slightly different versions of the above and still be fine.

## Usage

See `gen_filter_primers.sh`. The pipeline expects input genomes to be separate fasta files in one directory; file names should describe the sample/genome uniquely. This directory is the input to the first step of the pipeline.
