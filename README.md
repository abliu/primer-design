# PCR primer design tool

This pipeline (`gen_filter_primers.sh`) designs qPCR primers for a given set of input genomes, one for each genome. Each primer pair is specific to one genome (versus all others). The script chains together scripts for generating primers (`gen_primers.py`), checking cross-amplifications (`bowtie2 commands`), and filtering out non-specific primers (`filter_primers.R`).

This builds upon an existing set of tools for primer design, such as [SpeciesPrimer](https://peerj.com/articles/8544), [fdp](https://github.com/widdowquinn/find_differential_primers), and [RUCS](https://bitbucket.org/genomicepidemiology/rucs/src/master/). Our tool differs from these previous tools:

1. We achieve significant speed gains by replacing `primersearch` with `bowtie2`, repurposing the `bowtie2` aligner to identify cross-amplifications.
2. We focus on the use case of generating one primer per genome rather than per set of genomes; this allows us to simplify our pipeline by avoiding steps like identifying conserved genes.

## Dependencies

See `requirements.txt`.

## Usage

See `gen_filter_primers.sh`. The pipeline expects input genomes to be separate fasta files in one directory; file names should describe the sample/genome uniquely. This directory is the input to the first step of the pipeline.
