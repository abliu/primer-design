"""Script to generate primers.

This script uses primer3 to design primers for each of the FASTA genomes in
args.input_dir. By default, it excludes genomes for which primers have already
been generated in args.ignore_file. It requires Biopython, pandas and
primer3-py.

This file can also be imported as a module and contains the following
functions:

    * gen_primers - generates primers for a single genome
    * main - the main function of the script
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import logging
import os
import pandas as pd
import primer3

import shared_vars


def gen_primers(genome_filename, max_num_primers=1000,
                contig_max_length=10000000):
    """Generates primers for a single genome.

    Parameters
    ----------
    genome_filename: str
        Name of FASTA file containing genome (e.g. f'{
        args.input_dir}/LIB060784_GEN00271243_1_S1.ragtag.scaffold
        .fasta')
    max_num_primers : int, optional
        Maximum number of primer pairs to generate
    contig_max_length : int, optional
        The size of the largest contig in which to search for primers,
        which can be reduced to speed up search (specify None to have no limit)

    Returns
    -------
    list[tuple[str, str]]
        A list of primer pair sequences stored in length-2 tuples (first
        sequence is forward and second sequence is reverse primer)
    """

    seqs = SeqIO.parse(genome_filename, "fasta")
    if contig_max_length is not None:
        # Iterate through contigs until we hit one shorter
        # than contig_max_length bp.
        seq = next(seqs)
        while len(seq.seq) > contig_max_length:
            try:
                seq = next(seqs)
            except StopIteration:
                seqs = [seq]
                break
    primers = []
    for seq in seqs:
        # Run primer3; use default PRIMER_PRODUCT_SIZE_RANGE but shorten if
        # sequence is shorter than minimum.
        primer3_out = primer3.bindings.design_primers({
            "SEQUENCE_ID": genome_filename,
            "SEQUENCE_TEMPLATE": str(seq.seq)
        }, {'PRIMER_PRODUCT_SIZE_RANGE': [min(len(str(seq.seq)), 100), 300]})
        # Extract primer pairs from primer3_out in list of (forward, reverse)
        # tuples.
        seq_primers = zip([primer['SEQUENCE'] for primer in
                           primer3_out['PRIMER_LEFT']],
                          [primer['SEQUENCE'] for primer in
                           primer3_out['PRIMER_RIGHT']])
        # Add to primers list and exit loop if we've found enough primers.
        primers += seq_primers
        if len(primers) >= max_num_primers:
            break
    return primers


# Convert primer lists from gen2primers into SeqRecords, for one and both
# orientations.
def _primers2records_oriented(primers, genome, is_forward):
    """Extracts forward or reverse primers from primer pairs in
    Bio.SeqRecord.SeqRecord format.
    """

    tuple_idx = 0 if is_forward else 1
    orientation = "forward" if is_forward else "reverse"
    return [SeqRecord(Seq(primer[tuple_idx]),
                      id=(f"{genome}{shared_vars.GENOME_CONTIG_SEP}"
                          f"{orientation}_primer_{idx}"),
                      name="", description="") for
            idx, primer in enumerate(primers)]


def primers2records(primers, genome):
    """Converts primer pairs from a list of sequence 2-tuples to SeqRecord list.

    Parameters
    ----------
    primers: list[tuple[str, str]]
        A list of primer pair sequences stored in length-2 tuples (first
        sequence is forward and second sequence is reverse primer)
    genome : str
        The genome name, which will be stored in the SeqRecords' descriptions

    Returns
    -------
    list[Bio.SeqRecord.SeqRecord]
        A list of the same primers converted to Bio.SeqRecord.SeqRecords with
        descriptions containing the genome name
    """

    return (_primers2records_oriented(primers, genome, True) +
            _primers2records_oriented(primers, genome, False))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_dir', metavar='input-dir',
                        help="Directory with input FASTA genomes.")
    parser.add_argument('output_file', metavar='output-file',
                        help="Name for FASTA output file with primers.")
    parser.add_argument('--ignore-file', type=str, default=None,
                        help="A csv of genomes to ignore because specific "
                             "primers have already been found for those "
                             "genomes; genome names should be in the "
                             "'primer_genome' column.")
    args = parser.parse_args()
    # If one_file=True, put all primers in one file; otherwise create one
    # primer file per genome.
    one_file = True
    target_files = os.listdir(args.input_dir)
    genome2file = {f.split(".")[0]: f for f in target_files}
    if args.ignore_file is not None:
        try:
            ignore_genomes = pd.read_csv(args.ignore_file)['primer_genome'].unique()
            [genome2file.pop(genome) for genome in ignore_genomes]
        except FileNotFoundError:
            logging.warning(f"Ignore-file {args.ignore_file} not found; "
                            f"will generate primers against all genomes "
                            f"without ignoring any")
    for i, (genome, target_file) in enumerate(genome2file.items()):
        genome_filename = os.path.join(args.input_dir, target_file)
        logging.info(f"Generating primers for {genome_filename}")
        primers = gen_primers(genome_filename)
        # Convert primers to SeqRecord objects and write to file.
        primer_records = primers2records(primers, genome)
        out_file = os.path.join(shared_vars.PRIMERS_DIR,
                                f"primers"
                                f"{'' if one_file else '_' + genome}.fasta")
        with (open(args.output_file, "a" if one_file and i > 0 else "w") as
              primers_file):
            SeqIO.write(primer_records, primers_file, "fasta")


if __name__ == "__main__":
    main()
