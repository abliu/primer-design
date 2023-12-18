#! /usr/bin/env RScript
# Convert a SAM file of primer-genome alignments into a csv of primers which do
# not have non-specific alignments (alignments to genomes besides the one that
# generated the primer).

library(optparse)
library(tidyverse)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input SAM file name", metavar="FILE"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Output csv file name", metavar="OUTFILE")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#' Reverse complement a list of DNA strings.
#' 
#' @param seqs Vector of DNA strings (character)
#' @returns Vector of reverse-complemented DNA strings (character)
#' @examples
#' reverse_complement(c("TTAGCA", "CGATAA"))
reverse_complement = function(seqs) {
  chartr("ATGC", "TACG", stringi::stri_reverse(seqs))
}

#' Convert and preprocess a SAM file of primer alignments into a data frame.
#' 
#' @param primer_sam_file Name of SAM file (character)
#' @returns Data frame where each row is an alignment of a primer to a reference
#' genome, with additional columns storing primer and reference genome names.
#' @examples
#' primers_sam2df("primers/primers.sam")
primers_sam2df = function(primer_sam_file) {
  read_tsv(primer_sam_file, col_names = F, col_select = c(1:3,10),
           comment="@") %>%
    set_names(c("primer_id","sam_flag","ref_genome_contig_id","primer_seq")) %>%
    mutate(primer_genome = str_split_i(primer_id, ":", 1),
           primer_idx = str_split_i(str_split_i(primer_id, ":", 2), "_", 3),
           ref_genome = str_split_i(ref_genome_contig_id, ":", 1),
           primer_seq = case_when(
             bitwAnd(sam_flag, 16) == 16 ~ reverse_complement(primer_seq),
             T ~ primer_seq))
}
alignments = primers_sam2df(opt$input)

# Group alignments into the set of specific primers. Each primer (uniquely
# identified by primer_id) has 0 or more alignments to each genome (uniquely
# identified by ref_genome). Only keep primers for which no alignments are
# non-specific, i.e. no alignments (i) are to genomes besides the one that
# generated the primer and (ii) have the unmapped SAM flag (bit 0x4) unset.
primers = alignments %>% group_by(primer_id, primer_seq) %>%
  filter(!any(bitwAnd(sam_flag, 4) == 0 & primer_genome != ref_genome)) %>%
  summarize(across(starts_with("primer_"), ~ first(.x)), .groups = "drop")

# Group primers into pairs and only keep the pairs where both primers passed the
# previous filter for specificity.
primer_pairs = primers %>% group_by(primer_genome, primer_idx) %>%
  filter(n() == 2) %>% ungroup(primer_idx) %>%
  filter(primer_idx == min(primer_idx)) %>% ungroup() %>%
  select(-c(primer_idx))

# Append these specific primer pairs to the existing list.
if (file.exists(opt$output)) {
  all_primer_pairs = read_csv(opt$output) %>% rbind(primer_pairs)
} else {
  all_primer_pairs = primer_pairs
}
all_primer_pairs %>% arrange(primer_genome) %>% write_csv(opt$output)
