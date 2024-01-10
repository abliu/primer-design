#! /usr/bin/env RScript
# Convert a SAM file of primer-genome alignments into a csv of primers which do
# not have non-specific alignments (alignments to genomes besides the one that
# generated the primer).

library(glue)
library(optparse)
library(tidyverse)

FORWARD = "forward"
REVERSE = "reverse"

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
    set_names(c("primer_id","sam_flag","off_target_genome_contig_id","primer_seq")) %>%
    mutate(primer_genome = str_split_i(primer_id, ":", 1),
           primer_idx = str_split_i(str_split_i(primer_id, ":", 2), "_", 3),
           primer_forward = str_split_i(str_split_i(primer_id, ":", 2), "_", 1) == FORWARD,
           off_target_genome = str_split_i(off_target_genome_contig_id, ":", 1),
           primer_seq = case_when(
             bitwAnd(sam_flag, 16) == 16 ~ reverse_complement(primer_seq),
             T ~ primer_seq))
}
alignments = primers_sam2df(opt$input)

#### NEW APPROACH

# For each primer pair, get the list of off-target genomes to which it maps
# (cross-amplifies) (`summarize_primer_off_targets`). Distinguish between half
# and full matches (i.e. matches of one or both primers in the pair). We can use
# this as follows: if there is a primer pair specific to a genome (i.e. it is
# generated for that genome and has no off-target genomes), then use that primer
# pair as before. Otherwise, choose the primer pair(s) with the minimum number
# of off-target genomes.

ot_genomes = function(alignments) {
  alignments %>% distinct(off_target_genome) %>% pull(off_target_genome)
}

summarize_primer_off_targets = function(alignments_one_primer_pair, ...) {
  by_ref = alignments_one_primer_pair %>%
    # Filter for alignments which (i) have the unmapped SAM flag (bit 0x4) unset
    # and (ii) are to genomes besides the one that generated the primer. Running
    # this filter on the whole alignments dataframe instead of
    # alignments_one_primer_pair could be more efficient, but I am doing it here
    # for now to centralize the logic for identifying cross-amplification.
    filter(bitwAnd(sam_flag, 4) == 0 & primer_genome != off_target_genome) %>%
    group_by(off_target_genome)
  # We could abstract this quosure creation out to a function, but not for now.
  fprimer_mapped = quo(any(grepl(FORWARD, primer_id)))
  rprimer_mapped = quo(any(grepl(REVERSE, primer_id)))
  full_match_genomes = by_ref %>%
    filter(!!fprimer_mapped & !!rprimer_mapped) %>% ot_genomes()
  half_match_genomes = by_ref %>%
    filter(xor(!!fprimer_mapped, !!rprimer_mapped)) %>% ot_genomes()
  # Extract forward and reverse primer sequences.
  primer_seqs = c(T, F) %>% set_names(FORWARD, REVERSE) %>%
    map_chr(function(is_f) {alignments_one_primer_pair %>%
        filter(primer_forward == is_f) %>% slice(1) %>% pull(primer_seq)})
  # Summarize primer pair in one row (e.g. separate columns for forward and
  # reverse primers). Approach if I wanted to keep more primer columns:
  # primer_row = alignments_one_primer_pair %>%
  #   summarize(across(c(primer_genome), ~ first(.x)))
  primer_row = tibble_row()
  for (primer_dir in names(primer_seqs)) {
    primer_row = primer_row %>%
      mutate(!! glue("{primer_dir}_primer_seq") := primer_seqs[primer_dir])
  }
  # Alternative way to write:
  # fprimer_seq = filter(alignments_one_primer_pair, primer_forward==T)[1,"primer_seq"]
  # rprimer_seq = filter(alignments_one_primer_pair, primer_forward==F)[1,"primer_seq"]
  # primer_row = alignments_one_primer_pair %>%
  #   summarize(across(c(primer_id, primer_genome, primer_idx), ~ first(.x))) %>%
  #   mutate(forward_primer_seq = fprimer_seq, reverse_primer_seq = rprimer_seq)
  primer_row %>% 
    mutate(other_full_match_genomes = list(full_match_genomes),
           other_half_match_genomes = list(half_match_genomes),
           # num_matches calculation could be vectorized later across all primer
           # pairs. I'm keeping the logic in one place for now, but I might try
           # this if it speeds things up.
           num_matches = length(full_match_genomes) + length(half_match_genomes)/2)
}

most_spec_primer_pair = function(alignments_one_genome, ...) {
  print(alignments_one_genome$primer_genome[1])
  primer_pairs = alignments_one_genome %>% group_by(primer_idx) %>%
    group_modify(summarize_primer_off_targets) %>% ungroup()
  primer_pairs %>% filter(num_matches == min(num_matches))
}

# An approach where we removed genomes that were non-specific:
# primer_genomes = unique(alignments$primer_genome)
# primer_pairs = tibble()
# while (length(primer_genomes) > 0) {
#   primer_genome...
# }
# TODO: reduce runtime by reducing num primers.
primer_pairs = alignments %>%
  # filter(primer_genome %in% unique(alignments$primer_genome)[1:3]) %>%
  group_by(primer_genome) %>%
  group_modify(most_spec_primer_pair, .keep=T) %>% ungroup()
primer_pairs %>%
  # Convert list of genomes into comma-separated string for writing to file.
  mutate(across(tidyselect::starts_with("other_"),
                ~ map_chr(.x, ~ paste(., collapse=",")))) %>%
  write_tsv("primers/non_specific_primers.tsv")

#### OLD APPROACH

# Group alignments into the set of specific primers. Each primer (uniquely
# identified by primer_id) has 0 or more alignments to each genome (uniquely
# identified by off_target_genome). Only keep primers for all alignments are
# specific, i.e. no alignments (i) are to genomes besides the one that
# generated the primer and (ii) have the unmapped SAM flag (bit 0x4) unset.
primers = alignments %>% group_by(primer_id, primer_seq) %>%
  filter(!any(bitwAnd(sam_flag, 4) == 0 & primer_genome != off_target_genome)) %>%
  summarize(across(starts_with("primer_"), ~ first(.x)), .groups = "drop")

# Group primers into pairs and only keep the pairs where both primers passed the
# previous filter for specificity.
primer_pairs = primers %>% group_by(primer_genome, primer_idx) %>%
  # Keep primer groups which have the same index and are length 2 (i.e. pairs).
  filter(n() == 2) %>% ungroup(primer_idx) %>%
  # Keep one primer pair per genome (the one with minimum index).
  filter(primer_idx == min(primer_idx)) %>% ungroup() %>%
  select(-c(primer_idx))

# Append these specific primer pairs to the existing list.
if (file.exists(opt$output)) {
  all_primer_pairs = read_csv(opt$output) %>% rbind(primer_pairs)
} else {
  all_primer_pairs = primer_pairs
}
all_primer_pairs %>% arrange(primer_genome) %>% write_csv(opt$output)
