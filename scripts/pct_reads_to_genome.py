# The metadata are at xa45-metadata-aligned-annotated.csv.
# Steps:
# 1. I loop through the table (maybe pandas?).
# 2. I need to call bowtie and concatenate.

# First read in data/xa45-metadata-aligned-annotated.csv, originally found at
# /n/groups/springer/amogh/analysis/xa45-metadata-aligned-annotated.csv.

import glob
import numpy as np
import os
import pandas as pd
import re
import subprocess

import shared_vars

DATA_DIR = "data"
POOLED_READS_DIR = "pooled_reads"
ASSEMBLIES_DIR = "unicycler_assemblies"
HOST_DATA_DIR = "abl19@o2.hms.harvard.edu:/n/groups/springer/amogh/analysis/xa45-redone"
BOWTIE_IDX = "bowtie_idx/strain_genome"

assembly_ids = [file.split(".")[0] for file in os.listdir(ASSEMBLIES_DIR)]
samples = pd.read_csv(os.path.join(DATA_DIR,
                                   "xa45-metadata-aligned-annotated.csv"))
samples = (samples[(samples["Source"] != "xanthobacter") &
                   (~(samples["Strain"].str.contains("SIMM")))]
           .rename(columns={"user_genome": "id"}))
samples = samples.assign(vial=samples["Strain"].str.extract(r"V(\d+).*").astype(int),
                         has_assembly=samples["id"].isin(assembly_ids))

# for each vial
# maybe a download / check files step. probably move download logic to
# another file. for now assume files may not exactly match. add sample ids in
# that other file too.
keep_cols = ["id", "vial", "Strain"]
pooled_samples = samples[samples["Source"] == "pool"][keep_cols]
lone_strain_samples = (samples.loc[samples["Source"] == "colony"]
    [keep_cols + ["has_assembly"]])
strains = pd.merge(pooled_samples, lone_strain_samples, how="inner",
                   on="vial", suffixes=("_pool", "_strain"))


# We added sample ids to sequence descriptions of assemblies: python3
# add_genomes_to_seq_descs.py unicycler_assemblies. Then we created bowtie
# index.
# export BOWTIE_IDX="bowtie_idx/strain_genome"
# bowtie2-build -f unicycler_assemblies/*.assembly.fasta $BOWTIE_IDX
# OK I might need to build indexes per vial.
def bowtie_build(vial_strains):
    # vial_strains, rather than pool_strains, because one vial can have
    # multiple pools for different time points.
    assemblies = [os.path.join(ASSEMBLIES_DIR, f"{id}.assembly.fasta")
                  for id in vial_strains["id_strain"].unique()]
    assemblies = [f for f in assemblies if os.path.isfile(f)]
    pattern = ",".join(assemblies)
    vial = vial_strains["vial"].iloc[0]
    # TODO: would be best to refactor f"{BOWTIE_IDX}_V{vial}" convention (
    #  also used in align).
    subprocess.run(["bowtie2-build", "-f", pattern, f"{BOWTIE_IDX}_V{vial}"])


strains.groupby("vial").apply(bowtie_build)


def align(pool_strains):
    pool_id = pool_strains["id_pool"].iloc[0]
    vial = pool_strains["vial"].iloc[0]
    # sam_file = os.path.join(POOLED_READS_DIR, f"{pool_id}.sam")
    sam_file_pattern = os.path.join(POOLED_READS_DIR, f"{pool_id}*.sam")
    sam_files = glob.glob(sam_file_pattern)
    # if not os.path.isfile(sam_file):
    # If you want to keep it simpler, we didn't need sam_file_pattern / glob
    # stuff, and we did not need to mv, and we didn't have logic in the else
    # statement.
    if not sam_files:
        result = subprocess.run(["bowtie2", "-q", "-a", "-x",
                                 f"{BOWTIE_IDX}_V{vial}", "-1",
                                 os.path.join(POOLED_READS_DIR, f"{pool_id}_R1.fastq.gz"),
                                 "-2",
                                 os.path.join(POOLED_READS_DIR, f"{pool_id}_R2.fastq.gz"),
                                 "-S", sam_file_pattern],
                                capture_output=True, text=True)
        # If alignment step is too slow, we might want -k instead of -a.
        # TODO: could make this more robust (e.g. requiring word 'reads').
        if pool_id == "LIB060784_GEN00271393_151_S151":
            breakpoint()
        num_pool_reads = int(re.match("\d+", result.stderr).group(0))
        sam_file = sam_file_pattern.replace("*.sam", f"_{num_pool_reads}R.sam")
        subprocess.run(["mv", sam_file_pattern, sam_file])
        return sam_file, num_pool_reads
    else:
        num_pool_reads = int(re.search("(\d+)R", sam_files[0]).group(1))
        return sam_files[0], num_pool_reads

# def filter2(strain_info, sam_df, num_pool_reads):
#     breakpoint()
#     strain_id = strain_info["id_strain"].iloc[0]
#     # Need to disambiguate multimapping?
#     num_strain_reads = sam_df["strain_genome_contig"].str.startswith(
#         strain_id).sum()
#     if num_strain_reads == 0 and strain_info["has_assembly"].iloc[0] is False:
#         num_strain_reads = np.nan()
#     return strain_info.assign(percent=num_strain_reads/num_pool_reads)


def filter(pool_strains, sam_file, num_pool_reads):
    # Alignments is the SAM file, a dataframe of alignments.
    aligns = pd.read_csv(sam_file, sep="\t", usecols=[0, 2], comment="@",
                             header=None)
    aligns.columns = ["read_id", "strain_genome_contig"]
    # Need to extract id_strain before grouping to avoid double-counting
    # reads that align to multiple contigs in the same genome.
    aligns = (aligns.loc[aligns["strain_genome_contig"] != "*"]
        .assign(id_strain = aligns["strain_genome_contig"].str.split(
            pat=shared_vars.GENOME_CONTIG_SEP, n=1, expand=True).iloc[:,0])
        .drop(columns=["strain_genome_contig"]))
    # Might need to dedupe within strain, which means extracting the genome.
    strain_reads = (aligns.drop_duplicates().groupby(
        "id_strain", as_index=False).size()
                    .rename(columns={"size": "num_reads"}))
    strain_reads["pct_reads"] = strain_reads["num_reads"] / num_pool_reads
    strain_reads.drop(columns=["num_reads"], inplace=True)
    # Optional step: join to pool_strains to assign NAs to strains with no
    # assemblies.
    strain_reads = pd.merge(strain_reads,
                            pool_strains[["id_strain","has_assembly"]],
                            how="inner", on="id_strain")
    strain_reads.loc[~strain_reads["has_assembly"], "pct_reads"] = np.nan
    return strain_reads


def calc_align_pcts(pool_strains):
    # Actually apply the following: 1. align to create the SAM file. 2. filter it
    # to just generate the dataframe of reads to each row. 3. also get the error
    # rate; potentially list that as a column (starting as same for the whole
    # group, could make it more granular). get as output from SAM stats. 4. delete
    # the SAM file.
    print(pool_strains["id_pool"].iloc[0])
    sam_file, num_pool_reads = align(pool_strains)
    align_pcts = filter(pool_strains, sam_file, num_pool_reads)
    # Remove if I don't want to lose data, but right now I need to save space.
    if any([id in sam_file for id in ["LIB060784_GEN00271346_104_S104",
                                      "LIB060784_GEN00271354_112_S112",
                                      "LIB060784_GEN00271362_120_S120",
                                      "LIB060784_GEN00271370_128_S128",
                                      "LIB060784_GEN00271378_136_S136",
                                      "LIB060784_GEN00271385_143_S143",
                                      "LIB060784_GEN00271386_144_S144",
                                      "LIB060784_GEN00271393_151_S151",
                                      "LIB060784_GEN00271394_152_S152",
                                      "LIB060784_GEN00271410_168_S168",
                                      "LIB060784_GEN00271417_175_S175"]]):
        subprocess.run(["rm", sam_file])
    return align_pcts


# Could parallelize instead to at least successful jobs run. Snakemake?
z = (strains.groupby("id_pool").apply(calc_align_pcts)
     .reset_index(level=0, names=['id_pool',_]))

# Align one pooled read to everything.
# Process immediately for all strains in the vial (could be another function):
# take the sam file, filter for the contig in question, and give out
# percentages. I.e. sam_file = f"{id_pool}.sam"; bowtie2 -1 f"{
# id_pool}_R1.fastq.gz" -2 f"{id_pool}_R2.fastq.gz" BOWTIE_IDX -o sam_file.
# Then read in sam_file with similar processing as in filter_primers.R. Let's
# check, but I'm guessing that the ref name will have the id with
# GENOME_CONTIG_SEP (colon). Split on that to get the id_strain. How to apply
# here? Could do a map step with a function (e.g. str.contains) 12 time; could
# also iterate through the rows and add to a counter dictionary. Maybe start
# with the latter and see what happens.
# results
# If I want to say NaN
# for a
# missing assembly,
# hm I could
# probably store that at an earlier step. Also have percentage that are
# missing (100% - others, possibly a bound if NA is there).
# Output after all of this should be a ~13-row dataframe.
# Each group will look like id_pool, id_strain

# MAYBE SKIP THIS BECAUSE OF DOWNLOAD SLOWNESS! Download pooled reads if not
# already
all_read_files = [f"{sample}_R{i}.fastq.gz" for i in [1, 2] for sample in
                  list(pooled_samples["id"])]
read_files2dl = set(all_read_files) - set(os.listdir(POOLED_READS_DIR))
host_paths = os.path.join(f"{HOST_DATA_DIR}", "concatenated",
                          f"{{{','.join(read_files2dl)}}}")
# O2 password needs to be input here; follow
# https://unix.stackexchange.com/questions/744852/scp-many-files-without-reentering-password-all-different-paths
# to avoid having to enter password with every file.
# subprocess.run(["scp", host_paths, POOLED_READS_DIR])
# gzip -d *.fastq.gz
# Could do the same with individual assemblies. Assembly files seem much
# smaller than read files. Assembly download commands:
# Remote:
# export LAB_DIR="/n/groups/springer"
# export NEW_DIR="$LAB_DIR/andrew/unicycler_assemblies"
# mkdir -p $NEW_DIR
# cd "$LAB_DIR/amogh/analysis/xa45-redone/unicycler"
# for genome in *; do
#   cp "$genome/assembly.fasta" "$NEW_DIR/$genome.assembly.fasta"
# done
# Locally:
# scp -r abl19@o2.hms.harvard.edu:/n/groups/springer/andrew/unicycler_assemblies .

# bowtie alignment step.

# Then the per-vial, per-strain step is still a function, but it could be
# much simpler? And we'll see if we can read the SAM file in as a dataframe (
# unless it's too big). If not, perhaps we can filter out reads not mapping
# to own vial; we should also only be mapping to own reads.
# Should probably distinguish NA (not having assembly) from 0 (had assembly
# but no reads mapped). Fuck, then I might redo the above step.

# run bowtie on the *reads* from that pooled thing to the *assemblies* from
# the individual ones. so we should have a function that takes the genome id
# from the pool and maps it to one strain's assembly, a function around that
# which does it 12 times per vial, and then is there somehow to apply this
# with pandas?

# ok wait more:
# Initially I would think to concatenate genomes per vial, or even everywhere
# (but then we get off-targets we don't care about), and name them
# appropriately, to reduce the number of bowtie indexes I need to create. The
# alternative would be to do them on each individual genome, which might be
# OK? It's actually the simpler approach but might take longer. I don't know
# what the bowtie creation runtime is. But if it takes 3 minutes per
# individual strain, then x 100 strains is 5 hours. If it's 10 minutes per
# strain, it's even worse. Whereas with 12 vials it's maybe 36-60 minutes. So
# yeah let's write the faster code. I could even consider doing it for all of
# them at once??
# This means I need to include the sample id in the sequence description of
# at least the strain genomes before concatenation.
# Then I create the bowtie indexes.
# Then OK on the group by side... first draw out assemblies belonging to each
# vial.

# and the resulting df would look like vial, within-vial strain, and percent
# of reads from the pool mapping to that strain, so still ~100 rows.
