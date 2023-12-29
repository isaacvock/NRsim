#!/usr/bin/env Rscript
### PURPOSE OF THIS SCRIPT
## Filter out transcripts which are low coverage in a real dataset.
## Will also generate normalized read counts that will be used for simulation

# Load dependencies ------------------------------------------------------------

library(dplyr)
library(rtracklayer)
library(GenomicRanges)
library(optparse)
library(readr)
library(polyester)

# Process parameters -----------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

option_list <- list(
  make_option(c("-c", "--counts", type = "character"),
              help = "Path to csv with read counts for simulation"),
  make_option(c("-o", "--output", type = "character"),
              help = "Path to directory to output simulated data"),
  make_option(c("-f", "--fasta", type = "character"),
              help = "Path to transcriptome fasta for simulation"),
  make_option(c("-e", "--error_rate", type = "double"),
              default = 0.001,
              help = "Base calling error rate to simulate"),
  make_option(c("-n", "--nreps", type = "double"),
              default = 3,
              help = "Number of replicates to simulate"),
  make_option(c("-s", "--seed", type = "double"),
              default = -1,
              help = "Seed to set for simulation. If < 0, no seed is set"),
  make_option(c("-p", "--singleend"),
              action = "store_false",
              default = TRUE,
              help = "Simulate pre-mRNA"),
  make_option(c("-m", "--premRNA", type = "character"),
              default = "True",
              help = "Simulate single end data rather than paired-end"),
  make_option(c("-l", "--librarysize", type = "double"),
              default = 10000000,
              help = "Total number of reads to simulate")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser) # Load options from command line.


# Modify annotation ------------------------------------------------------------



# Load read counts
normalized_reads <- read_csv(opt$counts)


# (Mostly) remove those full gene transcripts if necessary
if(opt$premRNA == "False"){
  
  normalized_reads <- normalized_reads %>%
    mutate(norm_reads = ifelse(grepl(".I", transcript_id), 0.001/opt$librarysize, norm_reads))
  
}

# vector of reads
reads_per_transcript <- ceiling(normalized_reads$norm_reads * opt$librarysize)



# Fold changes
fold_changes <- matrix(1, nrow = nrow(normalized_reads),
                       ncol = 1)


### Simulate
if (opt$seed < 0) {

  simulate_experiment(fasta = opt$fasta,
                      outdir = opt$output,
                      num_reps = opt$nreps,
                      reads_per_transcript = reads_per_transcript,
                      fold_changes = fold_changes,
                      paired = opt$singleend,
                      error_rate = opt$error_rate,
                      strand_specific = TRUE)

}else {

  simulate_experiment(fasta = opt$fasta,
                      outdir = opt$output,
                      num_reps = opt$nreps,
                      reads_per_transcript = reads_per_transcript,
                      fold_changes = fold_changes,
                      error_rate = opt$error_rate,
                      paired = opt$singleend,
                      strand_specific = TRUE,
                      seed = opt$seed)
}
