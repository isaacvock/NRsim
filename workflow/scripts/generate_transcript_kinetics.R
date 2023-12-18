#!/usr/bin/env Rscript
### PURPOSE OF THIS SCRIPT
## Filter out transcripts which are low coverage in a real dataset.
## Will also generate normalized read counts that will be used for simulation

# Load dependencies ------------------------------------------------------------

library(dplyr)
library(rtracklayer)
library(GenomicRanges)
library(optparse)
library(data.table)
library(readr)

# Process parameters -----------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

option_list <- list(
  make_option(c("-c", "--counts", type = "character"),
              help = "Path to RNA-seq read counts for simulation"),
  make_option(c("-o", "--output", type = "character"),
              help = "Path to modified annotation output"),
  make_option(c("-t", "--labeltime", type = "numeric"),
              default = 2,
              help = "Path to modified annotation output"),
  make_option(c("-p", "--pkdeg", type = "numeric"),
              default = 1,
              help = "Proportion of genes for which isoform abundance 
              differences are kdeg driven")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser) # Load options from command line.


# Add fraction news to simulate ------------------------------------------------

tl <- opt$labeltime
pkdeg_diff <- opt$pkdeg

inv_logit <- function(x){
  return(exp(x) / (1 + exp(x)))
}

normalized_reads <- fread(opt$counts)

genes <- unique(normalized_reads$gene_id)

kdeg_diff <- rbinom(ngenes, size = 1, prob = pkdeg_diff)

ngenes <- length(genes)

fn_gene <- inv_logit(rnorm(ngenes, 0, 1))

fns <- tibble(gene_id = genes,
              fn_gene = fn_gene,
              kdeg_diff = kdeg_diff)

normalized_reads <- normalized_reads %>%
  inner_join(fns, by = "gene_id") %>%
  mutate(kdeg = -log(1 - fn_gene)/tl) %>%
  group_by(gene_id) %>%
  mutate(kdeg_factor = ifelse(kdeg_diff == 1, TPM / max(TPM), 1)) %>%
  mutate(kdeg = kdeg*kdeg_factor,
         fn = 1 - exp(-kdeg*tl)) %>%
  mutate(fn = ifelse(grepl(".I", transcript_id), 1, fn))


write_csv(normalized_reads, file = opt$output)