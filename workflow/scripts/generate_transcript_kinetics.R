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
              help = "Label time in hours"),
  make_option(c("-p", "--pkdeg", type = "numeric"),
              default = 0.5,
              help = "Proportion of genes for which isoform abundance 
              differences are kdeg driven"),
  make_option(c("-m", "--maxkdeg", type = "numeric"),
              default = 2.5,
              help = "Maximum simulated kdeg allowed (in hr-1)"),
  make_option(c("-n", "--minkdeg", type = "numeric"),
              default = 0.005,
              help = "Minimum simulated kdeg allowed")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser) # Load options from command line.


# Add fraction news to simulate ------------------------------------------------

### Simulation details and helper functions

# Label time
tl <- opt$labeltime


# What are min and max fns?
min_fn <- 1 - exp(-opt$minkdeg*tl)
max_fn <- 1 - exp(-opt$maxkdeg*tl)

# Probability that an isoforms kdeg is different from major isoform
pkdeg_diff <- opt$pkdeg

# Sigmoid function
inv_logit <- function(x){
  return(exp(x) / (1 + exp(x)))
}

# Simulated relative abundances
normalized_reads <- fread(opt$counts)

# Gene IDs
genes <- unique(normalized_reads$gene_id)

# Number of genes
ngenes <- length(genes)


### Simulate

# Are differences in isoform abundances kdeg driven?
kdeg_diff <- rbinom(ngenes, size = 1, prob = pkdeg_diff)

# Draw a gene-wise fraction new (technicallly fraction new of most abundant isoform)
fn_gene <- inv_logit(rnorm(ngenes, 0, 1))
fn_gene <- ifelse(fn_gene > max_fn, max_fn,
                  ifelse(fn_gene < min_fn, min_fn, fn_gene))

# Create a table to organize all simulation details
fns <- tibble(gene_id = genes,
              fn_gene = fn_gene,
              kdeg_diff = kdeg_diff)

# Calculate transcript isoform fraction news
  # If kdeg_diff is 1 (TRUE), then kdeg for an isoform
  # is (kdeg of dominant isoform)*[(dominant isoform TPM) / (isoform TPM)].
  # Else, all isoforms from a gene have same kdeg
normalized_reads <- normalized_reads %>%
  inner_join(fns, by = "gene_id") %>%
  mutate(kdeg = -log(1 - fn_gene)/tl) %>%
  group_by(gene_id) %>%
  mutate(kdeg_factor = ifelse(kdeg_diff == 1, max(TPM) / TPM, 1)) %>%
  rowwise() %>%
  mutate(kdeg = ifelse(kdeg*kdeg_factor > opt$maxkdeg, runif(1, kdeg, opt$maxkdeg), kdeg*kdeg_factor),
         fn = 1 - exp(-kdeg*tl)) %>%
  mutate(fn = ifelse(grepl(".I", transcript_id), 1, fn))


write_csv(normalized_reads, file = opt$output)