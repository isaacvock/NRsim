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
              default = 4.5,
              help = "Label time in hours. NOTE: the distribution of simulated 
              fraction news is not affected by this parameter. It only affects
              the conversion of simulated fraction new to kdeg."),
  make_option(c("-p", "--pkdeg", type = "numeric"),
              default = 0.5,
              help = "Proportion of genes for which isoform abundance 
              differences are kdeg driven"),
  make_option(c("-m", "--maxfn", type = "numeric"),
              default = 0.99,
              help = "Maximum simulated fraction new allowed."),
  make_option(c("-n", "--minfn", type = "numeric"),
              default = 0.01,
              help = "Minimum simulated fraction new allowed."),
  make_option(c("-f", "--avglfn", type = "numeric"),
              default = 0,
              help = "Average logit(fraction new) of dominant isoforms."),
  make_option(c("-s", "--sdlfn", type = "numeric"),
              default = 1,
              help = "Standard deviation of logit(fraction new) 
              of dominant isoforms."),
  make_option(c("-r", "--random", type = "numeric"),
              default = 0,
              help = "If set, then all isoforms are given independent kdegs
              drawn at random from a logit-Normal distribution. There will
              be no correlation between isoform abundance and isoform kdeg
              in this case. 0 = FALSE; 1 = TRUE")  
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser) # Load options from command line.


# Add fraction news to simulate ------------------------------------------------

if(!(opt$random %in% c(0, 1))){
  warning("--random (-r) should be 0 or 1! If not 0, treated as if it is 1!!")
}

### Simulation details and helper functions

# Label time
tl <- opt$labeltime


# What are min and max fns?
min_fn <- opt$minfn
max_fn <- opt$maxfn

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

if(opt$random == 0){

  # Draw a gene-wise fraction new (technicallly fraction new of most abundant isoform)
  fn_gene <- inv_logit(rnorm(ngenes, opt$avglfn, opt$sdlfn))
  fn_gene <- ifelse(fn_gene > max_fn, max_fn,
                    ifelse(fn_gene < min_fn, min_fn, fn_gene))

  # Adjust pkdeg_diff so that expected value of kdeg diffs = original pkdeg_diff
  nwithinbounds <- sum(!(fn_gene %in% c(max_fn, min_fn)))

  if (nwithinbounds > 0) {

    pkdeg_diff <- pkdeg_diff * (ngenes / nwithinbounds)

    if (pkdeg_diff > 1) {
      pkdeg_diff <- 1
    }

  }else {

    stop("There are no gene-wide fraction news not equal to the upper or lower bound
    on the fraction new. Did you make sure that `avglfn` falls between `minfn` and `maxfn?")

  }


  # Are differences in isoform abundances kdeg driven?
  kdeg_diff <- ifelse(fn_gene %in% c(max_fn, min_fn),
                        rbinom(ngenes, size = 1, prob = pkdeg_diff),
                        0)


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
    mutate(kdeg = ifelse(1 - exp(-kdeg * kdeg_factor * tl) > opt$maxfn,
                          runif(1, (kdeg - log(1 - opt$maxfn)/tl)/2,  -log(1 - opt$maxfn)/tl),
                          kdeg * kdeg_factor),
          fn = 1 - exp(-kdeg*tl)) %>%
    mutate(fn = ifelse(grepl(".I", transcript_id), 1, fn))


}else{

  fn_isoforms <- inv_logit(rnorm(nrow(normalized_reads), opt$avglfn, opt$sdlfn))
  fn_isoforms <- ifelse(fn_isoforms > max_fn, max_fn,
                    ifelse(fn_isoforms < min_fn, min_fn, fn_isoforms))

  nwithinbounds <- sum(!(fn_isoforms %in% c(max_fn, min_fn)))

  if (nwithinbounds == 0) {

    stop("There are no fraction news not equal to the upper or lower bound
    on the fraction new. Did you make sure that `avglfn` falls between `minfn` and `maxfn?")

  }


  normalized_reads$fn <- fn_isoforms

  normalized_reads <- normalized_reads %>%
    mutate(fn = ifelse(grepl(".I", transcript_id), 1, fn))


}

write_csv(normalized_reads, file = opt$output)