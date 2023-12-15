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

# Process parameters -----------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

option_list <- list(
  make_option(c("-g", "--gtf", type = "character"),
              help = "Path to user-provided annotation"),
  make_option(c("-q", "--quant", type = "character"),
              help = "Path to RNA-seq quantification file"),
  make_option(c("-c", "--counts", type = "character"),
              help = "Path to RNA-seq read counts for simulation"),
  make_option(c("-o", "--output", type = "character"),
              help = "Path to modified annotation output"),
  make_option(c("-t", "--tpm", type = "double"),
              default = 3,
              help = "Minimum TPM required to keep transcript")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser) # Load options from command line.


# Modify annotation ------------------------------------------------------------



# Load annotation to map genes to transcripts
gtf <- as_tibble(rtracklayer::import(opt$gtf))

message("gtf looks like:")

head(gtf)

# Load salmon quantification
salmon_quant <- fread(opt$quant)

message("salmon_quant looks like:")

head(salmon_quant)

# Filter out low abundance transcripts
quant_filter <- salmon_quant %>%
  filter(TPM >= opt$tpm | grepl(".I", Name)) %>%
  mutate(transcript_id = Name) %>%
  dplyr::select(-Name)

message("1st quant_filter looks like:")


head(quant_filter)

# Map transcripts to gene_ids
gene_to_tscript <- gtf %>%
  dplyr::filter(type == "transcript") %>%
  dplyr::mutate(location = paste0(seqnames, ":", start, "-", end)) %>%
  dplyr::select(gene_id, transcript_id, location) %>%
  dplyr::distinct()

message("genes_to_tscript looks like:")


head(gene_to_tscript)

quant_filter <- quant_filter %>%
  inner_join(gene_to_tscript,
             by = "transcript_id")

message("3rd to last quant_filter looks like:")


head(quant_filter)

# filter out genes with overestimated intronic content
weird_genes <- quant_filter %>%
  dplyr::group_by(gene_id) %>%
  dplyr::filter(TPM == max(TPM)) %>%
  dplyr::filter(grepl(".I", transcript_id)) %>%
  dplyr::select(gene_id) %>%
  dplyr::distinct()

message("weird_genes looks like:")

head(weird_genes)

quant_filter <- quant_filter %>%
  filter(!(gene_id %in% weird_genes$gene_id))

message("2nd to last quant_filter looks like:")


head(quant_filter)

# Check out intronic content
intronic <- quant_filter %>%
  filter(grepl(".I", transcript_id)) %>%
  arrange(TPM)

message("intronic looks like:")

head(intronic)

# Add reads mapping to filtered out transcripts to their respective gene transcript
filtered_out <- salmon_quant %>%
  filter(TPM < opt$tpm & !grepl(".I", Name) & TPM > 0) %>%
  mutate(transcript_id = Name) %>%
  dplyr::select(-Name) %>%
  inner_join(gene_to_tscript, 
             by = "transcript_id") %>%
  group_by(gene_id) %>%
  summarise(extra_intronic = sum(NumReads))

message("filtered_out looks like:")


head(filtered_out)

quant_filter <- quant_filter %>%
  inner_join(filtered_out,
             by = "gene_id") %>%
  mutate(NumReads = ifelse(grepl(".I", transcript_id), NumReads + extra_intronic,
                           NumReads))

message("Final quant_filter looks like:")


head(quant_filter)

# Create filtered annotation
transcript_keep <- quant_filter %>%
  filter(!grepl(".I", transcript_id)) %>%
  dplyr::select(transcript_id) %>%
  dplyr::distinct() %>%
  unlist() %>%
  unname()

gene_keep <- quant_filter %>%
  dplyr::select(gene_id) %>%
  dplyr::distinct() %>%
  unlist() %>%
  unname()


gtf_filter_t <- gtf %>%
  filter(transcript_id %in% transcript_keep)

gtf_filter_i <- gtf %>%
  filter(grepl(".I", transcript_id) & gene_id %in% gene_keep)

gtf_filter <- bind_rows(gtf_filter_t, gtf_filter_i) %>%
  filter(type %in% c("transcript", "exon"))


message("GTFs looks like:")

head(gtf_filter_t)
head(gtf_filter_i)
head(gtf_filter)

# Export filtered annotation
final_gr <- GRanges(seqnames = Rle(gtf_filter$seqnames),
                    ranges = IRanges(gtf_filter$start, end = gtf_filter$end),
                    strand = Rle(gtf_filter$strand))

mcols(final_gr) <- gtf_filter %>%
  dplyr::select(-seqnames, -start, -end, -strand, -width)

rtracklayer::export(final_gr, con = opt$output)



# Generate read counts to be used in simulation

normalized_reads <- quant_filter %>%
  ungroup() %>%
  mutate(norm_reads = NumReads/sum(NumReads))

write_csv(normalized_reads, file = opt$counts)
