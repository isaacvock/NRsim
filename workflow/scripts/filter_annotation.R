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
  make_option(c("-g", "--gtf", type = "character"),
              help = "Path to user-provided annotation"),
  make_option(c("-q", "--quant", type = "character"),
              help = "Path to RNA-seq quantification file"),
  make_option(c("-c", "--counts", type = "character"),
              help = "Path to RNA-seq read counts for simulation"),
  make_option(c("-o", "--output", type = "character"),
              help = "Path to modified annotation output"),
  make_option(c("-n", "--nointron", type = "character"),
              help = "Path to modified annotation output without introns"),
  make_option(c("-t", "--tpm", type = "double"),
              default = 3,
              help = "Minimum TPM required to keep transcript"),
  make_option(c("-i", "--introndiff", type = double),
              default = 1,
              help = "Minimum fraction difference between intronic coverage
              and least abundant transcript.")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser) # Load options from command line.


# Modify annotation ------------------------------------------------------------



# Load annotation to map genes to transcripts
gtf <- as_tibble(rtracklayer::import(opt$gtf))

# Load salmon quantification
salmon_quant <- fread(opt$quant)


# Filter out low abundance transcripts
quant_filter <- salmon_quant %>%
  filter(TPM >= opt$tpm | grepl(".I", Name)) %>%
  mutate(transcript_id = Name) %>%
  dplyr::select(-Name)

# Map transcripts to gene_ids
gene_to_tscript <- gtf %>%
  dplyr::filter(type == "transcript") %>%
  dplyr::mutate(location = paste0(seqnames, ":", start, "-", end)) %>%
  dplyr::select(gene_id, transcript_id, location) %>%
  dplyr::distinct()

quant_filter <- quant_filter %>%
  inner_join(gene_to_tscript,
             by = "transcript_id")

# filter out genes with overestimated intronic content
weird_genes <- quant_filter %>%
  dplyr::group_by(gene_id) %>%
  dplyr::filter(TPM == max(TPM)) %>%
  dplyr::filter(grepl(".I", transcript_id)) %>%
  dplyr::select(gene_id) %>%
  dplyr::distinct()

quant_filter <- quant_filter %>%
  filter(!(gene_id %in% weird_genes$gene_id))

message("quant_filter after weird_genes looks like:")

head(quant_filter)


# Check out intronic content
intronic <- quant_filter %>%
  filter(grepl(".I", transcript_id)) %>%
  arrange(TPM)


# Add reads mapping to filtered out transcripts to their respective gene transcript
filtered_out <- salmon_quant %>%
  filter(TPM < opt$tpm & !grepl(".I", Name) & TPM > 0) %>%
  mutate(transcript_id = Name) %>%
  dplyr::select(-Name) %>%
  inner_join(gene_to_tscript, 
             by = "transcript_id") %>%
  group_by(gene_id) %>%
  summarise(extra_intronic = sum(TPM))



quant_filter <- quant_filter %>%
  left_join(filtered_out,
             by = "gene_id") %>%
  mutate(extra_intronic = ifelse(is.na(extra_intronic), 0, extra_intronic)) %>%
  mutate(TPM = ifelse(grepl(".I", transcript_id),
                           TPM + extra_intronic + 0.1,
                           TPM))



# Create filtered annotation
transcript_keep <- quant_filter %>%
  dplyr::select(transcript_id) %>%
  dplyr::distinct() %>%
  unlist() %>%
  unname()


gtf_filter <- gtf %>%
  filter(transcript_id %in% transcript_keep) %>%
  filter(type %in% c("transcript", "exon"))


# Export filtered annotation
final_gr <- GRanges(seqnames = Rle(gtf_filter$seqnames),
                    ranges = IRanges(gtf_filter$start, end = gtf_filter$end),
                    strand = Rle(gtf_filter$strand))

mcols(final_gr) <- gtf_filter %>%
  dplyr::select(-seqnames, -start, -end, -strand, -width)


rtracklayer::export(final_gr, con = opt$output)


# Export intronless annotation
final_gr <- GRanges(seqnames = Rle(gtf_filter_t$seqnames),
                    ranges = IRanges(gtf_filter_t$start, end = gtf_filter_t$end),
                    strand = Rle(gtf_filter_t$strand))

mcols(final_gr) <- gtf_filter_t %>%
  dplyr::select(-seqnames, -start, -end, -strand, -width)


rtracklayer::export(final_gr, con = opt$nointron)


# Generate read counts to be used in simulation

normalized_reads <- quant_filter %>%
  ungroup() %>%
  group_by(gene_id) %>%
  mutate(TPM_adj = ifelse(grepl(".I", transcript_id), 
                          TPM,
                          TPM + sum(TPM[grepl(".I", transcript_id)])*opt$introndiff )) %>%
  ungroup() %>%
  mutate(NumReads_adj = NumReads*(TPM_adj/TPM)) %>%
  mutate(norm_reads = (NumReads_adj + 1) / sum(NumReads_adj + 1))

write_csv(normalized_reads, file = opt$counts)
