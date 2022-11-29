## Chunyu Zhao 2021-12-13
#! /usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

library(tidyverse)
library(magrittr)


specify_decimal <- function(x, k) as.numeric(trimws(format(round(x, k), nsmall=k)))


my_species <- args[1]

datadir <- file.path("/mnt/cz/20220602_close_species", my_species)

my_ani_bin <- read_delim(file.path(datadir, "curr_fastani.tsv"), delim = "\t", show_col_types = F) %>%
    mutate(ani_bin = paste("ani.", ani_int, sep="")) %>%
    mutate(species_id = as.character(species_id)) %>%
    mutate(on_target_species = as.character(on_target_species))

ani_on_bins <- my_ani_bin %>% filter(species_id == on_target_species) %>% .$ani_bin

list_of_reads_summary <- list()
for (ani_on in ani_on_bins) {
  midas_dir <- file.path(datadir, ani_on, "6_midas_exp1")
  ani_offs <- list.files(midas_dir)
  for (ani_off in ani_offs) {
    summary_dir <- file.path(midas_dir, ani_off, "cov_20X/snps", "snps_summary.tsv")
    snps_summary <- read_delim(summary_dir, delim = "\t", show_col_types = F) %>%
      mutate(on_bin = ani_on, off_bin = ani_off, on_target_species = my_species)

    list_of_reads_summary[[paste(ani_on, ani_off, sep="-")]] <- snps_summary
  }
}

reads_summary <- bind_rows(list_of_reads_summary)
reads_summary %>% write.table(file.path(datadir, "exp1_snps_summary.tsv"), sep = "\t", quote=F, row.names = F)
