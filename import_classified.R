library(tidyverse)

args       = commandArgs(trailingOnly=TRUE)
type       = args[1]
basefolder = args[2]

if(type == "all" | type == "singleunmapped"){
  collated_single = tibble()
  for (file in dir(paste0(basefolder, "results/classified/single_unmapped/"), pattern = "*classifiedreads*")){
    tmp = read_tsv(paste0(basefolder, "results/classified/single_unmapped/", file), col_names = c("classtype", "read_name", "taxID", "length", "k-mers info"))
    sample = gsub(".kraken2.classifiedreads.txt", "", file)
    collated_single = collated_single %>%
      bind_rows(
        bind_cols(
          sample = sample,
          tmp
        )
      )
  }
  saveRDS(collated_single, file = "classified_reads_collated_single.rds")
}

if(type == "all" | type == "bothunmapped"){
  collated_both = tibble()
  for (file in dir("results/classified/both_unmapped/", pattern = "*classifiedreads*")){
    tmp = read_tsv(paste0("results/classified/both_unmapped/", file), col_names = c("classtype", "read_name", "taxID", "length", "k-mers info"))
    sample = gsub(".kraken2.classifiedreads.txt", "", file)
    collated_both = collated_both %>%
      bind_rows(
        bind_cols(
          sample = sample,
          tmp
        )
      )
  }
  saveRDS(collated_both, file = "classified_reads_collated_both.rds")
}

