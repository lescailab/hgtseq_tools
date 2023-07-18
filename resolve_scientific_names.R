#!/usr/bin/env Rscript

args            = commandArgs(trailingOnly = TRUE)
input_data      = args[1]
output_data     = paste0(gsub(".rds","", input_data), "_resolved.rds")

writeLines("--- loading library")
library(taxonomizr)
library(tidyverse)

writeLines("######## R execution - arguments summary ######")
writeLines(paste0("## dataset to be parsed = ", input_data))
writeLines(paste0("## dataset to be written = ", output_data))

writeLines("## link local SQL database")
sqldb = "/home/lescailab/local_REFS/taxonomy_sql/accessionTaxa.sql"

writeLines("## setup parsing functions")
getGenusFromSpecies <- function(taxaid){
  taxonomy = getTaxonomy(taxaid, sqldb)
  genus = as_tibble(taxonomy)$genus[[1]]
  return(genus)
}

getLowestRankName <- function(taxaid, result="name"){
  taxonomy = as_tibble(getTaxonomy(taxaid, sqldb))
  if(!is.na(taxonomy$species[[1]])){
    scientific_name = taxonomy$species[[1]]
    rank_level = "species"
  } else if (!is.na(taxonomy$genus[[1]])){
    scientific_name = taxonomy$genus[[1]]
    rank_level = "genus"
  } else if (!is.na(taxonomy$family[[1]])){
    scientific_name = taxonomy$family[[1]]
    rank_level = "family"
  } else if (!is.na(taxonomy$order[[1]])){
    scientific_name = taxonomy$order[[1]]
    rank_level = "order"
  } else if (!is.na(taxonomy$class[[1]])){
    scientific_name = taxonomy$class[[1]]
    rank_level = "class"
  } else if (!is.na(taxonomy$phylum[[1]])){
    scientific_name = taxonomy$phylum[[1]]
    rank_level = "phylum"
  } else if (!is.na(taxonomy$superkingdom[[1]])){
    scientific_name = taxonomy$superkingdom[[1]]
    rank_level = "superkingdom"
  } else {
    scientific_name = NA
    rank_level = NA
  }
  if(result=="name"){
    return(scientific_name)
  }
  else {
    return(rank_level)
  }
}

writeLines("--- reading the data")
collated_reads = readRDS(input_data)

writeLines("--- resolving reads")
collated_reads_resolved = tibble(
  taxID = unique(collated_reads$taxID),
  sciname = unlist(lapply(unique(collated_reads$taxID), getLowestRankName)),
  rank = unlist(lapply(unique(collated_reads$taxID), getLowestRankName, "rank")),
  genus = unlist(lapply(unique(collated_reads$taxID), getGenusFromSpecies))
)



#####################################################
## creating original read-based resolved datasets ###
#####################################################


writeLines("--- joining resolved and original data")
collated_reads = collated_reads %>%
  left_join(collated_reads_resolved, by = "taxID")
writeLines("--- saving results")
saveRDS(collated_reads, file = output_data)

writeLines("--- JOB COMPLETED")